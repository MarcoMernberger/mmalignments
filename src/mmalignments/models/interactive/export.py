"""
export.py
──────────
Bundles the framework (spec.py, plots.py, orchestration.py, app.py) PLUS
all project-specific classes (Analysis subclass, custom PlotState subclasses)
and the real data into ONE standalone .py file.  The recipient runs:

    pip install -r requirements.txt
    panel serve standalone_app.py --show

No mmalignments import, no repo checkout, no path setup needed on their end.

How it works
────────────
1. The four framework modules are concatenated after stripping cross-module
   relative imports that would be broken in a flat file.
2. User-defined classes referenced by PlotConfig (type(config.analysis) and
   every class in config.extra_state / extra_state_bottom that is NOT already
   part of the framework) are extracted via inspect.getsource() and appended.
3. Every DataFrame in `data` is base64-encoded as Parquet and decoded at
   startup — one blob per key, so multi-source apps (raw + mageck, etc.) work.
4. PlotConfig is reconstructed as a readable Python literal so the recipient
   can open the file, read it, and tweak it if needed.
5. A .panel().servable() call is appended so `panel serve <file>.py` just works.

Limitations
───────────
- User-defined classes must be inspectable (defined in a .py file, not
  generated dynamically).  inspect.getsource() is used.
- Dependencies of user classes beyond pandas / param / plotly / panel must be
  listed in `extra_imports` and `extra_requirements`.
- Pass the ALREADY-TRANSFORMED data (post any initial_raw_data_transformation).
  The standalone app embeds it as-is, with no preprocessing step.
"""

from __future__ import annotations

import base64
import inspect
import io
import textwrap
from pathlib import Path
from typing import Type

import pandas as pd

from . import spec as _spec_mod
from . import plots as _plots_mod
from . import orchestration as _orch_mod
from . import app as _app_mod
from .spec import PlotConfig

_FRAMEWORK_MODULES = [_spec_mod, _plots_mod, _orch_mod, _app_mod]

# Short module names: "spec", "plots", "orchestration", "app"
_FW_SHORT_NAMES: set[str] = {m.__name__.split(".")[-1] for m in _FRAMEWORK_MODULES}

# Every class defined in the framework — used to distinguish user classes
_FW_CLASSES: set[type] = {
    obj for m in _FRAMEWORK_MODULES for _, obj in inspect.getmembers(m, inspect.isclass)
}


# ─────────────────────────────────────────────────────────────────────────────
# Framework source extraction
# ─────────────────────────────────────────────────────────────────────────────


def _strip_local_imports(src: str) -> str:
    """
    Remove from a module's source lines that must not appear in a flat file:
      - `from __future__ import annotations`  (only one allowed per file)
      - `pn.extension(...)` calls at module level  (deduped in template)
      - relative cross-module imports (`from .spec import ...`, etc.)
      - absolute package imports pointing at framework modules
    """
    kept: list[str] = []
    for line in src.splitlines():
        s = line.strip()
        if s == "from __future__ import annotations":
            continue
        if s.startswith("pn.extension("):
            continue
        # from .spec import ..., from .plots import ..., etc.
        if any(s.startswith(f"from .{name} import") for name in _FW_SHORT_NAMES):
            continue
        # from mmalignments.models.interactive.spec import ...
        if any(s.startswith(f"from {m.__name__} import") for m in _FRAMEWORK_MODULES):
            continue
        kept.append(line)
    return "\n".join(kept)


def _framework_source() -> str:
    sep = "\n\n# " + "─" * 79 + "\n\n"
    return sep.join(
        _strip_local_imports(inspect.getsource(m)) for m in _FRAMEWORK_MODULES
    )


# ─────────────────────────────────────────────────────────────────────────────
# User class detection & extraction
# ─────────────────────────────────────────────────────────────────────────────


def _is_user_class(cls: type) -> bool:
    return cls not in _FW_CLASSES


def _collect_user_classes(
    config: PlotConfig,
    extra_plot_classes: list[type],
) -> list[type]:
    """
    Returns, in stable deduplicated order:
      1. the Analysis subclass  (type(config.analysis))
      2. user-defined classes from config.extra_state / extra_state_bottom
      3. any explicitly passed extra_plot_classes
    Framework classes are silently skipped — they're already in _framework_source().
    """
    candidates = (
        [type(config.analysis)]
        + list(config.extra_state or [])
        + list(config.extra_state_bottom or [])
        + list(extra_plot_classes)
    )
    seen: set[type] = set()
    result: list[type] = []
    for cls in candidates:
        if cls not in seen and _is_user_class(cls):
            seen.add(cls)
            result.append(cls)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Data embedding
# ─────────────────────────────────────────────────────────────────────────────


def _embed_data(data: dict[str, pd.DataFrame]) -> str:
    """
    Renders the _DATA dict as Python source with each DataFrame embedded as
    a base64-encoded Parquet blob split into 76-char adjacent string literals
    (Python auto-concatenates them at compile time, so b64decode receives the
    full string while the file remains readable).
    """
    lines = ["_DATA: dict[str, pd.DataFrame] = {"]
    for key, df in data.items():
        buf = io.BytesIO()
        df.to_parquet(buf, index=False)
        b64 = base64.b64encode(buf.getvalue()).decode()
        chunks = [repr(b64[i : i + 76]) for i in range(0, len(b64), 76)]
        chunk_src = "\n        ".join(chunks)
        lines.append(f"    {key!r}: pd.read_parquet(io.BytesIO(base64.b64decode(")
        lines.append(f"        {chunk_src}")
        lines.append(f"    ))),")
    lines.append("}")
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# PlotConfig literal rendering
# ─────────────────────────────────────────────────────────────────────────────


def _render_config(config: PlotConfig) -> str:
    """
    Reconstructs a PlotConfig as readable Python source.  Class references
    in extra_state are written as bare names — those classes are defined
    earlier in the script (framework section or user-classes section).
    """
    analysis_name = type(config.analysis).__name__
    parts = ["config = PlotConfig(", f"    analysis={analysis_name}(),", "    layers=["]
    for layer in config.layers:
        parts += [
            "        PlotLayer(",
            f"            plot_type={layer.plot_type!r},",
            f"            alias={layer.alias!r},",
            f"            roles={layer.roles!r},",
            f"            source={layer.source!r},",
            f"            hover={layer.hover!r},",
            "        ),",
        ]
    parts.append("    ],")

    if config.labels:
        parts.append(f"    labels={config.labels!r},")
    if config.extra_state:
        names = ", ".join(cls.__name__ for cls in config.extra_state)
        parts.append(f"    extra_state=[{names}],")
    if config.extra_state_bottom:
        names = ", ".join(cls.__name__ for cls in config.extra_state_bottom)
        parts.append(f"    extra_state_bottom=[{names}],")
    if config.select_param is not None:
        parts.append(f"    select_param={config.select_param!r},")
    if config.shared_roles:
        parts.append(f"    shared_roles={config.shared_roles!r},")
    if config.sidebar_order is not None:
        parts.append(f"    sidebar_order={config.sidebar_order!r},")
    if config.panel_kwargs:
        parts.append(f"    panel_kwargs={config.panel_kwargs!r},")
    if config.explicit_color_maps:
        parts.append(f"    explicit_color_maps={config.explicit_color_maps!r},")
    parts.append(")")
    return "\n".join(parts)


# ─────────────────────────────────────────────────────────────────────────────
# Main export function
# ─────────────────────────────────────────────────────────────────────────────


def export_standalone_script(
    data: dict[str, pd.DataFrame],
    config: PlotConfig,
    outdir: Path | str = ".",
    stem: str = "standalone_app",
    extra_imports: list[str] = (),
    extra_requirements: list[str] = (),
    extra_plot_classes: list[Type] = (),
    title: str = "Interactive Analysis",
) -> Path:
    """
    Writes <outdir>/<stem>.py and <outdir>/requirements.txt.

    Parameters
    ----------
    data               : All DataFrames the app needs, keyed by source name
                         (same keys used as AnalysisResult.data, e.g.
                         {"raw": raw_df, "mageck": mageck_df}).
                         Pass the ALREADY-TRANSFORMED data.
    config             : The PlotConfig used to build the app.
    outdir             : Output directory (created if it does not exist).
    stem               : Base name for the output file (no .py extension).
    extra_imports      : Additional import lines for user-class dependencies
                         beyond pandas/param/plotly/panel.
    extra_requirements : Additional package names appended to requirements.txt.
    extra_plot_classes : Custom BasePlot subclasses used in config.layers that
                         are NOT in plots.PLOT_REGISTRY.
    title              : Browser-tab / window title for the app.

    Returns
    -------
    Path to the written .py file.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    user_classes = _collect_user_classes(config, list(extra_plot_classes))
    user_classes_src = "\n\n".join(
        textwrap.dedent(inspect.getsource(cls)) for cls in user_classes
    )

    # Build the script by simple sentinel replacement — avoids str.format()
    # clashing with { } in Python source, and string.Template clashing with $
    # in regex patterns or other string literals.
    script = _TEMPLATE
    script = script.replace("%%EXTRA_IMPORTS%%", "\n".join(extra_imports))
    script = script.replace("%%FRAMEWORK_SRC%%", _framework_source())
    script = script.replace("%%USER_CLASSES_SRC%%", user_classes_src)
    script = script.replace("%%DATA_SRC%%", _embed_data(data))
    script = script.replace("%%CONFIG_SRC%%", _render_config(config))
    script = script.replace("%%TITLE_REPR%%", repr(title))

    out_path = outdir / f"{stem}.py"
    out_path.write_text(script, encoding="utf-8")

    base_reqs = ["panel>=1.3", "plotly", "pandas", "param", "pyarrow"]
    all_reqs = base_reqs + [r for r in extra_requirements if r not in base_reqs]
    req_path = outdir / "requirements.txt"
    req_path.write_text("\n".join(all_reqs) + "\n", encoding="utf-8")

    kb = out_path.stat().st_size // 1024
    print(f"Wrote  {out_path}  ({kb} KB)")
    print(f"Wrote  {req_path}")
    print()
    print("Share both files.  Recipient runs:")
    print(f"  pip install -r requirements.txt")
    print(f"  panel serve {out_path.name} --show")
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# Script template  (sentinels %%NAME%% replaced via str.replace — safe against
# { } and $ characters that appear in the bundled Python source)
# ─────────────────────────────────────────────────────────────────────────────

_TEMPLATE = '''\
"""
Standalone interactive app — auto-generated, do not edit by hand.
Regenerate via the source project's export_standalone_script() function.

Usage
-----
    pip install -r requirements.txt
    panel serve standalone_app.py --show

    # remote / server access
    panel serve standalone_app.py --address 0.0.0.0 --port 5008 \\
        --allow-websocket-origin "*"
"""
from __future__ import annotations

import base64
import io
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
import panel as pn
import param
import plotly.express as px
import plotly.graph_objects as go
from param import Parameterized
%%EXTRA_IMPORTS%%

pn.extension("plotly")

# ── Framework (spec + plots + orchestration + app, bundled) ──────────────────
%%FRAMEWORK_SRC%%

# ── User-defined state classes and analysis ───────────────────────────────────
%%USER_CLASSES_SRC%%

# ── Embedded data (base64-encoded Parquet) ────────────────────────────────────
%%DATA_SRC%%

# ── App configuration ─────────────────────────────────────────────────────────
%%CONFIG_SRC%%

# ── Entry point ───────────────────────────────────────────────────────────────
pn.extension("plotly")
_app = build_app(_DATA, config, title=%%TITLE_REPR%%)
_app.panel().servable(title=%%TITLE_REPR%%)
'''
