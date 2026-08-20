"""interactive/apps.py
======================
DAG-integrated factory for Panel/Plotly interactive-app config files.

Usage inside a pipeline
-----------------------

    from mmalignments.models.interactive import Interactive

    interactive = Interactive(panel_app_script="/path/to/panel_app.py")

    # Rank plot: element depends on a long-format TableElement
    rank_el = interactive.app_element(
        long_table_element,
        analysis="count_rank",
        plot="rank",
        title="sgRNA Rank – HCT D4",
        port=5009,
        state_defaults={
            "analysis_topn": 200,
            "appearance_highlight_n": 30,
            "appearance_log_y": True,
        },
    )

    # Horizontal bar chart from the same data
    hbar_el = interactive.app_element(
        long_table_element,
        analysis="count_hbar",
        plot="hbar",
        title="Top 30 sgRNAs – HCT D4",
        port=5010,
    )

Each element's only artifact is a ``<name>.app.json`` config file.
No Panel server is started in the pipeline; the app is served separately:

    python panel_app.py --config path/to/element.app.json --serve

Or from within Python:

    from panel_app import build_app_from_config
    app = build_app_from_config("path/to/element.app.json")
    app.panel().show()
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from mmalignments.models.elements import (
    CallSpec,
    Element,
    Runnable,
    TableElement,
    element,
    generate_element_key_name,
)
from mmalignments.models.tags import (
    Method,
    State,
)
from mmalignments.models.overlay import (
    ElementTag,
    PartialElementTag,
    from_prior,
)


class Interactive:
    """Factory that turns a TableElement into an app-config Element.

    Parameters
    ----------
    panel_app_script:
        Absolute path to ``panel_app.py``.  Stored in the JSON config so
        users know how to serve the app.  Does not have to exist at build
        time.
    """

    def __init__(self, panel_app_script: str | Path | None = None) -> None:
        self.panel_app_script = Path(panel_app_script) if panel_app_script else None

    # ──────────────────────────────────────────────────────────────────────
    # Public API
    # ──────────────────────────────────────────────────────────────────────

    @element
    def app_element(
        self,
        data_element: TableElement,
        *,
        analysis: str = "count_rank",
        plot: str = "rank",
        title: str = "Interactive Analysis",
        port: int = 5009,
        address: str = "localhost",
        allow: str = "",
        figure_height: int = 560,
        state_defaults: dict[str, Any] | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: str | None = None,
    ) -> Element:
        """Create a DAG element whose artifact is a JSON app-config file.

        Parameters
        ----------
        data_element:
            The upstream ``TableElement`` that provides the data (long format).
        analysis:
            Analysis plugin registered in ``panel_app.analysis_registry``
            (e.g. ``"count_rank"`` or ``"count_hbar"``).
        plot:
            Plot-renderer plugin registered in ``panel_app.plot_registry``
            (e.g. ``"rank"`` or ``"hbar"``).
        title:
            Window title shown in the Panel app.
        port:
            Default port stored in the config (can be overridden at serve
            time with ``--port``).
        address:
            Default bind address stored in the config.
        allow:
            Extra allowed WebSocket origins (comma-separated).
        figure_height:
            Height in pixels of the Plotly pane.
        state_defaults:
            Dict of ``{param_name: value}`` written into the ``state``
            section of the JSON config, applied on app start.
        tag / outdir / filename:
            Override the auto-derived tag, output directory, or filename.
        """
        resolved_tag = from_prior(
            data_element.tag,
            tag,
            state=State.CONFIG,
            method=Method.CUSTOM,
            ext="json",
            param=f"{analysis}.{plot}",
        )

        out_dir = Path(outdir or data_element.file.parent)
        out_name = filename or resolved_tag.default_output
        output_file = out_dir / out_name

        run = self._write_config(
            data_path=data_element.file,
            output_path=output_file,
            analysis=analysis,
            plot=plot,
            title=title,
            port=port,
            address=address,
            allow=allow,
            figure_height=figure_height,
            state_defaults=state_defaults or {},
            panel_app_script=self.panel_app_script,
        )

        key, name = generate_element_key_name(
            resolved_tag, "Interactive", "app_element"
        )

        return Element(
            key,
            run,
            tag=resolved_tag,
            inputs=(data_element.file,),
            artifacts={"json": output_file},
            pres=(data_element,),
            name=name,
        )

    # ──────────────────────────────────────────────────────────────────────
    # Internal helpers
    # ──────────────────────────────────────────────────────────────────────

    @staticmethod
    def _write_config(
        data_path: Path,
        output_path: Path,
        analysis: str,
        plot: str,
        title: str,
        port: int,
        address: str,
        allow: str,
        figure_height: int,
        state_defaults: dict[str, Any],
        panel_app_script: Path | None,
    ) -> Runnable:
        """Return a Runnable that writes the JSON config when called."""

        def _run() -> None:
            config: dict[str, Any] = {
                "app": {
                    "data": str(data_path.resolve()),
                    "analysis": analysis,
                    "plot": plot,
                    "title": title,
                    "port": port,
                    "address": address,
                    "allow": allow,
                    "figure_height": figure_height,
                },
                "state": state_defaults,
            }
            if panel_app_script is not None:
                config["_meta"] = {
                    "serve_command": (
                        f"python {panel_app_script} "
                        f"--config {output_path} --serve --port {port}"
                    )
                }
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as fh:
                json.dump(config, fh, indent=2)

        spec = CallSpec(
            path=("interactive", "write_config"),
            kwargs={
                "data": str(data_path),
                "output": str(output_path),
                "analysis": analysis,
                "plot": plot,
            },
        ).render()

        return Runnable(_run, display=spec)
