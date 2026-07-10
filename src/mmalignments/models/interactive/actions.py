"Callable functions for action buttons in the interactive app."

from param import Parameterized  # type: ignore[import]
from typing import Callable
from pathlib import Path
from mmalignments.services.io import parents, write_frame


def action_reset(param_names: list) -> Callable[[Parameterized], None]:
    def __reset(interactive: Parameterized) -> None:
        print("Resetting PCAPlotState2")
        state = interactive.state
        for name in param_names:
            p = state.param[name]
            setattr(state, name, p.default)

    return __reset


def action_export_figure(fmt: str) -> Callable[[Parameterized], None]:

    def _export(interactive: Parameterized) -> None:
        if interactive._last_fig is None:
            return
        folder = getattr(interactive.state, "export_folder", "results")
        filename = getattr(interactive.state, "export_filename", "figure")
        path = Path(folder) / f"{filename}.{fmt}"
        path.parent.mkdir(parents=True, exist_ok=True)
        interactive._last_fig.write_image(str(path))
        print(f"Saved {path}")

    return _export


def action_export_data(
    key: str | None = None, format: str = "tsv"
) -> Callable[[Parameterized], None]:

    def _export(interactive: Parameterized) -> None:
        if interactive._last_result is None:
            return
        folder = getattr(interactive.state, "export_folder", "results")
        filename = getattr(interactive.state, "export_filename", "figure")
        if key:
            df = interactive._last_result.data[key]
            path = Path(folder) / f"{filename}.{key}.{format}"
            parents(path)
            write_frame(df, path)
            print(f"Saved table {key} to {path}.")
        else:
            for suffix, df in interactive._last_result.data.items():
                path = Path(folder) / f"{filename}.{suffix}.{format}"
                parents(path)
                write_frame(df, path)
                print(f"Saved table {suffix} to {path}.")

    return _export


def action_export_app() -> Callable[[Parameterized], None]:
    """
    Generate a self-contained standalone_app.py (+ requirements.txt) that
    runs without the mmalignments package.

    Uses a lazy relative import so that:
      - in the full package, export.py is imported at call time (avoids a
        circular dependency at module load: export.py imports app.py).
      - in a standalone script (which doesn't ship the export module), the
        ImportError / SystemError is caught silently and the button is
        harmlessly inert.
    """
    try:
        from . import export as _exp
    except (ImportError, SystemError):
        print("[export_app] Export module not available in this environment.")
        raise

    def __export(interactive: Parameterized) -> None:
        folder = getattr(interactive.state, "export_folder", "results")
        stem = getattr(interactive.state, "export_filename", "standalone_app")
        _exp.export_standalone_script(
            data=interactive.data,
            config=interactive.config,
            outdir=folder,
            stem=stem,
            title=interactive.title,
        )

    return __export
