from pathlib import Path

from matplotlib.figure import Figure  # type: ignore

from .io import parents

FIGURE_SUFFIXES = [".png", ".pdf", ".svg"]


def save_figure(
    fig: Figure, outfile: Path | str, suffixes: list[str] | None = None
) -> None:
    suffixes = suffixes or FIGURE_SUFFIXES
    base = Path(outfile)
    parents(base)
    for ext in suffixes:
        fig.savefig(str(base.with_suffix(ext)), bbox_inches="tight")
