from pathlib import Path
from typing import Literal

import pandas as pd
from mmalignments.models.elements import (
    CallSpec,
    Element,
    Runnable,
    TableElement,
    element,
    from_prior,
)
from mmalignments.models.parameters import Params
from mmalignments.models.tags import ElementTag, Method, PartialElementTag, State
from mmalignments.services.io import read_frame, write_frames


class Tables:

    @element
    def concat(
        self,
        *elements: Element,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        mode: Literal["tsv", "parquet", "both"] = "both",
        params: Params | None = None,
    ) -> TableElement:

        first_element = elements[0]
        root = root or first_element.tag.root
        tag = from_prior(
            first_element.tag,
            state=State.MERGED,
            method=Method.TABLES,
        )
        outdir = Path(outdir or first_element.file.parent)
        filename = filename or tag.default_output
        output_file = Path(outdir) / filename
        frame_paths = [element.file for element in elements]
        source = self.concat_tables(frame_paths, output_file, mode=mode, params=params)
        return TableElement(
            source=source,
            tag=tag,
            outdir=outdir,
            filename=filename,
            root=root,
            mode=mode,
            determinants=params.determinants() if params else None,
            pres=tuple(elements),
        )

    def concat_tables(
        self,
        paths: list[Path],
        output_path: Path,
        *,
        mode: Literal["tsv", "parquet", "both"] = "both",
        params: Params | None = None,
    ) -> Runnable:
        params = params or Params()

        def __run():
            # Concatenate the DataFrames
            frames = [read_frame(path) for path in paths]
            combined_df = pd.concat(frames, ignore_index=True)
            # Save the combined DataFrame to a file
            write_frames(combined_df, output_path, mode=mode, **params)

        spec = CallSpec(
            path=("tables", "concat_tables"),
            kwargs={
                "paths": paths,
                "output_path": output_path,
                "mode": mode,
                "params": params,
            },
        ).render()
        return Runnable(__run, display=spec)
