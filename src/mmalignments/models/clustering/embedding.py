from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping, Sequence

from pandas import DataFrame
from sklearn.decomposition import PCA

from mmalignments.core.annotations import View
from mmalignments.model.elements import (
    CallSpec,
    Element,
    Runnable,
    element,
    generate_element_key_name,
)
from mmalignments.models.artifacts import ArtifactSet, OutputSpec
from mmalignments.models.externals import (
    ExternalRunConfig,
)
from mmalignments.models.parameters import Params
from mmalignments.models.tables.frames import (
    TableSource,
    as_table_source,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.services.io import write_frames, write_json

logger = logging.getLogger(__name__)


class SKlearn:
    """
    A class for SKlearn data analysis methods, including PCA and other
    dimensionality reduction techniques."""

    def __init__(self):
        self.default_output_spec = OutputSpec(
            ext="parquet", additional_extensions=["tsv"]
        )
        self.__suffix = "sklearn"

    @property
    def base_dir(self) -> Path:
        return Path("results") / self.__suffix

    @element
    def pca(
        self,
        source: Element,
        *,
        n_components: int | None = 5,
        whiten: bool = True,
        sample_columns: Sequence[str] | None = None,
        view: View | None = None,
        propagation: bool = False,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Principal Component Analysis (PCA)."""
        if sample_columns is None and view is None:
            raise ValueError(
                "Either sample_columns or view must be provided to select the counts columns."  # noqa: E501
            )
        tag = from_prior(
            source.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.PCA,
            state=State.NORMAL,
            ext="parquet",
            param="pca",
        )
        params = params or Params(
            n_components=n_components, whiten=whiten, random_state=24
        )
        infile = source.primary.resolve()
        artifacts, output = ArtifactSet.generate_file_artifacts(
            tag=tag,
            infile=infile,
            spec=output_spec or self.default_output_spec,
            # column_schema=new_schema,
        )
        artifacts["json"] = artifacts["parquet"].with_suffix(".json")
        determinants = params.determinants()
        key, name = generate_element_key_name(tag, self.version_name, "vst")
        runner = self.call_vst(
            source=as_table_source(source),
            output=output,
            parameter=params.to_dict(),
            sample_columns=sample_columns,
            view=view,
            propagation=source.primary.load if propagation else None,
            cfg=cfg,
        )
        return Element(
            key=key,
            run=runner,
            tag=tag,
            determinants=determinants + (str(view),),
            inputs=(infile,),
            artifacts=artifacts,
            pres=(source,),
            name=name,
        )

    def call_pca(
        self,
        source: TableSource,
        output: Path | list[Path],
        *,
        sample_columns: Sequence[str] | None = None,
        view: View | None = None,
        n_components: int | None = 5,
        whiten: bool = False,
        parameter: Mapping[str, Any] | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Runnable:

        outputs = output if isinstance(output, list) else [output]
        cfg = cfg or ExternalRunConfig(threads=1)

        def __run():
            count_df = (
                source.view(view)[sample_columns]
                if sample_columns
                else source.view(view)
            )
            model = PCA(**parameter).fit(count_df)
            explained_variance = model.explained_variance_ratio_
            df_pca = DataFrame(
                model.transform(count_df),
                index=count_df.index,
                columns=[f"PC{i+1}" for i in range(model.n_components_)],
            )
            write_json(
                {"explained_variance_ratio": explained_variance.tolist()},
                outputs[0].with_suffix(".json"),
            )
            write_frames(df_pca, outputs)

        callspec = CallSpec(
            paths=("SKlearn", "call_pca"),
            kwargs=parameter.to_dict() if parameter else {},
        )
        return Runnable(run=__run, display=callspec.render())
