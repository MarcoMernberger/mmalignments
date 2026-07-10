from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping, Sequence

from pandas import DataFrame  # type: ignore[import]
from sklearn.decomposition import PCA  # type: ignore[import]

from mmalignments.core.annotations import View
from mmalignments.models.artifacts import ArtifactSet, OutputSpec, TableArtifact
from mmalignments.models.elements import (
    CallSpec,
    Element,
    Runnable,
    element,
    generate_element_key_name,
)
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
from mmalignments.services.io import write_frames

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
        self.name = "SKlearn"

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
        axis: int = 0,
        sample_columns: Sequence[str] | None = None,
        view: View | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
        params: Params | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Element:
        """Principal Component Analysis (PCA)."""
        # if sample_columns is None and view is None:
        #     raise ValueError(
        #         "Either sample_columns or view must be provided to select the counts columns."  # noqa: E501
        #     )
        tag = from_prior(
            source.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.SKLEARN,
            state=State.PCA,
            ext="parquet",
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
        artifacts_ev = {
            f"explained_variance{o.suffix}": TableArtifact(
                o.with_suffix(f".explained_variance{o.suffix}")
            )
            for o in output
        }
        artifacts = artifacts.with_extras(artifacts_ev)
        output = artifacts.output_files()
        determinants = params.determinants() + (
            str(n_components),
            str(whiten),
            str(axis),
        )
        key, name = generate_element_key_name(tag, self.name, subroutine="pca")
        runner = self.call_pca(
            source=as_table_source(source),
            output=output,
            sample_columns=sample_columns,
            view=view,
            parameter=params.to_dict(),
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
        output: Mapping[str, Path],
        *,
        sample_columns: Sequence[str] | None = None,
        axis: int = 0,
        view: View | None = None,
        parameter: dict[str, Any] | None = None,
        cfg: ExternalRunConfig | None = None,
    ) -> Runnable:

        cfg = cfg or ExternalRunConfig(threads=1)
        params = parameter or {"n_components": 5, "whiten": False}

        def __run():
            count_df = (
                source.view(view)[sample_columns]
                if sample_columns
                else source.view(view)
            )
            if axis == 0:
                count_df = count_df.T
            logger.info(f"Running PCA on {count_df.shape} data matrix.")
            model = PCA(**params).fit(count_df)
            explained_variance = model.explained_variance_ratio_
            df_pca = DataFrame(
                model.transform(count_df),
                index=count_df.index,
                columns=[f"PC{i+1}" for i in range(model.n_components_)],
            )
            explained_variance_df = DataFrame(
                {
                    "Explained_Variance": explained_variance.tolist(),
                    "PC": [f"PC{i+1}" for i in range(model.n_components_)],
                }
            )
            pca_output = [
                output
                for ext, output in output.items()
                if not ext.startswith("explained_variance.")
            ]
            variance_output = [
                output
                for ext, output in output.items()
                if ext.startswith("explained_variance.")
            ]
            write_frames(explained_variance_df, variance_output)
            write_frames(df_pca, pca_output)

        callspec = CallSpec(
            path=("SKlearn", "call_pca"),
            kwargs=params,
        )
        return Runnable(__run, display=callspec.render())
