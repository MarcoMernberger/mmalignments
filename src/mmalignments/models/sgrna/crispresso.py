from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd  # type: ignore[import]
from pandas import DataFrame, Series  # type: ignore[import]

from mmalignments.models.artifacts import ArtifactSet, OutputSpec
from mmalignments.models.elements import (
    CallSpec,
    Element,
    Runnable,
    element,
    generate_element_key_name,
)
from mmalignments.models.parameters import (
    ParamSet,
)
from mmalignments.models.tables.frames import (  # type: ignore[import]
    Morph,
    Tables,
)
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)

from ..externals import External

logger = logging.getLogger(__name__)



class Crispresso(External):
    """Wrapper for the CRISPResso CLI (``CRISPResso`` binary)."""

    def __init__(
        self,
        name: str = "crispresso",
        primary_binary: str = "CRISPResso",
        version: str | None = None,
        source: str = "",
        parameters: Mapping[str, ParamSet] | ParamSet | None = None,
    ) -> None:
        resolved_parameters: Mapping[str, ParamSet] | ParamSet = (
            {}# if parameters is not None else _build_param_registry_as_mapping()
        )
        print(resolved_parameters)
        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=resolved_parameters,
        )
        # self.param_registry = _build_param_registry()

    def get_version(self, fallback: str | None = None) -> str | None:
        """Detect ``CRISPResso`` version from ``CRISPResso --version`` output."""
        if self._version:
            return self._version
        if not self.primary_binary or not self.ensure_binary():
            return fallback
        try:
            cp = subprocess.run(
                [self.primary_binary, "--version"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if cp.returncode == 0:
                out = (cp.stdout or cp.stderr or "").strip()
                if out:
                    m = re.search(r"(\d+\.\d+(?:\.\d+)?)", out)
                    if m:
                        return m.group(1)
                    return out.splitlines()[0]
        except Exception:
            pass
        return fallback

    def default_output_dir(self, root: str) -> Path:
        """Default output root for CRISPResso runs."""
        return Path("results") / "crispresso" / self.version_name / root

    def _resolve_path(
        self,
        data: Element | Path | str,
        *,
        preferred_keys: Sequence[str] = ("tsv", "count_table", "gene_summary", "path"),
    ) -> Path:
        """Resolve an input path from an Element or plain path."""
        if isinstance(data, Element):
            for key in preferred_keys:
                p = data.artifacts.get(key)
                if p:
                    return Path(p).absolute()
            if data.artifacts:
                return Path(next(iter(data.artifacts.values()))).absolute()
            raise ValueError("Input element has no artifacts to resolve a path")
        return Path(data).absolute()

    def _resolve_pres(self, data: Element | Path | str) -> tuple[Element, ...]:
        return (data,) if isinstance(data, Element) else ()


@element
def editframe(
    self,
    crispresso_folder: Path,
    file_suffix: str = ".Quantification_window_nucleotide_frequency_table.txt",
    *,
    root: str = "window.frequencies",
    tag: PartialElementTag | ElementTag | None = None,
    output_spec: OutputSpec | None = None,
) -> Element:
    """Generate a DataFrame of editing rates from CRISPResso output files."""
    tag = ElementTag(
        root=root,
        level=1,
        stage=Stage.ANALYSIS,
        method=Method.CRISPRESSO,
        state=State.PROCESSED,
        ext="tsv",
    ).merge(tag)
    key, name = generate_element_key_name(
        self.version_name,
        "editing_frame",
    )
    inputfile = crispresso_folder / file_suffix
    artifacts = {"tsv": inputfile} if inputfile.exists() else None

    artifacts, _ = ArtifactSet.generate_file_artifacts(
        tag=tag,
        infile=inputfile,
        spec=output_spec,
        index_column=None,
    )

    return Element(
        key=key,
        run=generate_editing_frame(crispresso_folder, file_suffix),
        tag=tag,
        determinants=None,
        inputs=(inputfile,),
        artifacts=artifacts,
        pres=None,
        name=name,
    )


def generate_editing_frame(
    crispresso_folder: Path,
    file_suffix: str = ".Quantification_window_nucleotide_frequency_table.txt",
) -> Runnable:

    def __call():
        files_meta = get_substitution_rates_files(crispresso_folder, file_suffix)
        df = load_substitution_files(files_meta)
        return df

    spec = CallSpec(
        path=("crispresso", "generate_editing_frame"),
        kwargs={
            "crispresso_folder": crispresso_folder,
            "file_suffix": file_suffix,
        },
    ).render()

    return Runnable(__call, display=spec)


def editingrate(
    self,
    source: Element,
    *,
    untreated: str = "None",
    root: str = "window.rates",
    tag: PartialElementTag | ElementTag | None = None,
    output_spec: OutputSpec | None = None,
) -> Element:

    tag = from_prior(
        source.tag,
        tag,
        root=root,
    )
    # key, name = generate_element_key_name(
    #     self.version_name,
    #     "editing_rate",
    # )
    # inputfile = source.artifacts.primary.resolve()
    # artifacts = ArtifactSet.generate_file_artifacts(
    #     tag=tag,
    #     infile=inputfile,
    #     spec=output_spec,
    #     index_column=None,
    # )[0]

    # return Element(
    #     key=key,
    #     run=calculate_editing_rate(untreated=untreated),
    #     tag=tag,
    #     determinants=None,
    #     inputs=(inputfile,),
    #     artifacts=artifacts,
    #     pres=None,
    #     name=name,
    # )
    tables = Tables()
    morphs = (calculate_editing_rate(untreated=untreated),)
    return tables.transform(source, morphs, tag=tag, output_spec=output_spec)


def calculate_editing_rate(untreated: str = "None") -> Morph:

    def call(data: DataFrame) -> DataFrame:
        data["editing_rate"] = editing_rate(data)
        concats = []
        for gene, df_gene in data.groupby("gene"):
            for mouse, df_mouse in df_gene.groupby("mouse"):
                df_mouse = df_mouse.sort_values(["treatment", "position"])
                df_ref = df_mouse[df_mouse["treatment"] == untreated].copy()
                if df_ref.empty:
                    raise ValueError(
                        f"No reference data for gene {gene}, mouse {mouse}"
                    )
                for treatment in df_mouse["treatment"].unique():
                    df_treat = df_mouse[
                        df_mouse["treatment"] == treatment
                    ].copy()  # noqa: E501
                    if treatment == untreated:
                        df_treat["delta_editing_rate"] = 0
                    else:
                        df_treat["delta_editing_rate"] = (
                            df_treat["editing_rate"].values
                            - df_ref["editing_rate"].values
                        )
                    concats.append(df_treat)
        df_merged = pd.concat(concats, ignore_index=True)
        return df_merged

    return Morph.from_callable(call)



def get_substitution_rates_files(
    crispresso_dir: Path = Path("results/crispresso"),
    file_suffix: str = ".Quantification_window_nucleotide_frequency_table.txt",
) -> pd.DataFrame:

    records = []

    for path in crispresso_dir.iterdir():

        if not path.is_dir():
            continue

        sample_name = path.name
        splits = sample_name.split("_")

        gene = splits[1]
        mouse = splits[2]
        treatment = f"{splits[3]}_{splits[4]}" if len(splits) >= 5 else splits[3]
        sub_file = path / f"CRISPResso_on_Report_{sample_name}" / f"{gene}{file_suffix}"

        records.append(
            {
                "sample": sample_name,
                "mouse": mouse,
                "gene": gene,
                "treatment": treatment,
                "path": sub_file,
            }
        )

    df = DataFrame.from_records(records)
    return df.set_index("sample")


def load_substitution_files(
    files_meta: DataFrame,
    sep: str = "\t",
) -> DataFrame:
    """
    Converts CRISPResso substitution matrices into long format.

    Input:
        files_meta DataFrame with:
            mouse, gene, treatment, path

    Output:
        sample, mouse, gene, treatment, position, base, count
    """

    concats = []
    ignore_bases = ["N", "-"]
    for sample, row in files_meta.iterrows():

        file = Path(row["path"])
        if not file.exists():
            raise FileNotFoundError(f"Missing file: {file}")

        # read the frame
        df = pd.read_csv(
            file,
            sep=sep,
            header=None,
        )
        # transpose and set columns
        df = df.T
        columns = ["base"] + df.iloc[0, 1:].to_list()
        df.columns = columns
        df.drop(index=0, inplace=True)
        for col in columns[1:]:
            df[col] = df[col].astype(float).astype(int)

        # calculate count of reference and non-reference bases
        reference = []
        non_reference = []
        for _, posrow in df.iterrows():
            ref = posrow["base"]
            count_non_ref = 0
            for base in columns[1:]:
                if base in ignore_bases:
                    continue
                if base == ref:
                    reference.append(posrow[base])
                else:
                    count_non_ref += posrow[base]
            non_reference.append(count_non_ref)
        df["reference"] = reference
        df["non_reference"] = non_reference

        # add meta columns
        meta = {
            "mouse": row["mouse"],
            "gene": row["gene"],
            "treatment": row["treatment"],
        }
        for key, value in meta.items():
            df[key] = value

        # add position and sample
        df["position"] = pd.RangeIndex(start=0, stop=len(df))
        df["sample"] = sample
        concats.append(df)
    df_merged = pd.concat(concats, ignore_index=True)
    df_merged = df_merged.sort_values(["sample", "position"])
    return df_merged


def editing_rate(
    data: DataFrame, reference: str = "reference", non_reference: str = "non_reference"
) -> Series:
    return data[non_reference] / (data[reference] + data[non_reference])
