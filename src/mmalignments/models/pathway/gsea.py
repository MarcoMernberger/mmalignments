"""Wrappers for gseapy-based pathway analysis (GSEA, prerank, ORA/enrichr)
and MSigDB gene-set management.

Three analysis modes are covered:

``GseaPy.gsea``
    Classical phenotype-based GSEA — requires an expression matrix (genes ×
    samples) **and** class labels.  Use when you have ≥3 samples per group.

``GseaPy.prerank``
    Weighted pre-ranked GSEA — the standard follow-up after DESeq2/edgeR.
    Input: a two-column TSV ``gene / stat`` sorted descending by stat.

``GseaPy.enrichr``
    Over-representation analysis (ORA) against Enrichr libraries — input is a
    plain gene list (significant genes only).

MSigDB gene-set management is handled by the separate :class:`MSigDB` helper,
which downloads GMT files and optionally translates HGNC gene symbols to
Ensembl IDs using a caller-supplied mapping table.

Design notes
------------
- All three analyses are pure-Python (no subprocess), so they use ``Runnable``
  directly, not ``@subroutine``.
- Every method follows the ``@element`` pattern: it constructs the
  ``Runnable`` closure, declares its inputs/artifacts/pres, and returns an
  ``Element`` that the :class:`~mmalignments.models.executor.Executor` can
  schedule, cache, and skip.
- ``MSigDB`` is a plain helper class (not an ``External`` subclass) because
  it makes network requests rather than calling a local binary.
"""

from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Sequence

import gseapy as gp  # type: ignore[import]
import pandas as pd  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]

from mmalignments.core.annotations import ColumnTag, View
from mmalignments.models.artifacts import ArtifactSet, OutputSpec
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FileElement,
    Runnable,
    element,
    generate_element_key_name,
)
from mmalignments.models.parameters import Params
from mmalignments.models.tags import (
    ElementTag,
    Method,
    PartialElementTag,
    Stage,
    State,
    from_prior,
)
from mmalignments.models.tags import PartialElementTag as PTag
from mmalignments.services.dependencies import depends
from mmalignments.services.io import (
    concat_files,
    exists,
    parents,
    read_frame,
    write_frames,
)

logger = logging.getLogger(__name__)


class Species(str, Enum):
    Homo_sapiens = "Homo_sapiens"
    Mus_musculus = "Mus_musculus"


###############################################################################
# GseaPy wrapper
###############################################################################


class Gseapy:
    """Wrapper for the ``gseapy`` Python library.

    Each method follows the ``@element`` pattern: it returns an
    :class:`~mmalignments.models.elements.Element` whose signature captures
    all analysis parameters, so re-runs are skipped when nothing changed.

    Typical usage::

        gp = GseaPy()

        # (1) Classical GSEA
        gsea_el = gp.gsea(
            expression=expr_element,   # TableElement: genes × samples
            cls=["WT","WT","KO","KO"],
            gene_sets=hallmark_el,     # FileElement (.gmt) from MSigDB
        )

        # (2) Pre-ranked GSEA (most common after DESeq2)
        prerank_el = gp.prerank(
            ranking=ranking_el,        # TableElement: gene / stat
            gene_sets=hallmark_el,
        )

        # (3) ORA / Enrichr
        ora_el = gp.enrichr(
            gene_list=sig_genes_el,    # FileElement: one gene per line
            gene_sets=hallmark_el,
        )
    """

    # ------------------------------------------------------------------
    # Defaults
    # ------------------------------------------------------------------

    DEFAULT_MIN_SIZE = 15
    DEFAULT_MAX_SIZE = 500
    DEFAULT_PERMUTATIONS = 1000
    DEFAULT_SEED = 42

    def default_outdir(self, method: str, root: str) -> Path:
        return Path("results") / "gsea" / method / root

    def default_params(self) -> Params:
        return Params(
            permutation_type="phenotype",
            permutation_num=self.DEFAULT_PERMUTATIONS,
            min_size=self.DEFAULT_MIN_SIZE,
            max_size=self.DEFAULT_MAX_SIZE,
            weight=1.0,
            gene_col="gene",
            stat_col="stat",
            seed=self.DEFAULT_SEED,
        )

    ####################################################################
    # GCT
    ####################################################################

    @element
    def gct(
        self,
        expression: Element,
        annotation: Element,
        *,
        name_col: str = "name [META]",
        view: View | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        output_spec: OutputSpec | None = None,
    ) -> Element:
        """Create a gct file.

        Parameters
        ----------
        expression : Element
            Expression matrix (genes × samples).  Must have one gene-name
            column (``gene_col``) and one sample column per sample.
        annotation : Element
            Annotation matrix (genes × annotation columns).  Must have one
            gene-name column (``gene_col``) and one annotation column per gene.
        name_col : str
            Name of the gene-name column in *annotation* (default
            ``"name [META]"``).
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        output_spec : OutputSpec | None
            Optional output specification.  Defaults to ``parquet`` + ``tsv``.

        Returns
        -------
        Element
            GCT file element.
        """
        pres = (expression, annotation)
        tag = from_prior(
            expression.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.GSEA,
            state=State.ENRICHMENT,
            param="expression",
        )
        output_spec = OutputSpec(
            tag.default_output,
            self.default_outdir("expr", expression.name),
            ext="parquet",
        )
        artifacts = ArtifactSet.generate_from_outspec(
            tag=tag,
            infile=expression.primary.resolve(),
            spec=output_spec or self.default_output_spec(),
            default_dir=self.default_outdir("expr", expression.root),
        )

        @depends(_run_gsea)
        def _run():
            expression_frame = expression.primary.view()
            index_name = expression_frame.index.name
            # rename = ColumnTag.generate_new_tag(
            #     old_tag=ColumnTag(...),  # replace with actual old_tag
            #     view=View(...),  # replace with actual view
            #     append=True,  # or False, depending on your logic
            # )
            expression_frame = expression_frame.reset_index()
            expression_frame = expression_frame.rename(
                columns={index_name: "Gene"}
            )  # rename the Ensembl IDs to Gene for now
            expression_frame.insert(
                0, "NAME", annotation.primary.view()[name_col].values
            )
            columns = ColumnTag.select_from_view(
                expression_frame.columns, view=view
            )  # select samples that are actually relevant for the comparison
            expression_frame = expression_frame[
                ["NAME", "Gene"] + list(columns)
            ]  # the final expression frame
            outfiles = list(artifacts.output_files().values())
            write_frames(expression_frame, outfiles)

        runner = Runnable(_run, display=f"gsea({expression.name})")
        key, name = generate_element_key_name(tag, "gseapy", subroutine="gsea")

        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            pres=tuple(pres),
            name=name,
        )

    # ------------------------------------------------------------------
    # gsea — classical phenotype-permutation GSEA
    # ------------------------------------------------------------------

    @element
    def gsea(
        self,
        gct: Element,
        name: str,
        classes: Sequence[str] | FileElement,
        gene_sets: FileElement | str | Sequence[str],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        params: Params | None = None,
    ) -> Element:
        """Classical GSEA with phenotype permutations.

        Parameters
        ----------
        gct : Element
            GCT file element containing the expression matrix (genes × samples).
            Must have one gene-name column (``gene_col``) and one sample column
            per sample.
        name : str
            Name of the analysis (used in output filenames).
        classes : Sequence[str] | FileElement
            Class labels — either a Python list of strings matching the
            sample columns, or a :class:`FileElement` pointing to a
            ``.cls`` file.
        gene_sets : FileElement | str | Sequence[str]
            Gene-set collection — a FileElement pointing to a GMT file, a
            single Enrichr library name, or a list of library names.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outdir : Path | str | None
            Output directory.  Defaults to ``results/gsea/gsea/<root>``.
        permutation_type : str
            ``"phenotype"`` (recommended) or ``"gene_set"``.
        permutation_num : int
            Number of permutations (default 1000).
        min_size / max_size : int
            Gene-set size filter.
        weight : float
            Enrichment score weighting (default 1.0 = weighted GSEA).
        gene_col : str
            Name of the gene-identifier column in *expression* (default
            ``"gene"``).
        seed : int
            Random seed for reproducibility.

        Returns
        -------
        Element
            ``tsv`` → summary results TSV; ``parquet`` → Parquet version.
        """
        params = params or self.default_params()
        pres: list[Element] = [gct]
        gene_sets_arg, pres = _resolve_gene_sets(gene_sets, pres)
        classes_arg, pres = _resolve_cls(classes, pres)
        params = self.default_params().update(params)
        tag = from_prior(
            gct.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.GSEA,
            state=State.ENRICHMENT,
            ext="tsv",
            param="gsea",
        )

        outdir_path = Path(
            outdir or self.default_outdir("gsea", gct.root)
        )  # noqa: E501

        # determinants = (
        #     str(permutation_type),
        #     str(permutation_num),
        #     str(min_size),
        #     str(max_size),
        #     str(weight),
        #     str(gene_col),
        #     str(seed),
        # )

        determinants = (str(classes_arg), str(gene_sets_arg))
        determinants = determinants + params.determinants()

        @depends(_run_gsea)
        def _run():
            gct_file = gct.primary.resolve()
            return _run_gsea(
                name=name,
                expression_file=gct_file,
                gene_sets_arg=gene_sets_arg,
                cls_arg=classes_arg,
                outdir=outdir_path,
                permutation_type=params.permutation_type,
                permutation_num=params.permutation_num,
                min_size=params.min_size,
                max_size=params.max_size,
                weight=params.weight,
                gene_col=params.gene_col,
                seed=params.seed,
            )

        runner = Runnable(_run, display=f"gsea({gct.name})")
        key, name = generate_element_key_name(tag, "gseapy", subroutine="gsea")
        artifacts = ArtifactSet(
            outdir_path / "gseapy.phenotype.gsea.report.csv",
            primary_name="csv",
            cls=outdir_path / "group.cls",
            rnk=outdir_path / "gsea_data.rnk",
            gmt=outdir_path / "gene_sets.gmt",
        )
        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            inputs=(),
            determinants=determinants,
            pres=tuple(pres),
            name=name,
        )

    # ------------------------------------------------------------------
    # prerank — weighted pre-ranked GSEA
    # ------------------------------------------------------------------

    @element
    def prerank(
        self,
        ranking: Element,
        gene_sets: FileElement | str | Sequence[str],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        gene_col: str = "id",
        stat_col: str = "stat",
        permutation_num: int = DEFAULT_PERMUTATIONS,
        min_size: int = DEFAULT_MIN_SIZE,
        max_size: int = DEFAULT_MAX_SIZE,
        weight: float = 1.0,
        seed: int = DEFAULT_SEED,
    ) -> Element:
        """Pre-ranked GSEA — the standard workflow after DESeq2/edgeR.

        Parameters
        ----------
        ranking : Element
            Element whose ``tsv`` artifact is a two-column TSV with columns
            *gene_col* and *stat_col*, sorted descending by *stat_col*.
            Typically produced by a DESeq2 or edgeR element.
        gene_sets : FileElement | str | Sequence[str]
            Gene-set collection.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outdir : Path | str | None
            Output directory.
        gene_col : str
            Name of the gene column (default ``"gene"``).
        stat_col : str
            Name of the ranking statistic column (default ``"stat"``).
        permutation_num : int
            Number of permutations.
        min_size / max_size : int
            Gene-set size filter.
        weight : float
            Enrichment score weighting exponent.
        seed : int
            Random seed.

        Returns
        -------
        Element
            ``tsv`` with GSEA result columns (es, nes, pval,
            fdr, …).
        """
        pres: list[Element] = [ranking]
        gene_sets_arg, pres = _resolve_gene_sets(gene_sets, pres)

        tag = from_prior(
            ranking.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.GSEA,
            state=State.ENRICHMENT,
            ext="tsv",
            root=ranking.tag.root,
            param="prerank",
        )

        output_dir = Path(
            outdir or self.default_outdir("prerank", ranking.root)
        )  # noqa: E501
        outfile = output_dir / "check_real_outout"
        determinants = (
            str(gene_col),
            str(stat_col),
            str(permutation_num),
            str(min_size),
            str(max_size),
            str(weight),
            str(seed),
        )

        @depends(_run_prerank)
        def _run():
            return _run_prerank(
                ranking_tsv=Path(ranking.tsv),
                gene_sets_arg=gene_sets_arg,
                gene_col=gene_col,
                stat_col=stat_col,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                weight=weight,
                outdir=output_dir,
                seed=seed,
            )

        runner = Runnable(_run, display=f"prerank({ranking.name})")
        artifacts = ArtifactSet(outfile)
        key, name = generate_element_key_name(
            tag, "gseapy", subroutine="prerank"
        )  # noqa: E501
        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            determinants=determinants,
            inputs=(ranking.tsv,),
            pres=tuple(pres),
            name=name,
        )

    # ------------------------------------------------------------------
    # enrichr — ORA against Enrichr/MSigDB libraries
    # ------------------------------------------------------------------

    @element
    def enrichr(
        self,
        gene_list: Element | FileElement,
        gene_sets: FileElement | str | Sequence[str],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
        gene_col: str | None = "id",
        organism: str = "Human",
        cutoff: float = 0.05,
    ) -> Element:
        """Over-representation analysis (ORA) via Enrichr.

        Parameters
        ----------
        gene_list : Element | FileElement
            Source of the gene list.  Two modes:

            * **FileElement** pointing to a plain text file (one gene per line)
              — used as-is.
            * Any **Element** with a ``tsv`` artifact — *gene_col* selects the
              column to use as the gene list (required in this case).

        gene_sets : FileElement | str | Sequence[str]
            Gene-set collection — FileElement (GMT), a single library name,
            or a list of names (e.g. ``["Hallmark_2020", "KEGG_2021_Human"]``).
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outdir : Path | str | None
            Output directory.
        gene_col : str | None
            Column to use as the gene list when *gene_list* is a
            ``TableElement``.  Not needed when *gene_list* is a plain
            ``FileElement``.
        organism : str
            ``"Human"`` or ``"Mouse"`` (default ``"Human"``).
        cutoff : float
            Adjusted p-value cutoff written into the result file (does not
            filter the output table).

        Returns
        -------
        TableElement
            ``tsv`` / ``parquet`` with ORA result columns (Term, Overlap,
            P-value, Adjusted P-value, Genes, …).
        """
        pres: list[Element] = [gene_list]
        gene_sets_arg, pres = _resolve_gene_sets(gene_sets, pres)

        tag = from_prior(
            gene_list.tag,
            tag,
            stage=Stage.ANALYSIS,
            method=Method.GSEA,
            state=State.ENRICHMENT,
            ext="tsv",
            param="enrichr",
        )

        outdir_path = Path(
            outdir or self.default_outdir("enrichr", gene_list.root)
        )  # noqa: E501

        determinants = (
            str(organism),
            str(cutoff),
            str(gene_col or ""),
        )

        @depends(_run_enrichr)
        def _run():
            return _run_enrichr(
                gene_list_element=gene_list,
                gene_sets_arg=gene_sets_arg,
                outdir=outdir_path,
                gene_col=gene_col,
                organism=organism,
                cutoff=cutoff,
            )

        run = Runnable(_run, display=f"enrichr({gene_list.name})")

        inputs = (Path(gene_list.tsv),) if hasattr(gene_list, "tsv") else ()
        key, name = generate_element_key_name(
            tag, "gseapy", subroutine="enrichr"
        )  # noqa: E501
        artifacts = ArtifactSet(outdir_path / "check_real_outout")
        return Element(
            key=key,
            run=run,
            tag=tag,
            artifacts=artifacts,
            inputs=inputs,
            determinants=determinants,
            pres=tuple(pres),
            name=name,
        )


###############################################################################
# MSigDB gene-set manager
###############################################################################


class MSigDB:
    """Download and manage MSigDB gene sets as pipeline :class:`Element`.

    MSigDB organises gene sets into *collections*:

    ======  ================================================
    ``H``   Hallmark (50 coherent gene sets)
    ``C1``  Positional (chromosomal location)
    ``C2``  Curated (KEGG, Reactome, BioCarta, …)
    ``C3``  Regulatory (motif / miRNA targets)
    ``C4``  Computational
    ``C5``  GO (BP, MF, CC)
    ``C6``  Oncogenic signatures
    ``C7``  Immunologic signatures
    ``C8``  Cell-type signatures (single-cell)
    ======  ================================================

    Typical usage::

        msigdb = MSigDB()

        # Download Hallmark as HGNC gene symbols
        hallmark = msigdb.download("H", organism="Human")

        # Translate to Ensembl IDs (requires a gene-name → Ensembl mapping)
        hallmark_ens = msigdb.to_ensembl(hallmark, id_map=ensembl_map_element)

        # Feed into GSEA
        prerank_el = gp.prerank(ranking=ranking_el, gene_sets=hallmark_ens)

    The ``id_map`` element must point to a TSV with at least two columns:
    ``gene_name`` (HGNC symbol) and ``ensembl_id``.  You can produce this from
    an :class:`~mmalignments.models.genome.EnsemblGenome` or by downloading
    a BioMart table.
    """

    # # gseapy has a built-in GMT URL registry; we also support the msigdb.org
    # # REST API for full control.
    # _MSIGDB_GMT_URL = (
    #     "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs"
    # )

    def __init__(
        self,
        cache_dir: Path | str | None = None,
        species: Species = Species.Homo_sapiens,
        version: str = "2026.1",
    ) -> None:
        """
        Parameters
        ----------
        cache_dir : Path | str | None
            Directory where downloaded GMT files are stored.
            Defaults to ``cache/msigdb``.
        """
        self.cache_dir = Path(cache_dir or "cache/msigdb")
        self._version = version
        self._species = species
        self._species_short = "Hs" if species == "Homo_sapiens" else "Mm"

    @property
    def url(self) -> str:
        return f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{self._version}.{self._species_short}"  # noqa: E501

    @property
    def collections_files(self) -> dict[str, str]:
        if self._species == "Homo_sapiens":
            return {
                "H": f"h.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C1": f"c1.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C2": f"c2.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C3": f"c3.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C4": f"c4.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C5": f"c5.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C6": f"c6.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C7": f"c7.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C8": f"c8.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
                "C9": f"c9.all.v{self._version}.{self._species_short}.symbols.gmt",  # noqa: E501
            }
        else:
            # Mouse collections
            return {
                "MH": f"mh.all.v{self._version}.Mm.symbols.gmt",
                "M1": f"m1.all.v{self._version}.Mm.symbols.gmt",
                "M2": f"m2.all.v{self._version}.Mm.symbols.gmt",
                "M3": f"m3.all.v{self._version}.Mm.symbols.gmt",
                "M4": f"m4.all.v{self._version}.Mm.symbols.gmt",
                "M5": f"m5.all.v{self._version}.Mm.symbols.gmt",
                "M6": f"m6.all.v{self._version}.Mm.symbols.gmt",
                "M7": f"m7.all.v{self._version}.Mm.symbols.gmt",
            }

    def get_version(self) -> str:
        """
        Returns
        -------
        str
            Version string (e.g. ``"2026.1"``).
        """
        return self._version

    # ------------------------------------------------------------------
    # download
    # ------------------------------------------------------------------

    @element
    def download(
        self,
        collection: str,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
    ) -> Element:
        """Download a MSigDB gene-set collection as a GMT file.

        If the file already exists (signature matches) the download is
        skipped like any other pipeline element.

        Parameters
        ----------
        collection : str
            MSigDB collection identifier, e.g. ``"H"``, ``"C2"``, ``"C5"``.
            Append ``"_mm"`` for mouse (``"H_mm"``, ``"C2_mm"``).
        tag : PartialElementTag | ElementTag | None
            Optional tag override.

        Returns
        -------
        Element
            Points to the local GMT file.
        """

        collection_key = collection.upper()
        if collection_key not in self.collections_files:
            raise ValueError(
                f"Unknown MSigDB collection: {collection!r}.  "
                f"Known collections: {sorted(self.collections_files)}"
            )
        output_dir = Path(outdir or self.cache_dir)
        filename = self.collections_files[collection_key]
        gmt_path = output_dir / filename

        tag = (
            PTag(
                root=f"{collection_key.lower()}",
                level=0,
                stage=Stage.INPUT,
                method=Method.MSIGDB,
                state=State.DOWNLOAD,
                ext="gmt",
            )
            .merge(tag)
            .resolve()
        )

        @depends(_download_gmt)
        def _run():
            _download_gmt(
                collection_key=collection_key,
                filename=filename,
                gmt_path=gmt_path,
                base_url=self.url,
            )

        runner = Runnable(_run, display=f"msigdb.download({collection_key})")
        artifacts = ArtifactSet(gmt_path)
        key, name = generate_element_key_name(
            tag, "MSigDB", subroutine="download"
        )  # noqa: E501
        return Element(
            key,
            runner,
            tag=tag,
            artifacts=artifacts,
            inputs=(),
            determinants=(collection_key,),
            pres=(),
        )

    @element
    def merge(
        self,
        gmts: Iterable[Element],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
    ) -> Element:
        """Merge multiple MSigDB GMT files into a single GMT file.

        If the file already exists (signature matches) the merge is
        skipped like any other pipeline element.

        Parameters
        ----------
        gmts : Iterable[Element]
            List of GMT files to merge.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.
        outdir : Path | str | None
            Output directory for the merged GMT file.
        filename : Path | str | None
            Filename for the merged GMT file.

        Returns
        -------
        Element
            Points to the local GMT file.
        """
        first_element = next(iter(gmts), None)
        if first_element is None:
            raise ValueError("No GMT files provided for merging.")
        root = f"msigdb_{['_'.join(element.root for element in gmts)]}"
        tag = from_prior(
            first_element.tag,
            tag,
            stage=Stage.INPUT,
            state=State.MERGED,
            root=root,
            ext="gmt",
        )
        output_dir = Path(outdir or self.cache_dir)
        filename = filename or tag.default_output
        gmt_path = output_dir / filename

        @depends(_download_gmt, concat_files)
        def _run():
            concat_files(
                gmt_path, *(element.primary.resolve() for element in gmts)
            )  # noqa: E501

        spec = CallSpec(
            path=("MsigDB", "merge"),
            args=(gmts,),
        ).render()
        runner = Runnable(_run, display=spec)
        key, name = generate_element_key_name(tag, "MSigDB", subroutine="merge")
        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=ArtifactSet(gmt_path),
            inputs=tuple(element.primary.resolve() for element in gmts),
            pres=tuple(gmts),
            name=name,
        )

    @element
    def collections(
        self,
        collections: Iterable[str],
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
    ) -> Element:
        """
        Download and merge multiple MSigDB collections into a single GMT file.

        Parameters
        ----------
        collections : Iterable[str]
            List of MSigDB collections to download and merge.
        tag : PartialElementTag | ElementTag | None, optional
            Optional tag override, by default None
        outdir : Path | str | None, optional
            Output directory for the merged GMT file, by default None
        filename : Path | str | None, optional
            Filename for the merged GMT file, by default None

        Returns
        -------
        Element
            Points to the local GMT file.
        """
        gmts = [
            self.download(collection, outdir=outdir)
            for collection in collections  # noqa: E501
        ]  # noqa: E501

        return self.merge(gmts, tag=tag, outdir=outdir, filename=filename)

    def default_HS(self) -> Element:
        """
            Points to the default Homo_sapiens MSigDB Hallmark gene set
            collection we mostly use in our analyses.

        Returns
        -------
        Element
            An Element of H, C2, C6, C7 geneset in a single gmt.
        """
        return self.collections(
            ["H", "C2", "C6", "C7"],
            tag=PTag(
                stage=Stage.INPUT,
                method=Method.MSIGDB,
                state=State.MERGED,
                ext="gmt",
            ),
        )

    def default_Mm(self) -> Element:
        """
            Points to the default Mus_musculus MSigDB Hallmark gene set
            collection we mostly use in our analyses.

        Returns
        -------
        Element
            An Element of MH, M2, M6, M7 geneset in a single gmt.
        """
        return self.collections(
            ["MH", "M2", "M6", "M7"],
            tag=PTag(
                stage=Stage.INPUT,
                method=Method.MSIGDB,
                state=State.MERGED,
                ext="gmt",
            ),
        )

    @element
    def select(
        self,
        gmt: Element,
        genesets: Iterable[str],
        *,
        root: str | None = None,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
    ) -> Element:
        """
        Select and merge single MSigDB gene sets from a GMT file into a new GMT
        file.

        Gene set ids from MSigSB must be prefixed with the MSigSB collection
        name e.g. "H.HALLMARK_ADIPOGENESIS, C2.KEGG_P53_SIGNALLING". The input
        GMT file can be obtained from the download method or any custom GMT
        file.

        Parameters
        ----------
        gmt : Element
            Input GMT file containing the gene sets to select from.
        genesets : Iterable[str]
            List of MSigDB gene sets to select and merge.
        tag : PartialElementTag | ElementTag | None, optional
            Optional tag override, by default None
        outdir : Path | str | None, optional
            Output directory for the merged GMT file, by default None
        filename : Path | str | None, optional
            Filename for the merged GMT file, by default None

        Returns
        -------
        Element
            Points to the local GMT file.
        """
        if not genesets:
            raise ValueError("No gene sets provided for selection.")

        tag = from_prior(
            gmt.tag,
            tag,
            stage=Stage.INPUT,
            state=State.DOWNLOAD,
            root=root,
            ext="gmt",
        )
        output_dir = Path(outdir or self.cache_dir)
        filename = filename or tag.default_output
        gmt_path = output_dir / filename

        @depends(_select_gene_sets)
        def _run():
            _select_gene_sets(
                input_gmt=gmt.primary.resolve(),
                selected_genesets=genesets,
                output_gmt=gmt_path,
            )

        spec = CallSpec(
            path=("MsigDB", "select"),
            args=(gmt, genesets),
        ).render()
        runner = Runnable(_run, display=spec)
        key, name = generate_element_key_name(
            tag, "MSigDB", subroutine="select"
        )  # noqa: E501
        artifacts = ArtifactSet(gmt_path)
        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            inputs=(gmt.primary.resolve(),),
            determinants=tuple(genesets),
            pres=(gmt,),
            name=name,
        )

    @element
    def custom(
        self,
        genesets_msig: Iterable[str] | None,
        genesets_custom: FileElement | Iterable[FileElement] | None,
        *,
        tag: PartialElementTag | ElementTag | None = None,
        outdir: Path | str | None = None,
        filename: Path | str | None = None,
    ) -> Element:
        """
        Download and merge single MSigDB gene sets and self-compiled custom gene
        sets into a single GMT file.

        Gene set ids from MSigSB must be prefixed with the MSigSB collection
        name e.g. "H.HALLMARK_ADIPOGENESIS, C2.KEGG_P53_SIGNALLING". If
        self-defined custom gene sets are to be used, we assume that instead of
        ids to be retrieved, ready FileElements representing the input files are
        supplied.

        Parameters
        ----------
        genesets_msig : Iterable[str]
            List of MSigDB gene sets to download and merge.
        genesets_custom : Iterable[Element]
            List of custom gene sets to merge.
        tag : PartialElementTag | ElementTag | None, optional
            Optional tag override, by default None
        outdir : Path | str | None, optional
            Output directory for the merged GMT file, by default None
        filename : Path | str | None, optional
            Filename for the merged GMT file, by default None

        Returns
        -------
        Element
            Points to the local GMT file.
        """
        if genesets_msig is None and genesets_custom is None:
            raise ValueError(
                "At least one of genesets_msig or genesets_custom must be provided."  # noqa: E501
            )

        genesets_from_msig = {}
        if genesets_msig is not None:
            for geneset_id in genesets_msig:
                if "." not in geneset_id:
                    raise ValueError(
                        f"Invalid MSigDB gene set id: {geneset_id!r}.  Must be prefixed with collection name, e.g. 'H.HALLMARK_ADIPOGENESIS'."  # noqa: E501
                    )
                splits = geneset_id.split(".")
                collection, gene_set = splits[0], splits[1]
                genesets_from_msig.setdefault(collection, []).append(gene_set)
        # download the collections needed
        gmts = []
        for collection, gene_sets in genesets_from_msig.items():
            gmt = self.download(collection, outdir=outdir)
            selected = self.select(gmt, gene_sets, outdir=outdir)
            gmts.append(selected)

        # if we have custom sets as well, add them to the merge stack
        if genesets_custom:
            gmts.extend(
                [genesets_custom]
                if isinstance(genesets_custom, Element)
                else genesets_custom
            )
        # retrieve the relevant sets for each collection and merge them
        return self.merge(gmts, tag=tag, outdir=outdir, filename=filename)

    # ------------------------------------------------------------------
    # to_ensembl
    # ------------------------------------------------------------------

    @element
    def to_ensembl(
        self,
        gmt: Element,
        id_map: Element,
        *,
        gene_name_col: str = "gene_name",
        ensembl_col: str = "ensembl_id",
        tag: PartialElementTag | ElementTag | None = None,
    ) -> Element:
        """Translate HGNC gene symbols in a GMT file to Ensembl IDs.

        Parameters
        ----------
        gmt : Element
            Element whose ``gmt`` artifact points to the source GMT file
            (gene symbols).
            Typically the output of :meth:`download_as_element`.
        id_map : Element
            Element whose ``tsv`` artifact is a TSV with at minimum two
            columns:

            * ``gene_name_col`` (default ``"gene_name"``) — HGNC symbol
            * ``ensembl_col``   (default ``"ensembl_id"``) — Ensembl gene ID

            You can produce this table from an
            :class:`~mmalignments.models.genome.EnsemblGenome` or any BioMart
            export.
        gene_name_col : str
            Column name for gene symbols in *id_map*.
        ensembl_col : str
            Column name for Ensembl IDs in *id_map*.
        tag : PartialElementTag | ElementTag | None
            Optional tag override.

        Returns
        -------
        Element
            ``gmt`` artifact → GMT file with Ensembl IDs substituted for gene
            symbols.  Gene sets where no Ensembl ID is found are retained with
            the original symbol (a warning is logged).
        """
        src_gmt = Path(gmt.gmt)
        out_gmt = src_gmt.with_name(src_gmt.stem + ".ensembl.gmt")

        derived_tag = from_prior(
            gmt.tag,
            tag,
            stage=Stage.INPUT,
            method=Method.MSIGDB,
            state=State.TRANSLATE,
            ext="gmt",
        )

        determinants = (gene_name_col, ensembl_col)
        key, name = _build_key_name(derived_tag, "to_ensembl")

        @depends(_translate_gmt_to_ensembl)
        def _run():
            _translate_gmt_to_ensembl(
                src_gmt=src_gmt,
                id_map_tsv=Path(id_map.tsv),
                out_gmt=out_gmt,
                gene_name_col=gene_name_col,
                ensembl_col=ensembl_col,
            )

        runner = Runnable(_run, display=f"msigdb.to_ensembl({gmt.name})")

        return Element(
            key,
            runner,
            tag=derived_tag,
            artifacts={"gmt": out_gmt},
            determinants=determinants,
            inputs=(src_gmt, Path(id_map.tsv)),
            pres=(gmt, id_map),
            name=name,
        )


###############################################################################
# Private runtime helpers  (referenced by @depends for fingerprinting)
###############################################################################


def _run_gsea(
    *,
    name: str,
    expression_file: Path,
    gene_sets_arg,
    cls_arg,
    outdir: Path,
    permutation_type: str,
    permutation_num: int,
    min_size: int,
    max_size: int,
    weight: float,
    gene_col: str,
    seed: int,
) -> None:

    expr_df = read_frame(expression_file)
    # gseapy expects genes as the index
    if gene_col in expr_df.columns:
        expr_df = expr_df.set_index(gene_col)

    print(gene_sets_arg)
    print(cls_arg)
    res = gp.gsea(
        data=expr_df,
        gene_sets=gene_sets_arg,
        cls=cls_arg,
        permutation_type=permutation_type,
        permutation_num=permutation_num,
        min_size=min_size,
        max_size=max_size,
        weight=weight,
        outdir=str(outdir),
        name=name,
        seed=seed,
        verbose=True,
    )

    df = res.res2d
    logger.info("GSEA done — %d terms", len(df))
    return df


def _run_prerank(
    *,
    ranking_tsv: Path,
    gene_sets_arg,
    gene_col: str,
    stat_col: str,
    permutation_num: int,
    min_size: int,
    max_size: int,
    weight: float,
    outdir: Path,
    seed: int,
) -> None:

    rnk = pd.read_csv(ranking_tsv, sep="\t")
    rnk = rnk.sort_values(stat_col, ascending=False)

    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets_arg,
        permutation_num=permutation_num,
        min_size=min_size,
        max_size=max_size,
        weight=weight,
        outdir=str(outdir),
        seed=seed,
        verbose=True,
    )
    print(outdir)
    df = res.res2d
    logger.info("prerank done — %d terms in %s", df.shape[0], outdir)
    return df


def _run_enrichr(
    *,
    gene_list_element: Element,
    gene_sets_arg,
    outdir: Path,
    gene_col: str | None,
    organism: str,
    cutoff: float,
) -> DataFrame:

    exists(outdir)

    # Resolve gene list
    if gene_col is not None and hasattr(gene_list_element, "tsv"):
        df_src = pd.read_csv(Path(gene_list_element.tsv), sep="\t")
        genes = df_src[gene_col].dropna().tolist()
    elif hasattr(gene_list_element, "path"):
        # plain FileElement — one gene per line
        genes = Path(gene_list_element.path).read_text().splitlines()
        genes = [g.strip() for g in genes if g.strip()]
    else:
        raise ValueError(
            "enrichr: could not resolve gene list from element "
            f"{gene_list_element.name!r}.  Pass gene_col if the element has a tsv."  # noqa: E501
        )

    res = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets_arg,
        organism=organism,
        cutoff=cutoff,
        outdir=str(outdir),
        verbose=True,
    )

    df = res.results
    logger.info("enrichr done — %d terms in %s", len(df), outdir)
    return df
    # df.to_csv(out_tsv, sep="\t", index=False)
    # df.to_parquet(out_parquet, index=False)


def _download_gmt(
    *,
    collection_key: str,
    filename: str,
    gmt_path: Path,
    base_url: str,
) -> None:
    """Download a MSigDB GMT file from the Broad Institute server."""
    import urllib.request

    parents(gmt_path)
    if gmt_path.exists():
        logger.info("GMT already present: %s", gmt_path)
        return

    url = f"{base_url}/{filename}"

    logger.info("Downloading %s → %s", url, gmt_path)
    urllib.request.urlretrieve(url, gmt_path)
    logger.info("Downloaded %s (%d bytes)", gmt_path, gmt_path.stat().st_size)


def _translate_gmt_to_ensembl(
    *,
    src_gmt: Path,
    id_map_tsv: Path,
    out_gmt: Path,
    gene_name_col: str,
    ensembl_col: str,
) -> None:
    """Rewrite a GMT file replacing HGNC symbols with Ensembl IDs.

    Lines whose genes cannot be mapped are kept as-is (with a warning).
    The translation table is read from *id_map_tsv*.
    """
    parents(out_gmt)

    id_map = (
        pd.read_csv(id_map_tsv, sep="\t")
        .dropna(subset=[gene_name_col, ensembl_col])
        .set_index(gene_name_col)[ensembl_col]
        .to_dict()
    )
    logger.info("Loaded %d gene-name → Ensembl mappings", len(id_map))

    total_genes = 0
    unmapped = 0
    out_lines: list[str] = []

    for raw_line in src_gmt.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip():
            continue
        parts = raw_line.split("\t")
        name, url, *genes = parts
        translated: list[str] = []
        for g in genes:
            eid = id_map.get(g)
            if eid:
                translated.append(eid)
            else:
                translated.append(g)  # keep original if not found
                unmapped += 1
            total_genes += 1
        out_lines.append("\t".join([name, url] + translated))

    if unmapped:
        logger.warning(
            "to_ensembl: %d/%d gene entries had no Ensembl mapping "
            "(kept as symbols)",
            unmapped,
            total_genes,
        )

    out_gmt.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    logger.info(
        "Wrote translated GMT → %s (%d gene sets)", out_gmt, len(out_lines)
    )  # noqa: E501


###############################################################################
# Internal utilities
###############################################################################


def _resolve_gene_sets(
    gene_sets: "FileElement | str | Sequence[str]",
    pres: list[Element],
) -> tuple[Any, list[Element]]:
    """Normalise the *gene_sets* argument for gseapy.

    Returns
    -------
    gene_sets_arg
        Value to pass to gseapy (path string, library name, or list).
    """
    if isinstance(gene_sets, Element):
        gmt_path = Path(gene_sets.primary.resolve())
        pres.append(gene_sets)
        return str(gmt_path), pres
    elif isinstance(gene_sets, (list, tuple)):
        return list(gene_sets), pres
    # plain string — Enrichr library name
    return gene_sets, pres


def _resolve_cls(
    classes: Sequence[str] | FileElement,
    pres: list[Element],
) -> tuple[Any, list[Element]]:
    """Normalise the *classes* argument for gseapy.gsea.

    Returns
    -------
    cls_arg
        Value for gseapy (list of labels or path string).
    pres
        List of :class:`Element` objects.
    """
    if isinstance(classes, Element):
        cls_path = Path(classes.path)
        pres.append(classes)
        return str(cls_path), pres
    return list(classes), pres


def _build_key_name(tag: ElementTag, subcommand: str) -> tuple[str, str]:
    name = tag.default_name
    key = f"{name}_{subcommand}"
    return key, name


def _select_gene_sets(
    input_gmt: Path,
    selected_genesets: Iterable[str],
    output_gmt: Path,
) -> None:
    selected_set_names = set(selected_genesets)
    with input_gmt.open() as infile, output_gmt.open("w") as outfile:
        for line in infile:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            set_name = parts[0]
            if set_name in selected_set_names:
                outfile.write(line)
