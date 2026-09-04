"""Wrapper for CRISPResso in the mmalignments Element design.

The wrapper follows the same pattern as :class:`MmFqCount`:

- high-level ``@element`` methods that construct :class:`Element` objects
- low-level ``@subroutine`` methods that build runnable subprocess calls
"""

from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Mapping

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.elements import (
    Element,
    FileSource,
    NextGenSample,
    element,
)
from mmalignments.models.overlay import (
    CfgType,
    OutSpec,
    OutType,
    Params,
    ParType,
    TagType,
)
from mmalignments.models.parameters import ParamSet
from mmalignments.models.tags import Method, Stage, State
from mmalignments.services.io import read_frame

from ..externals import External, SubroutineIn, subroutine

logger = logging.getLogger(__name__)


def _overlay_signature_determinants(params: ParType | None) -> tuple[str, ...]:
    if params is None:
        return ()
    merged = params.to_dict()
    tokens: list[str] = []
    for key in sorted(merged.keys()):
        value = merged[key]
        if value is None or value is False or value == []:
            continue
        tokens.append(f"{key}={value}")
    return tuple(tokens)


def _params_to_cli(params: ParType | None) -> list[str]:
    """Convert overlay params to CLI flags for CRISPResso.

    Keys are converted from ``snake_case`` to ``--kebab-case``.
    ``True`` values are emitted as flag-only options.
    ``False``/``None``/empty lists are skipped.
    """
    if params is None:
        return []
    cli: list[str] = []
    for key in sorted(params.keys()):
        value = params.get(key)
        if value is None or value is False or value == []:
            continue
        flag = "--" + str(key)
        if value is True:
            cli.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            cli.extend([flag, ",".join(str(v) for v in value)])
            continue
        cli.extend([flag, str(value)])
    return cli


class Crispresso(External):
    """Wrapper for the ``CRISPResso`` CLI.

    The main entry point is :meth:`run`, which executes CRISPResso for one
    sample and records the HTML report as primary artifact.
    """

    def __init__(
        self,
        name: str = "crispresso",
        primary_binary: str = "CRISPResso",
        version: str | None = None,
        source: str = "https://github.com/pinellolab/CRISPResso2",
        parameters: Mapping[str, ParamSet] | ParamSet | str | Path | None = None,
    ) -> None:
        parameters_file = Path(__file__).parent / f"{Path(__file__).stem}.json"
        parameters = parameters or parameters_file

        super().__init__(
            name=name,
            primary_binary=primary_binary,
            version=version,
            source=source,
            parameters=parameters or {},
        )

    def get_version(self, fallback: str | None = None) -> str | None:
        "Detect ``CRISPResso`` version from ``CRISPResso --version`` output."
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

    def default_outdir(self, sample_name: str) -> Path:
        """Return default output directory for CRISPResso reports."""
        return Path("results") / "crispresso" / self.version_name / sample_name

    @element
    def run(
        self,
        sample: Element | NextGenSample,
        amplicons: Element | FileSource,
        *,
        target: str | None = None,
        report_name: str | None = None,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Run CRISPResso for one sample.

        Parameters
        ----------
        sample : Element | NextGenSample
            Input sample with FASTQ artifact(s).
        amplicons : Element | FileSource
            Comma-separated amplicon sequence(s) for ``--amplicon_seq``.
        report_name : str | None
            Report stem used by CRISPResso. Defaults to ``sample.root``.
        tag, out, par, cfg
            Standard mmalignments element overlays.
        """
        fastq_r1 = sample.primary.r1
        fastq_r2 = sample.primary.r2 if hasattr(sample.primary, "r2") else None
        sample_root = (
            sample.root if hasattr(sample, "root") else sample.tag.root
        )  # noqa: E501
        resolved_report_name = report_name or sample_root
        par = par or Params()
        tag = sample.tag.bump(
            stage=Stage.ANALYSIS,
            method=Method.CRISPRESSO,
            state=State.REPORT,
        ).resolve(
            tag,
        )
        out = OutSpec(
            stem=f"CRISPResso_on_{resolved_report_name}",
            folder=self.default_outdir(sample_root),
            ext="html",
        ).resolve(out)

        runner = self.run_crispresso(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            output_folder=out.folder,
            report_name=resolved_report_name,
            amplicons_path=amplicons.file,
            target=target,
            par=par,
            cfg=cfg,
        )

        key = Element.generate_key(
            tag,
            "run",
            report_name=resolved_report_name,
        )

        determinants = (
            f"amplicon_file={amplicons.file}",
            f"report_name={resolved_report_name}",
        ) + par.determinants
        if target:
            determinants += (f"target={target}",)
        inputs = (fastq_r1, fastq_r2) if fastq_r2 else (fastq_r1,)
        pres = (
            (sample, amplicons)
            if isinstance(sample, Element)
            else sample.pres + (amplicons,)
        )

        return Element(
            key,
            runner,
            tag=tag,
            artifacts=ArtifactSet(out.file, primary_name="html"),
            determinants=determinants,
            inputs=inputs,
            pres=pres,
        )

    @subroutine
    def run_crispresso(
        self,
        fastq_r1: Path,
        fastq_r2: Path | None,
        output_folder: Path,
        report_name: str,
        amplicons_path: Path,
        target: str | None = None,
        *,
        amplicon_col: str | None = None,
        target_col: str | None = None,
        guide_col: str | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> SubroutineIn:
        """Low-level wrapper for one CRISPResso invocation."""
        par = Params().resolve(par)
        amplicon_col = amplicon_col or "amplicon"
        target_col = target_col or "target"
        guide_col = guide_col or "guide"

        amplicons = read_frame(amplicons_path)
        if target:
            amplicons = amplicons[amplicons[target_col] == target]
        amplicon_seq = amplicons[amplicon_col].str.cat(sep=",")
        amplicon_name = amplicons[target_col].str.cat(sep=",")
        guide_seq = amplicons[guide_col].str.cat(sep=",")
        arguments = [
            "--fastq_r1",
            str(Path(fastq_r1).resolve()),
            "--amplicon_seq",
            str(amplicon_seq),
            "--guide_seq",
            str(guide_seq),
            "--amplicon_name",
            str(amplicon_name),
            "--output_folder",
            str(Path(output_folder).resolve()),
            "--name",
            str(report_name),
        ]

        if fastq_r2 is not None:
            arguments.extend(["--fastq_r2", str(Path(fastq_r2).resolve())])
        # arguments.extend(_params_to_cli(par))

        html_report = (
            Path(output_folder).resolve() / f"CRISPResso_on_{report_name}.html"
        )
        in_paths = [Path(fastq_r1).resolve()]
        if fastq_r2 is not None:
            in_paths.append(Path(fastq_r2).resolve())
        out_paths = [html_report]

        return (
            arguments,
            None,
            in_paths,
            out_paths,
            None,
            None,
            None,
        )
