from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import requests  # type: ignore[import]
from pandas import DataFrame  # type: ignore[import]
from requests import Session  # type: ignore[import]
from requests.adapters import HTTPAdapter  # type: ignore[import]
from urllib3.util.retry import Retry  # type: ignore[import]

from mmalignments.models.artifacts import ArtifactSet
from mmalignments.models.data import Genome
from mmalignments.models.elements import (
    CallSpec,
    Element,
    FileSource,
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
from mmalignments.models.tags import (
    Method,
    State,
)
from mmalignments.services.io import (
    read_frame,
    write_frame,
)
from mmalignments.services.genomic import reverse_complement
from ..externals import (
    Runnable,
)

logger = logging.getLogger(__name__)


@dataclass
class ReferenceSpec:
    name: str
    species: str
    assembly: str
    chromosome: str
    start: int
    stop: int
    strand: int = 1  # 1 for forward strand, -1 for reverse strand


class EnsemblAPI:

    @element
    def amplicons(
        self,
        flanks: FileSource | Element,
        genome: Genome,
        *,
        tag: TagType | None = None,
        out: OutType | None = None,
        par: ParType | None = None,
        cfg: CfgType | None = None,
    ) -> Element:
        """Add amplicon sequences to a TSV frame file.

        Parameters
        ----------
        flanks : FileSource | Element
            Input TSV file with columns ``chr``, ``start``, ``stop``.
        genome : Genome
            Genome description (species, revision, prebuild prefix).
        tag : TagType, optional
            Tag to assign to the output Element.
        out : OutType, optional
            Output specification (folder, extension).
        par : ParType, optional
            Additional parameters for the operation.
        cfg : CfgType, optional
            Configuration for the operation.

        Returns
        -------
        Element
            Element whose artifact ``"tsv"`` is the path to the output TSV file.
        """
        tag = flanks.tag.bump(
            method=Method.CUSTOM, state=State.ANNOTATED  # noqa: E501
        ).resolve(tag)
        out = OutSpec.from_tag(tag).resolve(out)
        artifacts = ArtifactSet.from_outspec(out)
        par = par or Params()
        runner = self.generate_amplicons(
            output_file=out.file,
            frame_file=flanks.file,
            genome=genome,
            par=par,
        )
        key = Element.generate_key(tag, "GenomeAPi", "amplicons")
        return Element(
            key=key,
            run=runner,
            tag=tag,
            artifacts=artifacts,
            determinants=par.determinants,
            inputs=(flanks.file,),
            pres=(flanks,),
        )

    def generate_amplicons(
        self,
        output_file: Path,
        frame_file: Path,
        genome: Genome,
        par: ParType = Params(),
    ) -> Runnable:

        def run():
            df = read_frame(frame_file)
            df = df.drop_duplicates(subset=["target"], keep="first")
            df = df[
                [
                    "target",
                    "flank_start",
                    "flank_end",
                    "chr",
                    "start",
                    "stop",
                    "strand",
                    "guide",
                ]
            ].as_frame()
            df_out = self.add_amplicon_sequences(
                df,
                par=par,
            )
            write_frame(df_out, output_file)

        callspec = CallSpec(
            path=("generate_amplicons",),
            kwargs={
                "output_file": output_file,
                "frame_file": frame_file,
                "genome": genome,
                "par": par,
            },
        ).render()
        return Runnable(
            run,
            display=callspec,
        )

    def add_amplicon_sequences(
        self,
        df: DataFrame,
        *,
        par: ParType | None = None,
    ) -> DataFrame:

        par = Params().resolve(par)
        species = cast(str, par.get("species", "Mus_musculus"))
        revision = cast(str, par.get("revision", "GRCm39"))
        chr_col = cast(str, par.get("chr_col", "chr"))
        start_col = cast(str, par.get("start_col", "start"))
        stop_col = cast(str, par.get("stop_col", "stop"))
        output_col = cast(str, par.get("output_col", "amplicon"))
        strand_col = cast(str, par.get("strand_col", "strand"))
        sleep = cast(float, par.get("sleep", 0.1))
        timeout = cast(float, par.get("timeout", 60.0))
        result = df.copy()

        # ---------------------------------------------------------
        # Session mit automatischen Retries
        # ---------------------------------------------------------
        session = requests.Session()

        retry = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET"]),
            respect_retry_after_header=True,
        )

        adapter = HTTPAdapter(
            max_retries=retry,
            pool_connections=10,
            pool_maxsize=10,
        )

        session.mount("https://", adapter)

        session.headers.update(
            {
                "Accept": "text/plain",
            }
        )

        amplicons = []

        for idx, row in result.iterrows():

            _species = str(row.get("species", species))
            _revision = str(row.get("revision", revision))

            chromosome = str(row[chr_col])
            start = int(row[start_col])
            stop = int(row[stop_col])
            strand = int(row.get(strand_col, 1))
            reference_spec = ReferenceSpec(
                "name", _species, _revision, chromosome, start, stop, strand
            )
            sequence = self.fetch_sequence(
                reference_spec, session, sleep=sleep, timeout=timeout, idx=idx
            )
            amplicons.append(sequence)

            if sleep > 0:
                time.sleep(sleep)

        result[output_col] = amplicons

        return result

    def open_session(self) -> Session:
        session = requests.Session()

        retry = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET"]),
            respect_retry_after_header=True,
        )

        adapter = HTTPAdapter(
            max_retries=retry,
            pool_connections=10,
            pool_maxsize=10,
        )

        session.mount("https://", adapter)

        session.headers.update(
            {
                "Accept": "text/plain",
            }
        )
        return session

    def fetch_sequence(
        self,
        reference: ReferenceSpec,
        session: requests.Session,
        timeout: float = 60.0,
        idx: Any = 0,
    ) -> str | None:
        # Ensembl needs start <= stop
        if reference.start > reference.stop:
            reference.start, reference.stop = reference.stop, reference.start

        url = (
            f"https://rest.ensembl.org/sequence/region/"
            f"{reference.species}/{reference.chromosome}:{reference.start}..{reference.stop}:1"
        )
        response = None
        try:
            response = session.get(
                url,
                params={
                    "coord_system_version": reference.assembly,
                },
                timeout=(10, timeout),  # connect timeout, read timeout
            )

            response.raise_for_status()

            sequence = response.text.strip().upper()
            if reference.strand == -1:
                sequence = reverse_complement(sequence)
            return sequence

        except requests.Timeout as e:
            print(
                f"Timeout bei {idx}: "
                f"{reference.species} {reference.chromosome}:{reference.start}-{reference.stop}: {e}"
            )
            return None

        except requests.HTTPError as e:
            print(
                f"HTTP error bei {idx}: "
                f"{reference.species} {reference.chromosome}:{reference.start}-{reference.stop}: "
            )
            if response:
                print(f"{response.status_code} {e}")
            return None
        except requests.RequestException as e:
            print(
                f"Request error bei {idx}: "
                f"{reference.species} {reference.chromosome}:{reference.start}-{reference.stop}: {e}"
            )
            return None
