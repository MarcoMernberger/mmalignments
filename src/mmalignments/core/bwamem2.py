import subprocess
from pathlib import Path

import pypipegraph as ppg  # type: ignore
from mbf.externals.aligners.base import Aligner  # type: ignore


class BWAMem2(Aligner):
    """BWA-MEM2 aligner - faster reimplementation of BWA-MEM with identical output

    Usage example:
        from mbf.externals.aligners import BWAMem2
        from mbf import align, genomes

        # Create sample
        sample = align.Sample(
            sample_name="my_sample",
            input_strategy=align.strategies.FASTQsFromPrefix("path/to/fastqs_"),
            reverse_reads=False,
            pairing="paired",
        )

        # Align with BWA-MEM2
        aligner = BWAMem2()
        genome = genomes.mm39()
        aligned_sample = sample.align(
            aligner=aligner,
            genome=genome,
            aligner_parameters={
                "-k": 19,      # seed length
                "-T": 30,      # minimum score
            }
        )
    Options:
    Algorithm options:
        -o STR        Output SAM file name
        -k INT        minimum seed length [19]
        -w INT        band width for banded alignment [100]
        -d INT        off-diagonal X-dropoff [100]
        -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
        -y INT        seed occurrence for the 3rd round seeding [20]
        -c INT        skip seeds with more than INT occurrences [500]
        -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
        -W INT        discard a chain if seeded bases shorter than INT [0]
        -m INT        perform at most INT rounds of mate rescues for each read [50]
        -S            skip mate rescue
        -o            output file name missing
        -P            skip pairing; mate rescue performed unless -S also in use
    Scoring options:
    -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
    -B INT        penalty for a mismatch [4]
    -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
    -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
    -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
    -U INT        penalty for an unpaired read pair [17]
    Input/output options:
    -p            smart pairing (ignoring in2.fq)
    -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
    -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
    -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)
    -5            for split alignment, take the alignment with the smallest coordinate as primary
    -q            don't modify mapQ of supplementary alignments
    -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []
    -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
    -T INT        minimum score to output [30]
    -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
    -a            output all alignments for SE or unpaired PE
    -C            append FASTA/FASTQ comment to SAM output
    -V            output the reference FASTA header in the XR tag
    -Y            use soft clipping for supplementary alignments
    -M            mark shorter split hits as secondary
    -I FLOAT[,FLOAT[,INT[,INT]]]
                    specify the mean, standard deviation (10% of the mean if absent), max
                    (4 sigma from the mean if absent) and min of the insert size distribution.
                    FR orientation only. [inferred]
    """

    @property
    def parameter_options(self):
        parameter_options = {
            "-k": {"help": "minimum seed length", "default": 19, "type": int},
            "-w": {
                "help": "band width for banded alignment",
                "default": 100,
                "type": int,
            },
            "-d": {"help": "off-diagonal X-dropoff", "default": 100, "type": int},
            "-r": {
                "help": "look for internal seeds inside a seed longer than {-k} * FLOAT",
                "default": 1.5,
                "type": float,
            },
            "-y": {
                "help": "seed occurrence for the 3rd round seeding",
                "default": 20,
                "type": int,
            },
            "-c": {
                "help": "skip seeds with more than INT occurrences",
                "default": 500,
                "type": int,
            },
            "-D": {
                "help": "drop chains shorter than FLOAT fraction of the longest overlapping chain",
                "default": 0.50,
                "type": float,
            },
            "-W": {
                "help": "discard a chain if seeded bases shorter than INT",
                "default": 0,
                "type": int,
            },
            "-m": {
                "help": "perform at most INT rounds of mate rescues for each read",
                "default": 50,
                "type": int,
            },
            "-S": {
                "help": "skip mate rescue",
                "action": "store_true",
            },
            "-P": {
                "help": "skip pairing; mate rescue performed unless -S also in use",
                "action": "store_true",
            },
            "-A": {
                "help": "score for a sequence match, which scales options -TdBOELU unless overridden",
                "default": 1,
                "type": int,
            },
            "-B": {"help": "penalty for a mismatch", "default": 4, "type": int},
            "-O": {
                "help": "gap open penalties for deletions and insertions",
                "default": "6,6",
                "type": str,
            },
            "-E": {
                "help": "gap extension penalty; a gap of size k cost '{-O} + {-E}*k'",
                "default": "1,1",
                "type": str,
            },
            "-L": {
                "help": "penalty for 5'- and 3'-end clipping",
                "default": "5,5",
                "type": str,
            },
            "-U": {
                "help": "penalty for an unpaired read pair",
                "default": 17,
                "type": int,
            },
            "-p": {"help": "smart pairing (ignoring in2.fq)", "action": "store_true"},
            "-R": {
                "help": "read group header line such as '@RG\\tID:foo\\tSM:bar'",
                "default": None,
                "type": str,
            },
            "-H": {
                "help": "insert STR to header if it starts with @; or insert lines in FILE",
                "default": None,
                "type": str,
            },
            "-j": {
                "help": "treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)",
                "action": "store_true",
            },
            "-5": {
                "help": "for split alignment, take the alignment with the smallest coordinate as primary",
                "action": "store_true",
            },
            "-q": {
                "help": "don't modify mapQ of supplementary alignments",
                "action": "store_true",
            },
            "-K": {
                "help": "process INT input bases in each batch regardless of nThreads (for reproducibility)",
                "default": None,
                "type": int,
            },
            "-v": {
                "help": "verbose level: 1=error, 2=warning, 3=message, 4+=debugging",
                "default": 3,
                "type": int,
            },
            "-T": {"help": "minimum score to output", "default": 30, "type": int},
            "-h": {
                "help": "if there are <INT hits with score >80% of the max score, output all in XA",
                "default": "5,200",
                "type": str,
            },
            "-a": {
                "help": "output all alignments for SE or unpaired PE",
                "action": "store_true",
            },
            "-C": {
                "help": "append FASTA/FASTQ comment to SAM output",
                "action": "store_true",
            },
            "-V": {
                "help": "output the reference FASTA header in the XR tag",
                "action": "store_true",
            },
            "-Y": {
                "help": "use soft clipping for supplementary alignments",
                "action": "store_true",
            },
            "-M": {
                "help": "mark shorter split hits as secondary",
                "action": "store_true",
            },
            "-I": {
                "help": "specify the mean, standard deviation (10% of the mean if absent), max (4 sigma from the mean if absent) and min of the insert size distribution. FR orientation only.",
                "default": None,
                "type": str,
            },
        }
        return parameter_options

    @property
    def description(self):
        return "BWA-MEM2 aligner - faster reimplementation of BWA-MEM with identical output"

    @property
    def name(self):
        return "bwa-mem2"

    @property
    def primary_binary(self):
        return "bwa-mem2"

    @property
    def multi_core(self):
        return True

    def parameter_options_help(self):
        return "\n".join(
            [
                f"{k} {v['help']} (default: {v.get('default', 'None')})"
                for k, v in self.parameter_options.items()
            ]
        )

    def _aligner_build_cmd(self, output_dir, ncores, arguments):
        """Add threading parameter for multi-core support"""
        return arguments + ["-t", str(ncores)]

    def _index_build_cmd(self, output_dir, ncores, arguments):
        """Build index command without threading parameter"""
        return arguments

    def build_index_func(self, fasta_files, gtf_input_filename, output_fileprefix):
        """Build BWA-MEM2 index from fasta file(s)"""
        if isinstance(fasta_files, (str, Path)):
            fasta_files = [fasta_files]
        if len(fasta_files) > 1:
            raise ValueError("BWA-MEM2 can only build from a single fasta file")

        cmd = [
            "FROM_ALIGNER",
            self.primary_binary,
            "index",
            "-p",
            str((output_fileprefix / "bwa_mem2_index").absolute()),
        ]
        cmd.append(str(Path(fasta_files[0]).absolute()))
        print("in build_index_func:", cmd)
        return self.get_run_func(output_fileprefix, cmd, cwd=output_fileprefix)

    def align_job(
        self,
        input_fastq,
        paired_end_filename,
        index_job,
        output_bam_filename,
        parameters,
    ):
        """Create alignment job using bwa-mem2 mem algorithm"""

        def build_cmd():
            # Get index path from different ppg job types
            if hasattr(index_job, "target_folder"):  # ppg2 sharedmultifilegenjob
                index_basename = index_job.target_folder
            elif hasattr(index_job, "output_path"):  # ppg1 PrebuildJob
                index_basename = index_job.output_path
            else:
                index_basename = Path(index_job.files[0]).parent

            cmd = [
                "FROM_ALIGNER",
                self.primary_binary,
                "mem",
            ]
            for k, v in parameters.items():
                if k in self.parameter_options:
                    # Add algorithmic options with their values
                    if isinstance(v, self.parameter_options[k].get("type", str)):
                        cmd.append(str(k))
                    else:
                        raise ValueError(
                            f"Invalid type for parameter {k}: expected {self.parameter_options[k].get('type', str)}, got {type(v)}"
                        )
                else:
                    raise ValueError(f"Unknown parameter {k} for BWA-MEM2")
            # Index and input files
            cmd.append(str(Path(index_basename) / "bwa_mem2_index"))
            cmd.append(str(Path(input_fastq).absolute()))
            if paired_end_filename:
                cmd.append(str(Path(paired_end_filename).absolute()))
            output_bam_filename_sam = Path(output_bam_filename).with_suffix(".sam")
            cmd.extend(["-o", str(output_bam_filename_sam.absolute())])
            return cmd

        def sam_to_bam():
            """Convert SAM output to BAM format"""

            sam_file = output_bam_filename.with_suffix(".sam")
            # Use samtools for conversion (more robust than pysam for streaming)
            tmp_bam = Path(output_bam_filename).parent / "tmp.bam"

            # Convert SAM to BAM
            subprocess.check_call(
                [
                    "samtools",
                    "view",
                    "-b",  # output BAM
                    "-o",
                    str(tmp_bam),
                    str(sam_file),
                ]
            )

            # Sort BAM
            subprocess.check_call(
                [
                    "samtools",
                    "sort",
                    "-o",
                    str(output_bam_filename),
                    str(tmp_bam),
                ]
            )

            # Clean up temporary files
            tmp_bam.unlink()
            sam_file.unlink()

        job = self.run(
            Path(output_bam_filename).parent,
            build_cmd,
            cwd=Path(output_bam_filename).parent,
            call_afterwards=sam_to_bam,
            additional_files_created=output_bam_filename,
        )
        job.depends_on(
            ppg.ParameterInvariant(output_bam_filename, sorted(parameters.items()))
        )
        if self.multi_core:
            job.cores_needed = -1

        return job

    def get_version(self):
        """Get BWA-MEM2 version"""
        try:
            res = (
                subprocess.check_output([self.primary_binary, "version"])
                .decode("utf-8")
                .strip()
            )
            return res
        except subprocess.CalledProcessError:
            # Try alternative version command
            res = subprocess.check_output(
                [self.primary_binary], stderr=subprocess.STDOUT
            ).decode("utf-8")
            # Extract version from output
            for line in res.split("\n"):
                if "Version:" in line:
                    return line.split("Version:")[1].strip()
            return "unknown"

    def get_index_filenames(self):
        """List of files created by index building"""
        return [
            "bwa_mem2_index.0123",
            "bwa_mem2_index.amb",
            "bwa_mem2_index.ann",
            "bwa_mem2_index.bwt.2bit.64",
            "bwa_mem2_index.pac",
        ]

    def build_cmd(self, output_dir, ncores, arguments):
        if (
            not isinstance(arguments, list)
            or len(arguments) < 2
            or arguments[0] != "FROM_ALIGNER"
        ):
            raise ValueError(
                "Please call one of the following functions instead: .align_job, .build_index_job"
                + str(arguments)
            )
        if "index" in arguments:
            # For index building, we don't add threading parameters
            return self._index_build_cmd(output_dir, ncores, arguments[1:])
        else:
            return self._aligner_build_cmd(output_dir, ncores, arguments[1:])
