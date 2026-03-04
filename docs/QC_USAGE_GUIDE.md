# Post-Mapping QC Functions - Usage Guide

This document provides usage examples for the newly added QC functions in the mmalignments package.

## Overview

The following QC functions have been added:

### In `mmalignments/models/callers/gatk.py`:

#### Utility Functions
- `bed_to_interval_list()` - Convert BED files to Picard interval_list format

#### Picard Metrics (low-level and high-level)
- `alignment_summary_metrics()` / `alignment_summary()` - Alignment summary metrics
- `insert_size_metrics()` / `insert_size()` - Insert size distribution
- `hs_metrics()` / `hs()` - Hybrid selection (capture) metrics

### In `mmalignments/models/qc/mosdepth.py`:
- `coverage()` / `depth()` - Coverage depth calculation on targets

### In `mmalignments/models/qc/samtools.py` (new file):
- `flagstat_bam()` / `flagstat()` - Basic alignment statistics
- `stats_bam()` / `stats()` - Comprehensive alignment statistics

### In `mmalignments/models/qc/multiqc.py`:
- `aggregate()` - Aggregate all QC reports into a single HTML report

### In `mmalignments/models/qc/qc.py` (new file):
- `post_mapping_qc()` - Convenience function to run all QC steps
- `post_mapping_qc_with_multiqc()` - Same as above, plus MultiQC aggregation

## Pattern

All functions follow the same low-level / high-level pattern:

**Low-level functions** (e.g., `alignment_summary_metrics()`, `coverage()`):
- Accept file paths and basic parameters
- Return a zero-argument callable (runner)
- Used when you need fine-grained control

**High-level functions** (e.g., `alignment_summary()`, `depth()`):
- Accept Elements (MappedElement, Genome, etc.)
- Return an Element with dependencies
- Automatically construct output paths
- Used in DAG workflows

## Usage Examples

### 1. Individual QC Steps (High-level)

```python
from mmalignments.models.callers.gatk import GATK
from mmalignments.models.qc.mosdepth import Mosdepth
from mmalignments.models.qc.samtools import SamtoolsQC
from mmalignments.models.data import Genome

# Assuming you have a MappedElement from alignment
# mapped = aligner.align(...)

genome = Genome(species="Mus_musculus", revision=115, prebuild_prefix="...")
gatk = GATK(primary_binary="code/gatk-4.6.2.0/gatk")
mosdepth = Mosdepth()
samtools_qc = SamtoolsQC()

# 1. Alignment summary metrics
alignment_qc = gatk.alignment_summary(mapped, genome)

# 2. Insert size metrics
insert_qc = gatk.insert_size(mapped)

# 3. Hybrid selection metrics
hs_qc = gatk.hs(
    mapped=mapped,
    reference=genome,
    bait_intervals="targets.interval_list",
    target_intervals="targets.padded.interval_list"
)

# 4. Coverage depth with mosdepth
depth_qc = mosdepth.depth(
    mapped=mapped,
    targets="targets.padded.bed",
    threads=8
)

# 5. Samtools flagstat
flagstat_qc = samtools_qc.flagstat(mapped)

# 6. Samtools stats
stats_qc = samtools_qc.stats(mapped)

# Execute any of them
alignment_qc.run()
insert_qc.run()
# ... etc
```

### 2. Using the Convenience Function

```python
from mmalignments.models.qc.qc import post_mapping_qc
from mmalignments.models.data import Genome

genome = Genome(species="Mus_musculus", revision=115, prebuild_prefix="...")

# Run all QC steps at once
qc_elements = post_mapping_qc(
    mapped=mapped_sample,
    reference=genome,
    targets="incoming/targets.padded.bed",
    bait_intervals="targets.interval_list",
    target_intervals="targets.padded.interval_list",
    threads_mosdepth=8,
    run_alignment_summary=True,
    run_insert_size=True,
    run_hs_metrics=True,
    run_mosdepth=True,
    run_samtools_stats=True,
    run_samtools_flagstat=True
)

# Execute all QC steps
for elem in qc_elements:
    elem.run()
```

### 3. With MultiQC Aggregation

```python
from mmalignments.models.qc.qc import post_mapping_qc_with_multiqc
from mmalignments.models.data import Genome

genome = Genome(species="Mus_musculus", revision=115, prebuild_prefix="...")

# Run all QC steps and create MultiQC report
qc_elems, multiqc_elem = post_mapping_qc_with_multiqc(
    mapped=mapped_sample,
    reference=genome,
    targets="incoming/targets.padded.bed",
    bait_intervals="targets.interval_list",
    target_intervals="targets.padded.interval_list",
    threads_mosdepth=8,
    multiqc_output_dir="results/qc/multiqc",
    multiqc_report_name="project_qc_report.html"
)

# Execute all QC steps
for elem in qc_elems:
    elem.run()

# Generate MultiQC report
multiqc_elem.run()
```

### 4. Converting BED to interval_list

```python
from mmalignments.models.callers.gatk import GATK
from mmalignments.models.data import Genome

gatk = GATK(primary_binary="code/gatk-4.6.2.0/gatk")
genome = Genome(species="Mus_musculus", revision=115, prebuild_prefix="...")

# Convert BED to interval_list for Picard
converter = gatk.bed_to_interval_list(
    input_bed="incoming/targets.bed",
    output_interval_list="incoming/targets.interval_list",
    sequence_dict=f"{genome.fasta.parent}/genome.dict"
)
converter()

# For padded targets
converter_padded = gatk.bed_to_interval_list(
    input_bed="incoming/targets.padded.bed",
    output_interval_list="incoming/targets.padded.interval_list",
    sequence_dict=f"{genome.fasta.parent}/genome.dict"
)
converter_padded()
```

### 5. Low-level Usage (Fine-grained Control)

```python
from mmalignments.models.callers.gatk import GATK
from pathlib import Path

gatk = GATK(primary_binary="code/gatk-4.6.2.0/gatk")

# Low-level call with explicit paths
runner = gatk.alignment_summary_metrics(
    reference=Path("genome.fa"),
    input_bam=Path("sample.bam"),
    output_metrics=Path("qc/alignment_summary.txt")
)
runner()
```

### 6. Integration with Executor/DAG

```python
from mmalignments.models.executor import Executor
from mmalignments.models.qc.qc import post_mapping_qc_with_multiqc

# Create QC workflow
qc_elems, multiqc_elem = post_mapping_qc_with_multiqc(
    mapped=mapped_sample,
    reference=genome,
    targets=targets_bed,
    bait_intervals=targets_interval_list,
    target_intervals=targets_padded_interval_list,
    threads_mosdepth=8
)

# Add to executor
executor = Executor()
for elem in qc_elems:
    executor.add(elem)
executor.add(multiqc_elem)

# Execute DAG
executor.run()
```

## Output Structure

By default, all QC outputs are organized under the BAM file's parent directory:

```
results/aligned/<genome>/<sample>/
├── <sample>_sorted.bam
├── <sample>_sorted.bam.bai
└── qc/
    ├── picard/
    │   ├── <sample>_sorted_alignment_summary.txt
    │   ├── <sample>_sorted_insert_size.txt
    │   ├── <sample>_sorted_insert_size.pdf
    │   └── <sample>_sorted_hs_metrics.txt
    ├── mosdepth/
    │   ├── <sample>_sorted.targets.mosdepth.summary.txt
    │   ├── <sample>_sorted.targets.regions.bed.gz
    │   └── <sample>_sorted.targets.mosdepth.global.dist.txt
    ├── samtools/
    │   ├── <sample>_sorted_flagstat.txt
    │   └── <sample>_sorted_stats.txt
    └── multiqc/
        ├── multiqc_report.html
        └── multiqc_data/
```

## Key Metrics to Monitor

### From Picard CollectAlignmentSummaryMetrics
- `PCT_PF_READS_ALIGNED` - Percentage of reads aligned
- `PF_MISMATCH_RATE` - Mismatch rate
- `MEAN_READ_LENGTH` - Average read length

### From Picard CollectInsertSizeMetrics
- `MEAN_INSERT_SIZE` - Average fragment size
- `STANDARD_DEVIATION` - Insert size variation
- Check histogram PDF for distribution

### From Picard CollectHsMetrics
- `PCT_SELECTED_BASES` - On-target percentage
- `MEAN_TARGET_COVERAGE` - Average coverage on targets
- `FOLD_80_BASE_PENALTY` - Coverage uniformity (lower is better)
- `PCT_TARGET_BASES_10X`, `PCT_TARGET_BASES_20X` - Target coverage breadth

### From mosdepth
- Check `*.mosdepth.summary.txt` for overall depth statistics
- Check `*.regions.bed.gz` for per-region depth
- Check `*.mosdepth.global.dist.txt` for depth distribution

### From samtools flagstat
- Total reads, mapped reads, properly paired
- Duplicates (if marked)

### From samtools stats
- Comprehensive alignment statistics
- Error rates, insert size details

## Notes

- **Sequence Dictionary**: GATK/Picard tools require a `.dict` file for the reference.
  Create it with `samtools dict genome.fa > genome.dict` or `gatk CreateSequenceDictionary`.

- **BED vs interval_list**: mosdepth uses BED format, Picard tools use interval_list.
  Use `gatk.bed_to_interval_list()` to convert.

- **Threading**: Only mosdepth benefits from multiple threads in these QC steps.
  Set `threads_mosdepth=8` or similar for faster coverage calculation.

- **MultiQC**: Automatically detects output files from supported tools. Place all
  QC outputs under a common directory (e.g., `results/qc/`) for best results.

- **Dependencies**: All high-level functions create Elements with proper dependency
  tracking. The DAG ensures QC runs after alignment completes.
