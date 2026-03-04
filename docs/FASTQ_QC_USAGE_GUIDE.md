# Pre-Alignment FASTQ Quality Control - Usage Guide

This document provides comprehensive usage examples for the FASTQ quality control workflow in the mmalignments package.

## Overview

The pre-alignment QC workflow consists of the following steps:

1. **Input Check** - Validates FASTQ structure and paired-end ID matching
2. **fastp** - Quality filtering, adapter trimming, and QC report generation
3. **FastQC (raw)** - Quality assessment of raw FASTQ files
4. **FastQC (cleaned)** - Quality assessment of fastp-cleaned FASTQ files
5. **QC Summary** - Parses and aggregates metrics from fastp and FastQC
6. **MultiQC (optional)** - Project-level aggregation of all QC reports

## Quick Start

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
from pathlib import Path

# Create sample
sample = Sample(
    name="sample1",
    pairing="paired",
    fastq_r1_path="incoming/sample1_R1.fastq.gz",
    fastq_r2_path="incoming/sample1_R2.fastq.gz"
)

# Run complete pre-alignment QC workflow
qc_elements = pre_alignment_qc(
    sample=sample,
    qc_root=Path("results/qc"),
    threads=8
)

# Execute all QC steps
for elem in qc_elements:
    elem.run()

print(f"QC completed! Generated {len(qc_elements)} QC elements")
```

## Detailed Workflow Example

### 1. Single Sample QC

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
from pathlib import Path

# Define sample
sample = Sample(
    name="tumor_sample_1",
    pairing="paired",
    fastq_r1_path="incoming/251215_A02023/tumor_1_R1.fastq.gz",
    fastq_r2_path="incoming/251215_A02023/tumor_1_R2.fastq.gz"
)

# Run QC with custom parameters
qc_elements = pre_alignment_qc(
    sample=sample,
    qc_root="results/qc",
    threads=8,
    parameters={
        "fastp": {
            "--qualified_quality_phred": 20,
            "--length_required": 36,
        },
        "fastqc": {
            "--nogroup": True,
        }
    }
)

# Execute workflow
for elem in qc_elements:
    print(f"Running: {elem.name}")
    elem.run()
    print(f"  Status: OK")
    print(f"  Artifacts: {list(elem.artifacts.keys())}")
```

### 2. Multiple Samples in a Loop

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
from pathlib import Path

# Define multiple samples
samples = [
    Sample(
        name="sample1",
        pairing="paired",
        fastq_r1_path="incoming/sample1_R1.fastq.gz",
        fastq_r2_path="incoming/sample1_R2.fastq.gz"
    ),
    Sample(
        name="sample2",
        pairing="paired",
        fastq_r1_path="incoming/sample2_R1.fastq.gz",
        fastq_r2_path="incoming/sample2_R2.fastq.gz"
    ),
    Sample(
        name="sample3",
        pairing="paired",
        fastq_r1_path="incoming/sample3_R1.fastq.gz",
        fastq_r2_path="incoming/sample3_R2.fastq.gz"
    ),
]

qc_root = Path("results/qc")
all_qc_elements = []

# Process each sample
for sample in samples:
    print(f"\nProcessing sample: {sample.name}")
    
    qc_elems = pre_alignment_qc(
        sample=sample,
        qc_root=qc_root,
        threads=8
    )
    
    # Execute QC for this sample
    for elem in qc_elems:
        elem.run()
    
    all_qc_elements.extend(qc_elems)
    print(f"  Completed {len(qc_elems)} QC steps")

print(f"\nAll samples processed! Total QC elements: {len(all_qc_elements)}")
```

### 3. With Project-Level MultiQC Aggregation

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
from mmalignments.models.qc.multiqc import MultiQC
from pathlib import Path

# Define samples
samples = [
    Sample(name="sample1", pairing="paired",
           fastq_r1_path="incoming/sample1_R1.fastq.gz",
           fastq_r2_path="incoming/sample1_R2.fastq.gz"),
    Sample(name="sample2", pairing="paired",
           fastq_r1_path="incoming/sample2_R1.fastq.gz",
           fastq_r2_path="incoming/sample2_R2.fastq.gz"),
]

qc_root = Path("results/qc")
all_qc_elements = []

# Process all samples
for sample in samples:
    qc_elems = pre_alignment_qc(
        sample=sample,
        qc_root=qc_root,
        threads=8
    )
    
    for elem in qc_elems:
        elem.run()
    
    all_qc_elements.extend(qc_elems)

# Create MultiQC aggregation report
multiqc = MultiQC()
multiqc_elem = multiqc.aggregate(
    qc_root=qc_root,
    qc_elements=all_qc_elements,
    output_dir=qc_root / "multiqc",
    report_name="pre_alignment_multiqc.html"
)

# Generate MultiQC report
multiqc_elem.run()

print(f"\nMultiQC report generated: {multiqc_elem.artifacts['multiqc_html']}")
```

### 4. Using Individual Components

You can also use the QC components individually for more control:

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import build_input_check
from mmalignments.models.qc.fastp import FastP
from mmalignments.models.qc.fastqc import FastQC
from mmalignments.models.qc.parsers import build_qc_summary
from pathlib import Path

sample = Sample(
    name="sample1",
    pairing="paired",
    fastq_r1_path="incoming/sample1_R1.fastq.gz",
    fastq_r2_path="incoming/sample1_R2.fastq.gz"
)

qc_dir = Path("results/qc") / sample.name

# 1. Input check
input_check = build_input_check(
    sample=sample,
    out_dir=qc_dir / "input_check",
    n_records=10_000
)
input_check.run()

# 2. fastp
fastp = FastP()
fastp_elem = fastp.fastp(
    sample=sample,
    folder=qc_dir / "fastp",
    threads=8,
    parameters={
        "--qualified_quality_phred": 20,
        "--length_required": 36,
    }
)
fastp_elem.run()

# 3. FastQC on raw
fastqc = FastQC()
fastqc_raw = fastqc.fastqc(
    sample=sample,
    folder=qc_dir / "fastqc_raw",
    label="raw",
    threads=4
)
fastqc_raw.run()

# 4. FastQC on cleaned (using fastp output)
fastqc_cleaned = fastqc.fastqc(
    sample=fastp_elem,  # Use fastp element (contains cleaned FASTQs)
    folder=qc_dir / "fastqc_cleaned",
    label="cleaned",
    threads=4
)
fastqc_cleaned.run()

# 5. QC summary
qc_summary = build_qc_summary(
    sample=sample,
    fastp_json=fastp_elem.artifacts["fastp_json"],
    fastqc_raw_r1_zip=fastqc_raw.artifacts["fastqc_zip_r1"],
    fastqc_raw_r2_zip=fastqc_raw.artifacts["fastqc_zip_r2"],
    fastqc_cleaned_r1_zip=fastqc_cleaned.artifacts["fastqc_zip_r1"],
    fastqc_cleaned_r2_zip=fastqc_cleaned.artifacts["fastqc_zip_r2"],
    out_json=qc_dir / f"{sample.name}_qc.json",
    out_tsv=qc_dir / f"{sample.name}_qc.tsv"
)
qc_summary.run()

print("All QC steps completed individually!")
```

### 5. Integration with Executor/DAG

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
from mmalignments.models.executor import Executor
from pathlib import Path

# Define samples
samples = [
    Sample(name="sample1", pairing="paired",
           fastq_r1_path="incoming/sample1_R1.fastq.gz",
           fastq_r2_path="incoming/sample1_R2.fastq.gz"),
    Sample(name="sample2", pairing="paired",
           fastq_r1_path="incoming/sample2_R1.fastq.gz",
           fastq_r2_path="incoming/sample2_R2.fastq.gz"),
]

# Create executor
executor = Executor()

# Add all QC elements to executor
for sample in samples:
    qc_elements = pre_alignment_qc(
        sample=sample,
        qc_root=Path("results/qc"),
        threads=8
    )
    
    for elem in qc_elements:
        executor.add(elem)

# Execute entire DAG
executor.run()

print("All QC completed via Executor!")
```

### 6. Single-End Reads

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
from pathlib import Path

# Define single-end sample
sample = Sample(
    name="single_end_sample",
    pairing="single",
    fastq_r1_path="incoming/single_sample.fastq.gz"
)

# Run QC (automatically handles single-end)
qc_elements = pre_alignment_qc(
    sample=sample,
    qc_root=Path("results/qc"),
    threads=8
)

for elem in qc_elements:
    elem.run()

print("Single-end QC completed!")
```

## Output Structure

After running the QC workflow, the output directory structure will be:

```
results/qc/
├── sample1/
│   ├── input_check/
│   │   ├── sample1_input_check.json
│   │   └── sample1_input_check.txt
│   ├── fastp/
│   │   ├── sample1_cleaned_R1.fastq.gz
│   │   ├── sample1_cleaned_R2.fastq.gz
│   │   ├── sample1_fastp.json
│   │   └── sample1_fastp.html
│   ├── fastqc_raw/
│   │   ├── sample1_R1_fastqc.html
│   │   ├── sample1_R1_fastqc.zip
│   │   ├── sample1_R2_fastqc.html
│   │   └── sample1_R2_fastqc.zip
│   ├── fastqc_cleaned/
│   │   ├── sample1_cleaned_R1_fastqc.html
│   │   ├── sample1_cleaned_R1_fastqc.zip
│   │   ├── sample1_cleaned_R2_fastqc.html
│   │   └── sample1_cleaned_R2_fastqc.zip
│   ├── sample1_qc.json
│   └── sample1_qc.tsv
├── sample2/
│   └── ... (same structure)
└── multiqc/
    ├── pre_alignment_multiqc.html
    └── multiqc_data/
```

## Key Metrics to Monitor

### From Input Check
- `status`: PASS/FAIL/ERROR
- `r1_records_checked`: Number of R1 records validated
- `r2_records_checked`: Number of R2 records validated (if paired)
- `ids_match`: Whether R1/R2 IDs match (paired-end only)
- `id_mismatches`: Number of ID mismatches found

### From fastp JSON
- `total_reads_before`/`total_reads_after`: Read counts
- `q20_rate_after`, `q30_rate_after`: Quality rates
- `gc_content_after`: GC content percentage
- `duplication_rate`: Duplication rate
- `adapter_trimmed_reads`: Number of reads with adapters

### From FastQC
- `per_base_sequence_quality`: PASS/WARN/FAIL
- `per_sequence_quality_scores`: PASS/WARN/FAIL
- `per_base_sequence_content`: PASS/WARN/FAIL
- `adapter_content`: PASS/WARN/FAIL
- `total_sequences`: Total number of sequences
- `gc_content`: %GC content

### From QC Summary Files
- `sample1_qc.json`: Machine-readable JSON with all metrics
- `sample1_qc.tsv`: Tab-delimited file for easy parsing/plotting

## Advanced Configuration

### Custom fastp Parameters

```python
qc_elements = pre_alignment_qc(
    sample=sample,
    qc_root="results/qc",
    threads=8,
    parameters={
        "fastp": {
            "--qualified_quality_phred": 25,      # Higher quality threshold
            "--length_required": 50,              # Minimum read length
            "--cut_mean_quality": 20,             # Sliding window quality
            "--cut_window_size": 4,               # Window size for sliding
            "--dedup": True,                       # Enable deduplication
            "--trim_poly_g": True,                 # Trim poly-G tails
        }
    }
)
```

### Custom FastQC Parameters

```python
qc_elements = pre_alignment_qc(
    sample=sample,
    qc_root="results/qc",
    threads=8,
    parameters={
        "fastqc": {
            "--nogroup": True,      # Don't group bases in output
            "--kmers": 7,           # K-mer size for overrepresentation
        }
    }
)
```

## Error Handling

```python
from mmalignments.models.data import Sample
from mmalignments.models.qc.qc import pre_alignment_qc
import logging

# Enable logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

sample = Sample(
    name="sample1",
    pairing="paired",
    fastq_r1_path="incoming/sample1_R1.fastq.gz",
    fastq_r2_path="incoming/sample1_R2.fastq.gz"
)

try:
    qc_elements = pre_alignment_qc(
        sample=sample,
        qc_root="results/qc",
        threads=8
    )
    
    for elem in qc_elements:
        try:
            logger.info(f"Running: {elem.name}")
            elem.run()
            logger.info(f"  Completed: {elem.name}")
        except Exception as e:
            logger.error(f"  Failed: {elem.name} - {e}")
            # Continue with other elements or break
            continue
            
except Exception as e:
    logger.error(f"Failed to create QC workflow: {e}")
```

## Best Practices

1. **Always run input_check first** - It validates FASTQ format before expensive QC
2. **Use adequate threading** - fastp and FastQC benefit from 4-8 threads
3. **Check QC summaries** - Review JSON/TSV files for outliers
4. **Use MultiQC for projects** - Aggregates results across samples
5. **Monitor disk space** - Cleaned FASTQs and QC reports can be large
6. **Archive raw FASTQs** - Keep originals after validation

## Troubleshooting

### Issue: ID mismatch warnings in paired-end data
**Solution**: Check if files are properly paired. Some sequencers use different ID formats. The input_check will report mismatches.

### Issue: fastp reports low quality
**Solution**: Adjust `--qualified_quality_phred` threshold or inspect with FastQC first.

### Issue: FastQC shows adapter content
**Solution**: fastp should trim adapters automatically. Check fastp JSON for `adapter_trimmed_reads`.

### Issue: MultiQC doesn't find all reports
**Solution**: Ensure all QC elements completed successfully and outputs are in the `qc_root` directory.

## Next Steps

After pre-alignment QC:

1. Review QC summaries and MultiQC report
2. Decide on filtering thresholds
3. Use cleaned FASTQs from fastp for alignment
4. Proceed with alignment workflow (e.g., BWA-MEM2)
5. Run post-alignment QC (see `POST_MAPPING_QC_USAGE.md`)
