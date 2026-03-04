"""Reporting services for generating final reports and summaries."""

import pypipegraph2 as ppg
from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd
import json
from ..models.variant import MutationalLoad
from ..models.qc_metrics import SampleQCReport


class ReportingService:
    """Service for generating analysis reports."""

    def __init__(self):
        """Initialize reporting service."""
        pass

    def generate_mutational_load_report(
        self,
        mutational_loads: List[MutationalLoad],
        output_file: Path,
        qc_reports: Optional[Dict[str, SampleQCReport]] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Generate comprehensive mutational load report.
        
        Args:
            mutational_loads: List of MutationalLoad objects
            output_file: Output TSV file
            qc_reports: Optional QC reports for each sample
            
        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def generate():
            rows = []
            for ml in mutational_loads:
                row = ml.to_dict()
                
                # Add QC metrics if available
                if qc_reports and ml.sample_name in qc_reports:
                    qc = qc_reports[ml.sample_name]
                    qc_dict = qc.to_dict()
                    row.update(qc_dict)
                
                rows.append(row)
            
            # Create DataFrame and save
            df = pd.DataFrame(rows)
            
            # Reorder columns for readability
            primary_cols = [
                'sample', 'normal', 'n_snv', 'n_insertion', 'n_deletion',
                'n_total', 'callable_mb', 'mutational_load'
            ]
            
            ordered_cols = [c for c in primary_cols if c in df.columns]
            remaining_cols = [c for c in df.columns if c not in ordered_cols]
            df = df[ordered_cols + remaining_cols]
            
            df.to_csv(output_file, sep='\t', index=False, float_format='%.4f')

        job = ppg.FileGeneratingJob(output_file, generate)
        return job

    def generate_sample_summary_report(
        self,
        sample_name: str,
        qc_report: SampleQCReport,
        mutational_load: Optional[MutationalLoad],
        output_file: Path,
    ) -> ppg.FileGeneratingJob:
        """
        Generate individual sample summary report.
        
        Args:
            sample_name: Sample name
            qc_report: QC report for sample
            mutational_load: Mutational load data (if tumor)
            output_file: Output JSON file
            
        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def generate():
            report = {
                "sample_name": sample_name,
                "qc_metrics": qc_report.to_dict(),
            }
            
            if mutational_load:
                report["mutational_load"] = mutational_load.to_dict()
            
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)

        job = ppg.FileGeneratingJob(output_file, generate)
        return job

    def generate_variant_summary(
        self,
        annotated_vcf: Path,
        output_file: Path,
        sample_name: str,
    ) -> ppg.FileGeneratingJob:
        """
        Generate variant summary from annotated VCF.
        
        Args:
            annotated_vcf: VEP-annotated VCF file
            output_file: Output summary file
            sample_name: Sample name
            
        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def generate():
            from ..services.annotation import AnnotationService
            
            ann_service = AnnotationService()
            counts = ann_service.count_variants_by_effect(annotated_vcf)
            
            # Write summary
            with open(output_file, 'w') as f:
                f.write(f"Sample: {sample_name}\n")
                f.write(f"VCF: {annotated_vcf}\n\n")
                f.write("Variant Counts by Effect:\n")
                f.write("-" * 40 + "\n")
                for effect, count in counts.items():
                    f.write(f"{effect:20s}: {count:8d}\n")

        job = ppg.FileGeneratingJob(output_file, generate)
        return job

    def create_final_report(
        self,
        output_dir: Path,
        mutational_loads: List[MutationalLoad],
        qc_reports: Dict[str, SampleQCReport],
        sample_pairs: Dict[str, str],
        project_name: str = "Exome_Analysis",
    ) -> Dict[str, ppg.Job]:
        """
        Create comprehensive final report with multiple outputs.
        
        Args:
            output_dir: Output directory
            mutational_loads: List of mutational load results
            qc_reports: QC reports for all samples
            sample_pairs: Tumor-normal pairs mapping
            project_name: Project name for report
            
        Returns:
            Dictionary of pypipegraph jobs
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        jobs = {}
        
        # Main mutational load report
        ml_report = output_dir / f"{project_name}_mutational_load.tsv"
        jobs['mutational_load'] = self.generate_mutational_load_report(
            mutational_loads, ml_report, qc_reports
        )
        
        # Summary statistics
        summary_file = output_dir / f"{project_name}_summary.txt"
        jobs['summary'] = self._generate_summary_stats(
            mutational_loads, qc_reports, summary_file
        )
        
        # Individual sample reports
        for sample_name, qc_report in qc_reports.items():
            # Find corresponding mutational load if tumor
            ml = None
            for ml_obj in mutational_loads:
                if ml_obj.sample_name == sample_name:
                    ml = ml_obj
                    break
            
            sample_report = output_dir / "samples" / f"{sample_name}_report.json"
            job_name = f'sample_report_{sample_name}'
            jobs[job_name] = self.generate_sample_summary_report(
                sample_name, qc_report, ml, sample_report
            )
        
        return jobs

    def _generate_summary_stats(
        self,
        mutational_loads: List[MutationalLoad],
        qc_reports: Dict[str, SampleQCReport],
        output_file: Path,
    ) -> ppg.FileGeneratingJob:
        """Generate summary statistics text file."""
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def generate():
            with open(output_file, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("EXOME SEQUENCING ANALYSIS SUMMARY\n")
                f.write("=" * 60 + "\n\n")
                
                # Sample counts
                f.write(f"Total samples analyzed: {len(qc_reports)}\n")
                f.write(f"Tumor samples: {len(mutational_loads)}\n\n")
                
                # Mutational load statistics
                if mutational_loads:
                    loads = [ml.mutational_load for ml in mutational_loads]
                    f.write("Mutational Load Statistics:\n")
                    f.write("-" * 40 + "\n")
                    f.write(f"  Mean: {sum(loads) / len(loads):.2f} mutations/Mb\n")
                    f.write(f"  Min:  {min(loads):.2f} mutations/Mb\n")
                    f.write(f"  Max:  {max(loads):.2f} mutations/Mb\n\n")
                
                # QC summary
                f.write("QC Summary:\n")
                f.write("-" * 40 + "\n")
                
                mapping_rates = []
                duplication_rates = []
                mean_coverages = []
                
                for sample_name, qc in qc_reports.items():
                    if qc.alignment:
                        mapping_rates.append(qc.alignment.mapping_rate)
                        duplication_rates.append(qc.alignment.duplication_rate)
                        mean_coverages.append(qc.alignment.mean_coverage)
                
                if mapping_rates:
                    f.write(f"  Mean mapping rate: {sum(mapping_rates) / len(mapping_rates):.1f}%\n")
                if duplication_rates:
                    f.write(f"  Mean duplication rate: {sum(duplication_rates) / len(duplication_rates):.1f}%\n")
                if mean_coverages:
                    f.write(f"  Mean coverage: {sum(mean_coverages) / len(mean_coverages):.1f}x\n")
                
                f.write("\n" + "=" * 60 + "\n")

        job = ppg.FileGeneratingJob(output_file, generate)
        return job

    def plot_mutational_load(
        self,
        mutational_loads: List[MutationalLoad],
        output_file: Path,
        group_by: Optional[str] = None,
    ) -> ppg.FileGeneratingJob:
        """
        Create mutational load plot.
        
        Args:
            mutational_loads: List of mutational loads
            output_file: Output plot file (PDF or PNG)
            group_by: Optional grouping variable (e.g., tissue type)
            
        Returns:
            pypipegraph job
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        def plot():
            try:
                import matplotlib.pyplot as plt
                import seaborn as sns
                
                # Prepare data
                samples = [ml.sample_name for ml in mutational_loads]
                loads = [ml.mutational_load for ml in mutational_loads]
                
                # Create plot
                plt.figure(figsize=(10, 6))
                plt.bar(range(len(samples)), loads)
                plt.xticks(range(len(samples)), samples, rotation=45, ha='right')
                plt.ylabel('Mutational Load (mutations/Mb)')
                plt.xlabel('Sample')
                plt.title('Tumor Mutational Burden')
                plt.tight_layout()
                plt.savefig(output_file, dpi=300)
                plt.close()
                
            except ImportError:
                # If matplotlib not available, create simple text output
                with open(output_file.with_suffix('.txt'), 'w') as f:
                    f.write("Mutational Load Plot (matplotlib not available)\n\n")
                    for ml in mutational_loads:
                        f.write(f"{ml.sample_name}: {ml.mutational_load:.2f} mutations/Mb\n")

        job = ppg.FileGeneratingJob(output_file, plot)
        return job
