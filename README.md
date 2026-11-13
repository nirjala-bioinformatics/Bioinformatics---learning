# Breast Cancer 44-Gene Targeted NGS Analysis

## Project Overview

This project focuses on the bioinformatics analysis of targeted Next-Generation Sequencing (NGS) data from a 44-gene breast cancer susceptibility panel.

The objective was to process sequencing reads, assess sequence quality, perform alignment and variant calling, evaluate target-region coverage, and annotate genetic variants.

## Dataset

- BioProject: PRJNA820693
- Study: SRP366275
- Reference genome: GRCh38
- Gene panel: 44 breast cancer-associated genes
- Samples analysed: 5 SRA runs
- Total target regions: 202,857 bases

### SRA Accessions

| Sample | SRA Accession |
|---|---|
| Sample 1 | SRR18502497 |
| Sample 2 | SRR18502498 |
| Sample 3 | SRR18502508 |
| Sample 4 | SRR18502509 |
| Sample 5 | SRR18502511 |

## Bioinformatics Workflow

1. **Quality Control:** FastQC was used to assess raw sequencing read quality.
2. **Read Trimming:** fastp was used for adapter trimming and quality filtering.
3. **Quality Assessment:** MultiQC was used to aggregate quality-control reports.
4. **Read Alignment:** BWA-MEM was used to align reads against the GRCh38 reference genome.
5. **BAM Processing:** SAMtools was used for sorting and indexing alignment files.
6. **Duplicate Marking:** Picard MarkDuplicates was used to identify duplicate reads.
7. **Base Quality Score Recalibration:** GATK BaseRecalibrator and ApplyBQSR were used.
8. **Variant Calling:** GATK HaplotypeCaller was used for per-sample variant calling.
9. **Variant Filtering:** SNP and INDEL variants were processed using GATK filtering steps.
10. **Target Region Analysis:** BEDTools and coverage analysis were used to assess coverage across the 44-gene panel.
11. **Variant Annotation:** Ensembl VEP and Jannovar were used for variant annotation.
12. **Variant Review:** Franklin-compatible VCF files were prepared for downstream variant review.

## Coverage and Depth Analysis

Target-region coverage was evaluated across five samples.

| SRA Accession | Mean Depth | 1x Coverage | 10x Coverage | 20x Coverage | 30x Coverage |
|---|---:|---:|---:|---:|---:|
| SRR18502497 | 18.05x | 97.32% | 62.29% | 38.00% | 19.20% |
| SRR18502498 | 9.41x | 97.27% | 41.55% | 5.52% | 0.79% |
| SRR18502508 | 12.16x | 97.22% | 53.94% | 17.38% | 4.27% |
| SRR18502509 | 10.16x | 97.35% | 46.05% | 8.51% | 1.36% |
| SRR18502511 | 13.18x | 97.32% | 56.30% | 21.23% | 6.61% |

Coverage at 1x represents the percentage of target bases covered by at least one sequencing read. Higher depth thresholds provide additional information about sequencing coverage adequacy.

## Variant Annotation

Variant annotation was performed using:

- Ensembl Variant Effect Predictor (VEP)
- Jannovar

Separate annotation reports were generated for the five samples.

## Tools and Technologies

- Linux / Ubuntu (WSL)
- Bash scripting
- FastQC
- fastp
- MultiQC
- BWA-MEM
- SAMtools
- Picard
- GATK
- BEDTools
- Ensembl VEP
- Jannovar
- Franklin-compatible VCF preparation

## Repository Structure

```text
breast-cancer-44-gene-panel/
|
|-- scripts/
|   |-- fastq_to_vcf.sh
|   |-- run_remaining_4_samples.sh
|   |-- run_jannovar_5samples.sh
|   |-- extract_jannovar_separate_annotations.py
|
|-- reference/
|   |-- genes_44.txt
|   |-- breast_cancer_44genes_GRCh38_50bp_accession.bed
|   |-- NC_to_chr.txt
|
|-- results/
|   |-- all_samples_coverage.tsv
|   |-- fastp_summary.txt
|
|-- qc_reports/
|   |-- multiqc_report.html
|   |-- sample_fastp_reports/
|
|-- annotation_reports/
|   |-- VEP annotation summary reports
|
|-- README.md
```

## Learning Outcomes

- Practical experience with targeted NGS data analysis.
- Understanding of FASTQ quality control and read trimming.
- Read alignment and BAM processing.
- Duplicate marking and base quality score recalibration.
- Variant calling and SNP/INDEL filtering.
- Target-region coverage and sequencing depth assessment.
- Variant annotation using VEP and Jannovar.
- Bash scripting and Linux command-line workflows.

## Limitations

- Analysis was performed on five selected sequencing runs.
- Coverage varied between samples and across depth thresholds.
- Coverage at 1x alone does not establish clinical-grade sequencing adequacy.
- Variant annotation does not constitute clinical interpretation or diagnostic validation.

## Data Source

NCBI Sequence Read Archive (SRA)

BioProject: PRJNA820693

Study: SRP366275

Reference genome: GRCh38

---

**Project Type:** Targeted Next-Generation Sequencing (NGS) Analysis  
**Focus:** Breast Cancer Susceptibility Gene Panel
