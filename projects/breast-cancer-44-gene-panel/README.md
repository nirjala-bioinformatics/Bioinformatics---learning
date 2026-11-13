
# Breast Cancer 44-Gene Targeted NGS Analysis

## Project Overview

This project focuses on the bioinformatics analysis of targeted Next-Generation Sequencing (NGS) data from a 44-gene breast cancer susceptibility panel.

The objective was to perform quality control, read preprocessing, genome alignment, variant calling, target-region coverage assessment, variant annotation, and comparison of identified variants with findings reported in a published study.

## Dataset Information

| Parameter | Details |
|---|---|
| BioProject | PRJNA820693 |
| Study | SRP366275 |
| Reference Genome | GRCh38 |
| Gene Panel | 44 breast cancer-associated genes |
| Samples Analysed | 5 selected SRA runs |
| Total Target Region | 202,857 bases |

### SRA Accessions

- SRR18502497
- SRR18502498
- SRR18502508
- SRR18502509
- SRR18502511

## Bioinformatics Workflow

1. **Quality Control:** FastQC was used to assess raw sequencing read quality.
2. **Read Trimming:** fastp was used for adapter trimming and quality filtering.
3. **Quality Report Aggregation:** MultiQC was used to summarize quality-control reports.
4. **Read Alignment:** BWA-MEM was used to align reads against the GRCh38 reference genome.
5. **BAM Processing:** SAMtools was used for sorting and indexing alignment files.
6. **Duplicate Marking:** Picard MarkDuplicates was used to identify duplicate reads.
7. **Base Quality Score Recalibration:** GATK BaseRecalibrator and ApplyBQSR were used.
8. **Variant Calling:** GATK HaplotypeCaller was used for per-sample variant calling.
9. **Variant Filtering:** SNP and INDEL variants were processed using GATK filtering steps.
10. **Target Coverage Analysis:** BED-based coverage analysis was performed across the 44-gene panel.
11. **Variant Annotation:** Ensembl VEP and Jannovar were used for variant annotation.
12. **Variant Review:** Franklin-compatible VCF files were prepared for downstream variant review.

## Target Region Coverage and Depth Analysis

Coverage was evaluated across the five selected samples.

| SRA Accession | Mean Depth | 1× Coverage | 10× Coverage | 20× Coverage | 30× Coverage |
|---|---:|---:|---:|---:|---:|
| SRR18502497 | 18.05× | 97.32% | 62.29% | 38.00% | 19.20% |
| SRR18502498 | 9.41× | 97.27% | 41.55% | 5.52% | 0.79% |
| SRR18502508 | 12.16× | 97.22% | 53.94% | 17.38% | 4.27% |
| SRR18502509 | 10.16× | 97.35% | 46.05% | 8.51% | 1.36% |
| SRR18502511 | 13.18× | 97.32% | 56.30% | 21.23% | 6.61% |

**Coverage interpretation:** 1× coverage represents the percentage of target bases covered by at least one sequencing read. Higher depth thresholds provide additional information about sequencing coverage. Coverage at 1× alone does not establish clinical-grade sequencing adequacy.

## Variant Annotation and Franklin Review

Variant annotation was performed using:

- Ensembl Variant Effect Predictor (VEP)
- Jannovar

Separate VEP annotation summary reports were generated for all five samples.

Franklin-compatible VCF files were prepared for downstream variant review.

### Exact Variant Match with the Published Study

During Franklin-based variant review, one exact match was identified with the variant reported in the reference publication:

| Parameter | Variant Information |
|---|---|
| Gene | MUTYH |
| Nucleotide Change | c.1103G>A |
| Protein Change | p.(Gly368Asp) |
| dbSNP ID | rs36053993 |

This variant was reported in the reference study (P21/P39, Figure 6).

The finding represents an exact match with the published variant record. It should not be interpreted as independent clinical validation or a diagnostic conclusion.

## Tools and Technologies

- Linux / Ubuntu (WSL)
- Bash scripting
- FastQC
- fastp
- MultiQC
- BWA-MEM
- SAMtools
- Picard
