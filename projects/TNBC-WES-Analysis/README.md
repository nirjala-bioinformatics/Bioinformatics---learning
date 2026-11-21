# TNBC Whole Exome Sequencing (WES) Analysis

## Project Overview

This project focuses on tumour-only Whole Exome Sequencing (WES) data analysis of a Triple-Negative Breast Cancer (TNBC) sample using bioinformatics tools and the GATK best-practice-based workflow.

## Dataset Information

- **BioProject:** PRJNA1422845
- **SRA Study:** SRP676624
- **SRA Run:** SRR37211335
- **Sample:** WES_Tumor_01
- **Reference Genome:** GRCh38

## Bioinformatics Workflow

1. FASTQ Quality Control – FastQC
2. Read Preprocessing – fastp
3. Quality Report Aggregation – MultiQC
4. Read Alignment – BWA-MEM
5. BAM Processing – SAMtools
6. Duplicate Marking – Picard
7. Base Quality Score Recalibration – GATK BQSR
8. Somatic Variant Calling – GATK Mutect2
9. Contamination Estimation – GetPileupSummaries and CalculateContamination
10. Variant Filtering – FilterMutectCalls
11. Functional Annotation – Ensembl VEP
12. Candidate Variant Review – IGV and Franklin
## Workflow Diagram

![TNBC WES Analysis Workflow](tnbc-wes-workflow.png)

## Key Results

| Parameter | Result |
|---|---:|
| Total sequencing reads | 44,320,552 |
| Mapping rate | 99.98% |
| Duplicate rate | 6.92% |
| Mean target depth | 54.41× |
| Target bases covered at ≥1× | 99.15% |
| PASS variant records | 983 |

Six candidate variants were prioritised for further review involving KMT2C, TP53, BRIP1, GNAS, ACADVL and CENPF.

## Tools and Technologies

FastQC, fastp, MultiQC, BWA-MEM, SAMtools, Picard, GATK, Ensembl VEP, IGV, Franklin, Linux and Bash.

## Important Limitation

This analysis was performed using tumour-only sequencing data because matched normal sequencing data was unavailable. Therefore, the candidate variants cannot be confirmed as somatic variants without matched normal comparison and further validation.

