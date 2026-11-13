# NGS/WGS Quality Control Analysis


## Project Overview

This project demonstrates a practical bioinformatics workflow for quality control and preprocessing of paired-end sequencing reads from sample ERR1485610.

## Tools Used

- **FastQC:** Initial sequencing quality assessment
- **fastp:** Quality filtering and adapter trimming
- **MultiQC:** Combined quality-control report

## Workflow

1. FastQC analysis of paired-end reads
2. Read filtering and trimming using fastp
3. Post-trimming quality assessment using FastQC and MultiQC

## Reports

- [FastQC Report — Read 1](ERR1485610_1_fastqc.html)
- [FastQC Report — Read 2](ERR1485610_2_fastqc.html)
- [fastp Report](ERR1485610_fastp.html)
- [MultiQC Report — ERR1485610 (Single Sample)](ERR1485610_multiqc_report.html)

# NGS/WGS Quality Control Analysis
## Results Summary

| Metric | Result |
|---|---:|
| Sample | ERR1485610 |
| Reads after filtering | ~0.9 million |
| GC content | 38.0% |
| Reads passing filter (% PF) | 90.3% |
| Adapter content | 1.4% |
| Read length | 126 bp |

The FastQC and fastp reports were consolidated using MultiQC for quality assessment. These results summarize the preprocessing stage before downstream alignment and variant analysis.

### Alignment Summary

Alignment was performed using BWA-MEM against the GRCh38 reference genome. Alignment statistics were generated using `samtools flagstat` from the recalibrated BAM file.

- **Total reads:** 895,980
- **Mapped reads:** 99.99%
- **Properly paired reads:** 99.22%
- **Duplicate reads:** 33,357

**Output file:** [`alignment_summary.txt`](alignment_summary.txt)



## Project Objective

To gain practical experience in NGS data quality control and preprocessing using public sequencing data.

## Project Status

Learning and portfolio project.
## Downstream Analysis

### Reference Genome and Alignment

- Reference genome: GRCh38
- Alignment tool: BWA-MEM
- BAM processing: SAMtools and Picard
- Base quality score recalibration: GATK

### Coverage and Depth Summary

| Metric | Result |
|---|---:|
| Mean depth | 0.0335207× |
| Coverage at 1× | 0.1318% |
| Coverage at 10× | 0.0544% |
| Coverage at 20× | 0.0410% |
| Coverage at 30× | 0.0325% |

See [coverage_depth_summary.txt](coverage_depth_summary.txt).

### Variant Calling

Variant calling was performed using GATK HaplotypeCaller and GenotypeGVCFs.

| Metric | Result |
|---|---:|
| PASS variant records | 3,211 |
| SNPs | 2,476 |
| INDELs | 735 |
| Multiallelic sites | 64 |

See [variant_calling_summary.txt](variant_calling_summary.txt).

### Functional Annotation

Variants were annotated using Ensembl Variant Effect Predictor (VEP) version 116.3.

- Input PASS records: 3,211
- Annotated records: 3,210
- HIGH-impact candidate records: 2
- MODERATE-impact candidate records: 12

See [vep_annotation_summary.txt](vep_annotation_summary.txt).

### Candidate Variants

Two HIGH-impact candidate variants were identified:

- **HRH2:** Predicted splice acceptor variant.
- **TROAP:** Predicted frameshift variant.

See [candidate_variants.txt](candidate_variants.txt) for details.

## Limitations

The sample showed extremely low mean sequencing depth (approximately 0.0335×) and low coverage across positions assessed in the BAM reference contigs.

Therefore, variant calls and predicted consequences should be treated as preliminary computational findings. They are not clinically validated or confirmed pathogenic variants.

## Conclusion

This project provided practical experience in a human WGS bioinformatics workflow, including quality control, read alignment, BAM processing, variant calling, functional annotation, and interpretation of sequencing quality limitations.

**Project type:** Independent secondary analysis of public WGS data for learning and portfolio development.
