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

## Project Overview

This project demonstrates a practical bioinformatics workflow for quality control and preprocessing of paired-end sequencing reads from sample ERR1485610.

## Tools Used

- **FastQC:** Initial sequencing quality assessment
- **fastp:** Quality filtering and adapter trimming

## Project Objective

To gain practical experience in NGS data quality control and preprocessing using public sequencing data.

## Project Status

Learning and portfolio project.
