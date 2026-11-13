#!/usr/bin/env bash
set -Eeuo pipefail

# ==========================================
# PROJECT CONFIGURATION
# ==========================================

PROJECT="$HOME/gene_panel_5samples"
FASTQ="$PROJECT/fastq"
REF="$HOME/NGS/reference/GRCh38.fa"
MILLS="$HOME/NGS/known_sites/Mills_and_1000G_gold_standard.indels.hg38.NC.main.vcf.gz"

THREADS=4

SAMPLES=(
    SRR18502498
    SRR18502508
    SRR18502509
    SRR18502511
)

ALL_SAMPLES=(
    SRR18502497
    SRR18502498
    SRR18502508
    SRR18502509
    SRR18502511
)

cd "$PROJECT"

# ==========================================
# CHECK REQUIRED TOOLS
# ==========================================

for TOOL in fastqc fastp multiqc bwa samtools gatk bcftools; do
    command -v "$TOOL" >/dev/null 2>&1 || {
        echo "ERROR: $TOOL not found"
        exit 1
    }
done

[[ -f "$REF" ]] || { echo "Reference missing"; exit 1; }
[[ -f "$MILLS" ]] || { echo "Mills known-sites VCF missing"; exit 1; }

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -f "$FASTQ/${SAMPLE}_1.fastq" ]] || {
        echo "Missing R1: $SAMPLE"
        exit 1
    }
    [[ -f "$FASTQ/${SAMPLE}_2.fastq" ]] || {
        echo "Missing R2: $SAMPLE"
        exit 1
    }
done

# ==========================================
# PREPARE REFERENCE
# ==========================================

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

DICT="${REF%.*}.dict"
[[ -f "$DICT" ]] || \
    gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

[[ -f "${REF}.bwt" ]] || bwa index "$REF"

# ==========================================
# PROCESS FOUR REMAINING SAMPLES
# ==========================================

for SAMPLE in "${SAMPLES[@]}"; do

    echo "======================================"
    echo "PROCESSING SAMPLE: $SAMPLE"
    echo "======================================"

    R1="$FASTQ/${SAMPLE}_1.fastq"
    R2="$FASTQ/${SAMPLE}_2.fastq"

    TRIM1="trimmed/${SAMPLE}_R1.trim.fastq.gz"
    TRIM2="trimmed/${SAMPLE}_R2.trim.fastq.gz"

    SAM="alignment/${SAMPLE}.sam"
    UNSORTED="bam/${SAMPLE}_unsorted.bam"
    SORTED="bam/${SAMPLE}_sorted.bam"
    RG_BAM="bam/${SAMPLE}_RG.bam"
    MARKED="bam/${SAMPLE}_marked.bam"
    RECAL="bqsr/${SAMPLE}_recal.bam"

    # Remove previous intermediate outputs for this sample
    rm -f "$SAM" "$UNSORTED" "$SORTED" "$SORTED.bai"
    rm -f "$RG_BAM" "$RG_BAM.bai"
    rm -f "$MARKED" "$MARKED.bai"
    rm -f "$RECAL" "$RECAL.bai"
    rm -f "output/${SAMPLE}"*.vcf.gz*
    rm -f "bqsr/${SAMPLE}_recal.table"

    # STEP 1: RAW FASTQC

    fastqc -t "$THREADS" \
        -o qc/raw "$R1" "$R2"

    # STEP 2: FASTP

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$TRIM1" \
        -O "$TRIM2" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --html "fastp_reports/${SAMPLE}.html" \
        --json "fastp_reports/${SAMPLE}.json"

    # STEP 3: TRIMMED FASTQC

    fastqc -t "$THREADS" \
        -o qc/trimmed "$TRIM1" "$TRIM2"

    # STEP 4: BWA ALIGNMENT

    bwa mem -t "$THREADS" "$REF" "$TRIM1" "$TRIM2" \
        > "$SAM"

    # STEP 5: SAM TO BAM

    samtools view -@ "$THREADS" -b \
        -o "$UNSORTED" "$SAM"

    # STEP 6: SORT BAM

    samtools sort -@ "$THREADS" \
        -o "$SORTED" "$UNSORTED"

    # STEP 7: INDEX BAM

    samtools index "$SORTED"

    # STEP 8: ADD READ GROUPS

    gatk AddOrReplaceReadGroups \
        -I "$SORTED" \
        -O "$RG_BAM" \
        -RGID "$SAMPLE" \
        -RGLB "$SAMPLE" \
        -RGPL ILLUMINA \
        -RGPU "${SAMPLE}.unit1" \
        -RGSM "$SAMPLE"

    samtools index "$RG_BAM"

    # STEP 9: MARK DUPLICATES

    gatk MarkDuplicates \
        -I "$RG_BAM" \
        -O "$MARKED" \
        -M "logs/${SAMPLE}_duplicate_metrics.txt"

    samtools index "$MARKED"

    # STEP 10: BASE RECALIBRATOR

    gatk BaseRecalibrator \
        -R "$REF" \
        -I "$MARKED" \
        --known-sites "$MILLS" \
        -O "bqsr/${SAMPLE}_recal.table"

    # STEP 11: APPLY BQSR

    gatk ApplyBQSR \
        -R "$REF" \
        -I "$MARKED" \
        --bqsr-recal-file "bqsr/${SAMPLE}_recal.table" \
        -O "$RECAL"

    samtools index "$RECAL"

    # STEP 12: HAPLOTYPE CALLER
    # Individual sample VCF; no joint genotyping

    gatk HaplotypeCaller \
        -R "$REF" \
        -I "$RECAL" \
        -O "output/${SAMPLE}_raw.vcf.gz"

    # STEP 13: SELECT SNPs

    gatk SelectVariants \
        -R "$REF" \
        -V "output/${SAMPLE}_raw.vcf.gz" \
        --select-type-to-include SNP \
        -O "output/${SAMPLE}_SNP.vcf.gz"

    # STEP 14: SELECT INDELs

    gatk SelectVariants \
        -R "$REF" \
        -V "output/${SAMPLE}_raw.vcf.gz" \
        --select-type-to-include INDEL \
        -O "output/${SAMPLE}_INDEL.vcf.gz"

    # STEP 15: SNP FILTRATION

    gatk VariantFiltration \
        -R "$REF" \
        -V "output/${SAMPLE}_SNP.vcf.gz" \
        --filter-name "SNP_FILTER" \
        --filter-expression \
        "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        -O "output/${SAMPLE}_SNP_filtered.vcf.gz"

    # STEP 16: INDEL FILTRATION

    gatk VariantFiltration \
        -R "$REF" \
        -V "output/${SAMPLE}_INDEL.vcf.gz" \
        --filter-name "INDEL_FILTER" \
        --filter-expression \
        "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        -O "output/${SAMPLE}_INDEL_filtered.vcf.gz"

    # STEP 17: MERGE SNP AND INDEL

    gatk MergeVcfs \
        -I "output/${SAMPLE}_SNP_filtered.vcf.gz" \
        -I "output/${SAMPLE}_INDEL_filtered.vcf.gz" \
        -O "output/${SAMPLE}_filtered.vcf.gz"

    # STEP 18: VERIFY FILTERED VCF

    bcftools norm \
        -f "$REF" \
        -c e \
        "output/${SAMPLE}_filtered.vcf.gz" \
        -Ou > /dev/null

    [[ -s "output/${SAMPLE}_filtered.vcf.gz" ]] || {
        echo "ERROR: Filtered VCF missing for $SAMPLE"
        exit 1
    }

    echo "SUCCESS: $SAMPLE filtered VCF completed"

done

# ==========================================
# MULTIQC REPORT FOR ALL SAMPLES
# ==========================================

multiqc \
    qc/raw \
    qc/trimmed \
    fastp_reports \
    -o multiqc_reports \
    --force

# ==========================================
# FRANKLIN CONVERSION FOR ALL FIVE SAMPLES
# ==========================================

cat > output/chr_mapping.txt <<'MAPEOF'
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
NC_000004.12 chr4
NC_000005.10 chr5
NC_000006.12 chr6
NC_000007.14 chr7
NC_000008.11 chr8
NC_000009.12 chr9
NC_000010.11 chr10
NC_000011.10 chr11
NC_000012.12 chr12
NC_000013.11 chr13
NC_000014.9 chr14
NC_000015.10 chr15
NC_000016.10 chr16
NC_000017.11 chr17
NC_000018.10 chr18
NC_000019.10 chr19
NC_000020.11 chr20
NC_000021.9 chr21
NC_000022.11 chr22
NC_000023.11 chrX
NC_000024.10 chrY
NC_012920.1 chrM
MAPEOF

echo '##reference=GRCh38' > output/reference_header.txt

STANDARD_CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

for SAMPLE in "${ALL_SAMPLES[@]}"; do

    echo "Franklin conversion: $SAMPLE"

    INPUT="output/${SAMPLE}_filtered.vcf.gz"

    [[ -s "$INPUT" ]] || {
        echo "ERROR: Missing filtered VCF for $SAMPLE"
        exit 1
    }

    # STEP 1: Rename chromosome accessions

    bcftools annotate \
        --rename-chrs output/chr_mapping.txt \
        -Ov \
        -o "output/${SAMPLE}_chr.vcf" \
        "$INPUT"

    # STEP 2: Keep standard chromosomes only

    bcftools view \
        -t "$STANDARD_CHROMS" \
        -Ov \
        -o "output/${SAMPLE}_chr_standard.vcf" \
        "output/${SAMPLE}_chr.vcf"

    # STEP 3: Add GRCh38 reference information

    bcftools annotate \
        -h output/reference_header.txt \
        -Ov \
        -o "output/${SAMPLE}_Franklin.vcf" \
        "output/${SAMPLE}_chr_standard.vcf"

    [[ -s "output/${SAMPLE}_Franklin.vcf" ]] || {
        echo "ERROR: Franklin VCF missing for $SAMPLE"
        exit 1
    }

    echo "Franklin VCF created: ${SAMPLE}_Franklin.vcf"

done

echo "======================================"
echo "ALL FIVE SAMPLES PROCESSED"
echo "Franklin VCFs are in output/"
echo "======================================"

