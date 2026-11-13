
#!/usr/bin/env bash
set -Eeuo pipefail

# ==========================================
# 1. PROJECT CONFIGURATION
# ==========================================

PROJECT="$HOME/gene_panel_5samples"
FASTQ="$PROJECT/fastq"
REF="$HOME/NGS/reference/GRCh38.fa"
KNOWN_DIR="$HOME/NGS/known_sites"

MILLS="$KNOWN_DIR/Mills_and_1000G_gold_standard.indels.hg38.NC.main.vcf.gz"

THREADS=4

SAMPLES=(
    SRR18502497
    SRR18502498
    SRR18502508
    SRR18502509
    SRR18502511
)

cd "$PROJECT"

# ==========================================
# 2. CHECK REQUIRED TOOLS AND FILES
# ==========================================

for TOOL in fastqc fastp multiqc bwa samtools gatk bcftools tabix; do
    command -v "$TOOL" >/dev/null 2>&1 || {
        echo "ERROR: $TOOL is not installed or not in PATH"
        exit 1
    }
done

[[ -f "$REF" ]] || {
    echo "ERROR: Reference FASTA not found: $REF"
    exit 1
}

[[ -f "$MILLS" ]] || {
    echo "ERROR: Mills known-sites VCF not found: $MILLS"
    exit 1
}

for SAMPLE in "${SAMPLES[@]}"; do
    [[ -f "$FASTQ/${SAMPLE}_1.fastq" ]] || {
        echo "ERROR: Missing R1 FASTQ for $SAMPLE"
        exit 1
    }

    [[ -f "$FASTQ/${SAMPLE}_2.fastq" ]] || {
        echo "ERROR: Missing R2 FASTQ for $SAMPLE"
        exit 1
    }
done

# ==========================================
# 3. PREPARE REFERENCE
# ==========================================

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

DICT="${REF%.*}.dict"

[[ -f "$DICT" ]] || \
    gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

[[ -f "${REF}.bwt" ]] || bwa index "$REF"

# ==========================================
# 4. FRANKLIN CONTIG MAPPING
# ==========================================

MAP="$PROJECT/reference/NC_to_chr.txt"

cat > "$MAP" <<'MAPEOF'
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

# ==========================================
# 5. PROCESS EACH SAMPLE INDEPENDENTLY
# ==========================================

for SAMPLE in "${SAMPLES[@]}"; do

    echo "=========================================="
    echo "Starting sample: $SAMPLE"
    echo "=========================================="

    R1="$FASTQ/${SAMPLE}_1.fastq"
    R2="$FASTQ/${SAMPLE}_2.fastq"

    TRIM1="$PROJECT/trimmed/${SAMPLE}_R1.trim.fastq.gz"
    TRIM2="$PROJECT/trimmed/${SAMPLE}_R2.trim.fastq.gz"

    SAM="$PROJECT/alignment/${SAMPLE}.sam"
    UNSORTED="$PROJECT/bam/${SAMPLE}_unsorted.bam"
    SORTED="$PROJECT/bam/${SAMPLE}_sorted.bam"
    RG_BAM="$PROJECT/bam/${SAMPLE}_RG.bam"
    MARKED="$PROJECT/bam/${SAMPLE}_marked.bam"
    RECAL="$PROJECT/bqsr/${SAMPLE}_recal.bam"

    # --------------------------------------
    # STEP 1: RAW FASTQC
    # --------------------------------------

    fastqc -t "$THREADS" \
        -o qc/raw "$R1" "$R2"

    # --------------------------------------
    # STEP 2: FASTP
    # --------------------------------------

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$TRIM1" \
        -O "$TRIM2" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --html "fastp_reports/${SAMPLE}.html" \
        --json "fastp_reports/${SAMPLE}.json"

    # --------------------------------------
    # STEP 3: TRIMMED FASTQC
    # --------------------------------------

    fastqc -t "$THREADS" \
        -o qc/trimmed "$TRIM1" "$TRIM2"

    # --------------------------------------
    # STEP 4: BWA ALIGNMENT
    # --------------------------------------

    bwa mem -t "$THREADS" "$REF" "$TRIM1" "$TRIM2" \
        > "$SAM"

    # --------------------------------------
    # STEP 5: SAM TO BAM
    # --------------------------------------

    samtools view -@ "$THREADS" -b \
        -o "$UNSORTED" "$SAM"

    # --------------------------------------
    # STEP 6: SORT BAM
    # --------------------------------------

    samtools sort -@ "$THREADS" \
        -o "$SORTED" "$UNSORTED"

    # --------------------------------------
    # STEP 7: INDEX SORTED BAM
    # --------------------------------------


    samtools index "$SORTED"

    # --------------------------------------
    # STEP 8: ADD READ GROUPS
    # --------------------------------------

    gatk AddOrReplaceReadGroups \
        -I "$SORTED" \
        -O "$RG_BAM" \
        -RGID "$SAMPLE" \
        -RGLB "$SAMPLE" \
        -RGPL ILLUMINA \
        -RGPU "${SAMPLE}.unit1" \
        -RGSM "$SAMPLE"

    samtools index "$RG_BAM"

    # --------------------------------------
    # STEP 9: MARK DUPLICATES
    # --------------------------------------

    gatk MarkDuplicates \
        -I "$RG_BAM" \
        -O "$MARKED" \
        -M "logs/${SAMPLE}_duplicate_metrics.txt"

    samtools index "$MARKED"

    # --------------------------------------
    # STEP 10: BASE RECALIBRATOR
    # --------------------------------------

    gatk BaseRecalibrator \
        -R "$REF" \
        -I "$MARKED" \
        --known-sites "$MILLS" \
        -O "bqsr/${SAMPLE}_recal.table"

    # --------------------------------------
    # STEP 11: APPLY BQSR
    # --------------------------------------

    gatk ApplyBQSR \
        -R "$REF" \
        -I "$MARKED" \
        --bqsr-recal-file "bqsr/${SAMPLE}_recal.table" \
        -O "$RECAL"

    samtools index "$RECAL"

    # --------------------------------------
    # STEP 12: HAPLOTYPE CALLER
    # Individual sample VCF; no joint genotyping
    # --------------------------------------

    gatk HaplotypeCaller \
        -R "$REF" \
        -I "$RECAL" \
        -O "output/${SAMPLE}_raw.vcf.gz"

    # --------------------------------------
    # STEP 13: SELECT SNPs AND INDELS
    # --------------------------------------

    gatk SelectVariants \
        -R "$REF" \
        -V "output/${SAMPLE}_raw.vcf.gz" \
        --select-type-to-include SNP \
        -O "output/${SAMPLE}_SNP.vcf.gz"

    gatk SelectVariants \
        -R "$REF" \
        -V "output/${SAMPLE}_raw.vcf.gz" \
        --select-type-to-include INDEL \
        -O "output/${SAMPLE}_INDEL.vcf.gz"

    # --------------------------------------
    # STEP 14: SNP FILTRATION
    # --------------------------------------

    gatk VariantFiltration \
        -R "$REF" \
        -V "output/${SAMPLE}_SNP.vcf.gz" \
        --filter-name "SNP_FILTER" \
        --filter-expression \
        "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        -O "output/${SAMPLE}_SNP_filtered.vcf.gz"

    # --------------------------------------
    # STEP 15: INDEL FILTRATION
    # --------------------------------------

    gatk VariantFiltration \
        -R "$REF" \
        -V "output/${SAMPLE}_INDEL.vcf.gz" \
        --filter-name "INDEL_FILTER" \
        --filter-expression \
        "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        -O "output/${SAMPLE}_INDEL_filtered.vcf.gz"

    # --------------------------------------
    # STEP 16: MERGE FILTERED SNP AND INDEL
    # --------------------------------------

    gatk MergeVcfs \
        -I "output/${SAMPLE}_SNP_filtered.vcf.gz" \
        -I "output/${SAMPLE}_INDEL_filtered.vcf.gz" \
        -O "output/${SAMPLE}_filtered.vcf.gz"

    # --------------------------------------
    # STEP 17: VALIDATE ORIGINAL FILTERED VCF
    # --------------------------------------

    bcftools norm \
        -f "$REF" \
        -c e \
        "output/${SAMPLE}_filtered.vcf.gz" \
        -Ou > /dev/null

    # --------------------------------------
    # STEP 18: REMOVE NW AND NT CONTIGS
    # AND CONVERT NC ACCESSIONS TO chr FORMAT
    # --------------------------------------

    EXCLUDE="output/${SAMPLE}_exclude_contigs.txt"

    bcftools view -h "output/${SAMPLE}_filtered.vcf.gz" |
        sed -n \
        -e 's/^##contig=<ID=\(NW_[^,>]*\).*/\1/p' \
        -e 's/^##contig=<ID=\(NT_[^,>]*\).*/\1/p' \
        > "$EXCLUDE"

    FRANKLIN_TMP="output/${SAMPLE}_Franklin_tmp.vcf.gz"
    FRANKLIN="output/${SAMPLE}_Franklin.vcf.gz"

    if [[ -s "$EXCLUDE" ]]; then
        bcftools view -T "^$EXCLUDE" -Ou \
            "output/${SAMPLE}_filtered.vcf.gz" |
        bcftools annotate \
            --rename-chrs "$MAP" \
            -Oz -o "$FRANKLIN_TMP"
    else
        bcftools annotate \
            --rename-chrs "$MAP" \
            -Oz -o "$FRANKLIN_TMP" \
            "output/${SAMPLE}_filtered.vcf.gz"
    fi

    # Remove NW/NT contigs from the VCF header too
    bcftools view -h "$FRANKLIN_TMP" |
        sed -E '/^##contig=<ID=(NW_|NT_)/d' \
        > "output/${SAMPLE}_Franklin_header.txt"

    bcftools reheader \
        -h "output/${SAMPLE}_Franklin_header.txt" \
        -o "$FRANKLIN" \
        "$FRANKLIN_TMP"

    rm -f "$FRANKLIN_TMP"

    # Index final Franklin VCF
    bcftools index -f -t "$FRANKLIN"

    # --------------------------------------
    # STEP 19: VERIFY FRANKLIN VCF
    # --------------------------------------

    echo "Franklin VCF contigs for $SAMPLE:"
    bcftools query -f '%CHROM\n' "$FRANKLIN" |
        sort -u | head -30

    echo "Finished sample: $SAMPLE"

done

# ==========================================
# 6. MULTIQC REPORT FOR ALL SAMPLES
# ==========================================

multiqc \
    qc/raw \
    qc/trimmed \
    fastp_reports \
    -o multiqc_reports \
    --force

echo "=========================================="
echo "ALL FIVE SAMPLES COMPLETED"
echo "Franklin VCF files are in output/"
echo "=========================================="

