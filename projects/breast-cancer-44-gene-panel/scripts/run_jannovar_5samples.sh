#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT="$HOME/gene_panel_5samples"
JAR="$HOME/gene_panel/jannovar/jannovar-cli/target/jannovar-cli-0.42-SNAPSHOT.jar"
DB="$HOME/gene_panel/jannovar/data/refseq_109_hg38.ser"

cd "$PROJECT"

SAMPLES=(
    SRR18502497
    SRR18502498
    SRR18502508
    SRR18502509
    SRR18502511
)

mkdir -p output/jannovar

for SAMPLE in "${SAMPLES[@]}"; do

    INPUT="output/${SAMPLE}_filtered.vcf.gz"
    ANNOTATED="output/jannovar/${SAMPLE}_Jannovar_annotated.vcf"

    echo "===================================="
    echo "Jannovar annotation: $SAMPLE"
    echo "===================================="

    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: Missing input: $INPUT"
        exit 1
    fi

    java -jar "$JAR" annotate-vcf \
        -i "$INPUT" \
        -o "$ANNOTATED" \
        -d "$DB"

    echo "Jannovar completed: $SAMPLE"

done

echo "All five Jannovar annotations completed."

python3 - <<'PY'
import os
import csv

samples = [
    "SRR18502497",
    "SRR18502498",
    "SRR18502508",
    "SRR18502509",
    "SRR18502511"
]

outdir = "output/jannovar"

columns = [
    "Sample", "CHROM", "POS", "REF", "ALT", "QUAL",
    "FILTER", "GT", "DP", "AD_REF", "AD_ALT",
    "VAF_percent", "ANN"
]

for sample in samples:
    infile = f"{outdir}/{sample}_Jannovar_annotated.vcf"
    tsvfile = f"{outdir}/{sample}_Jannovar_VAF.tsv"
    txtfile = f"{outdir}/{sample}_Jannovar_VAF.txt"

    if not os.path.isfile(infile):
        raise FileNotFoundError(infile)

    with open(infile) as vcf, \
         open(tsvfile, "w", newline="") as tsv, \
         open(txtfile, "w", newline="") as txt:

        writer = csv.writer(tsv, delimiter="\t")
        txtwriter = csv.writer(txt, delimiter="\t")

        writer.writerow(columns)
        txtwriter.writerow(columns)

        sample_index = None

        for line in vcf:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
                if len(header) < 10:
                    raise ValueError(f"No sample FORMAT column in {infile}")
                sample_index = 9
                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip().split("\t")

            chrom, pos, _, ref, alts, qual, filt, info = fields[:8]
            format_keys = fields[8].split(":")
            sample_values = fields[sample_index].split(":")

            fmt = dict(zip(format_keys, sample_values))

            ad_values = fmt.get("AD", ".").split(",")
            dp = fmt.get("DP", ".")
            gt = fmt.get("GT", ".")

            ad_ref = ad_values[0] if ad_values else "."
            try:
                ad_numbers = [int(x) for x in ad_values]
            except ValueError:
                ad_numbers = []

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info_dict[key] = value

            ann = info_dict.get("ANN", ".")

            for alt_index, alt in enumerate(alts.split(","), start=1):
                ad_alt = (
                    ad_values[alt_index]
                    if len(ad_values) > alt_index else "."
                )

                vaf = "."
                if ad_numbers and len(ad_numbers) > alt_index:
                    total_ad = sum(ad_numbers)
                    if total_ad > 0:
                        vaf = f"{(ad_numbers[alt_index] / total_ad) * 100:.2f}"

                row = [
                    sample, chrom, pos, ref, alt, qual,
                    filt, gt, dp, ad_ref, ad_alt, vaf, ann
                ]

                writer.writerow(row)
                txtwriter.writerow(row)

    print(f"Created: {tsvfile}")
    print(f"Created: {txtfile}")

print("All five TXT and TSV files created.")
PY

echo "===================================="
echo "WORKFLOW COMPLETED"
echo "===================================="
