import csv
import os

samples = [
    "SRR18502497",
    "SRR18502498",
    "SRR18502508",
    "SRR18502509",
    "SRR18502511"
]

outdir = "output/jannovar"

ann_columns = [
    "Allele", "Annotation", "Annotation_Impact", "Gene_Name",
    "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType",
    "Rank", "HGVS.c", "HGVS.p", "cDNA.pos_cDNA.length",
    "CDS.pos_CDS.length", "AA.pos_AA.length", "Distance",
    "ERRORS_WARNINGS_INFO"
]

columns = [
    "Sample", "CHROM", "POS", "REF", "ALT", "QUAL",
    "FILTER", "GT", "DP", "AD_REF", "AD_ALT", "VAF_percent"
] + ann_columns

for sample in samples:
    infile = f"{outdir}/{sample}_Jannovar_annotated.vcf"

    if not os.path.isfile(infile):
        raise FileNotFoundError(infile)

    tsvfile = f"{outdir}/{sample}_Jannovar_VAF_separated.tsv"
    txtfile = f"{outdir}/{sample}_Jannovar_VAF_separated.txt"

    with open(infile) as vcf, \
         open(tsvfile, "w", newline="") as tsv, \
         open(txtfile, "w", newline="") as txt:

        tw = csv.writer(tsv, delimiter="\t")
        xw = csv.writer(txt, delimiter="\t")

        tw.writerow(columns)
        xw.writerow(columns)

        for line in vcf:
            if line.startswith("#"):
                continue

            f = line.rstrip().split("\t")
            chrom, pos, _, ref, alts, qual, filt, info = f[:8]

            format_keys = f[8].split(":")
            sample_values = f[9].split(":")
            fmt = dict(zip(format_keys, sample_values))

            ad = fmt.get("AD", ".").split(",")
            gt = fmt.get("GT", ".")
            dp = fmt.get("DP", ".")

            try:
                ad_nums = [int(x) for x in ad]
            except ValueError:
                ad_nums = []

            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info_dict[key] = value

            ann_entries = info_dict.get("ANN", ".").split(",")

            for alt_index, alt in enumerate(alts.split(","), start=1):
                ad_ref = ad[0] if ad else "."
                ad_alt = ad[alt_index] if len(ad) > alt_index else "."

                vaf = "."
                if len(ad_nums) > alt_index and sum(ad_nums) > 0:
                    vaf = f"{ad_nums[alt_index] / sum(ad_nums) * 100:.2f}"

                annotations = []
                for entry in ann_entries:
                    parts = entry.split("|")
                    parts += [""] * (16 - len(parts))
                    if parts[0] == alt:
                        annotations.append(parts[:16])

                if not annotations:
                    annotations = [[""] * 16]

                for ann in annotations:
                    row = [
                        sample, chrom, pos, ref, alt, qual,
                        filt, gt, dp, ad_ref, ad_alt, vaf
                    ] + ann

                    tw.writerow(row)
                    xw.writerow(row)

    print(f"Completed: {sample}")
    print(f"  {tsvfile}")
    print(f"  {txtfile}")

print("All five samples completed.")
