#!/usr/bin/env python3

import pandas as pd
import os

input_dir = "./"
diagnostic_dir = "./sets_rapture"
output_dir = "counts_output"
os.makedirs(output_dir, exist_ok=True)

sets = ["Omy20"]

# Define debug loci for checking
DEBUG_LOCI = ["Omy01_3302505", "Omy02_20694331"]

for set_name in sets:
    gt_file = os.path.join(input_dir, f"{set_name}_all_common_input.GT.FORMAT")
    diag_file = os.path.join(diagnostic_dir, f"{set_name}.txt")
    if not os.path.exists(gt_file) or not os.path.exists(diag_file):
        print(f"Skipping {set_name}: missing GT or diagnostic allele file")
        continue

    print(f"Processing {set_name}...")

    # Load genotype file
    df = pd.read_csv(gt_file, sep="\t")

    # Create unique site identifier
    df["Locus"] = df["CHROM"].astype(str) + "_" + df["POS"].astype(str)

    # Load diagnostic allele info (3-column: CHROM POS ref/alt)
    diag_df = pd.read_csv(diag_file, sep="\t", header=None, names=["CHROM", "POS", "DIAG"])
    diag_df["Locus"] = diag_df["CHROM"].astype(str) + "_" + diag_df["POS"].astype(str)

    # Create mapping from Locus to DIAG allele (ref or alt)
    diag_map = dict(zip(diag_df["Locus"], diag_df["DIAG"]))

    # Drop CHROM, POS and reorder
    df = df.drop(columns=["CHROM", "POS"])
    df = df.set_index("Locus")

    # Function to recode genotype based on diagnostic allele
    def gt_to_diag_score(gt, diag_allele):
        if gt in {"0/0", "0|0"}:
            return 2 if diag_allele == "ref" else 0
        elif gt in {"1/1", "1|1"}:
            return 2 if diag_allele == "alt" else 0
        elif gt in {"0/1", "1/0", "0|1", "1|0"}:
            return 1
        else:
            return "NA"

    # Recode all genotype columns
    for sample in df.columns:
        df[sample] = [
            gt_to_diag_score(gt, diag_map.get(locus, None))
            for locus, gt in zip(df.index, df[sample])
        ]

    # Debug print for diagnostic SNPs
    print(f"\nDebug output for {set_name}:")
    debug_subset = df.loc[df.index.intersection(DEBUG_LOCI)]
    print(debug_subset.iloc[:, :5])  # show only first 5 individuals

    # Save counts
    df.reset_index(inplace=True)
    df.to_csv(os.path.join(output_dir, f"{set_name}_counts.tsv"), sep="\t", index=False)

    # === NEW: Calculate diagnostic proportions ===
    df = df.set_index("Locus")  # reset index for processing
    proportions = {}
    for sample in df.columns:
        values = df[sample]
        counts = values.value_counts(dropna=True)

        count_0 = counts.get(0, 0)
        count_1 = counts.get(1, 0)
        count_2 = counts.get(2, 0)
        total_sites = count_0 + count_1 + count_2

        if total_sites == 0:
            proportions[sample] = "NA"
        else:
            proportions[sample] = ((2 * count_2) + count_1) / (2 * total_sites)

    # Save diagnostic proportion output
    pd.DataFrame.from_dict(proportions, orient="index", columns=[f"{set_name}_diagnostic_proportion"]) \
        .to_csv(os.path.join(output_dir, f"{set_name}_proportions.tsv"), sep="\t")

print("Done.")
