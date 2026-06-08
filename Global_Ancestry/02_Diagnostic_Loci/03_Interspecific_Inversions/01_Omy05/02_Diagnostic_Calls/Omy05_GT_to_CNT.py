#!/usr/bin/env python3

import pandas as pd
import os

# === Directories ===
input_dir = "./"
diagnostic_dir = "./sets_rapture"
output_dir = "counts_output_inversions"
os.makedirs(output_dir, exist_ok=True)

# === Species sets ===
sets = ["RBT", "WCT", "YCT"]

# === Inversion regions (Omy05 positions in bp) ===
inversion_regions = {
    "Omy05_RBT_Inversion_1": (30102000, 60368999),
    "Omy05_RBT_Inversion_2": (60369000, 87502000),
    "Omy05_RBT_Inversion_full": (30102000, 87502000)
}

# === Debug loci ===
DEBUG_LOCI = ["Omy05_33025000", "Omy05_70000000"]

# === Helper to recode genotype ===
def gt_to_diag_score(gt, diag_allele):
    if diag_allele is None:
        return "NA"
    if gt in {"0/0", "0|0"}:
        return 2 if diag_allele == "ref" else 0
    elif gt in {"1/1", "1|1"}:
        return 2 if diag_allele == "alt" else 0
    elif gt in {"0/1", "1/0", "0|1", "1|0"}:
        return 1
    else:
        return "NA"

# === Process each species set ===
for set_name in sets:
    gt_file = os.path.join(input_dir, f"{set_name}_all_common_input.GT.FORMAT")
    diag_file = os.path.join(diagnostic_dir, f"{set_name}.txt")

    if not os.path.exists(gt_file) or not os.path.exists(diag_file):
        print(f"Skipping {set_name}: missing GT or diagnostic allele file")
        continue

    print(f"\nProcessing {set_name}...")

    # Load GT and diagnostic info
    df = pd.read_csv(gt_file, sep="\t")
    diag_df = pd.read_csv(diag_file, sep="\t", header=None, names=["CHROM", "POS", "DIAG"])

    # Create unique locus IDs
    df["Locus"] = df["CHROM"].astype(str) + "_" + df["POS"].astype(str)
    diag_df["Locus"] = diag_df["CHROM"].astype(str) + "_" + diag_df["POS"].astype(str)

    diag_map = dict(zip(diag_df["Locus"], diag_df["DIAG"]))

    # Filter to only Omy05
    df = df[df["CHROM"] == "Omy05"].copy()

    # === Process each inversion region ===
    for inv_name, (start, end) in inversion_regions.items():
        print(f"{inv_name} ({start:,}-{end:,})")

        # Subset to inversion window
        region_df = df[(df["POS"] >= start) & (df["POS"] <= end)].copy()
        if region_df.empty:
            print(f"No loci found in {inv_name}")
            continue

        # Drop unnecessary columns
        region_df = region_df.drop(columns=["CHROM", "POS"])
        region_df = region_df.set_index("Locus")

        # Recode genotypes to 0/1/2
        for sample in region_df.columns:
            region_df[sample] = [
                gt_to_diag_score(gt, diag_map.get(locus, None))
                for locus, gt in zip(region_df.index, region_df[sample])
            ]

        # Debug check
        debug_subset = region_df.loc[region_df.index.intersection(DEBUG_LOCI)]
        if not debug_subset.empty:
            print("Debug loci (first 5 individuals):")
            print(debug_subset.iloc[:, :5])

        # Save counts per site
        region_df.reset_index(inplace=True)
        region_df.to_csv(os.path.join(output_dir, f"{set_name}_{inv_name}_counts.tsv"),
                         sep="\t", index=False)

        # === Calculate per-sample diagnostic proportions ===
        region_df = region_df.set_index("Locus")
        proportions = {}
        for sample in region_df.columns:
            values = region_df[sample]
            counts = values.value_counts(dropna=True)

            c0 = counts.get(0, 0)
            c1 = counts.get(1, 0)
            c2 = counts.get(2, 0)
            total = c0 + c1 + c2

            if total == 0:
                proportions[sample] = "NA"
            else:
                proportions[sample] = ((2 * c2) + c1) / (2 * total)

        # Save proportions
        pd.DataFrame.from_dict(proportions, orient="index",
                               columns=[f"{set_name}_{inv_name}_prop"]) \
          .to_csv(os.path.join(output_dir, f"{set_name}_{inv_name}_proportions.tsv"),
                  sep="\t")

print("\n Done. Inversion count and proportion tables saved to:", output_dir)
