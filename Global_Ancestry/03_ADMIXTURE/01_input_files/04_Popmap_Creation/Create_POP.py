# make_pop_file.py
import pandas as pd

# File paths
sample_file = "sample_ids.txt"   # one sample ID per line, in same order as .bed
pop_file = "pops.txt"        # tab-delimited, with headers: Sample_ID   Pop
output_file = "output.pop"

# Read files
samples = pd.read_csv(sample_file, header=None, names=["Sample_ID"])
pop_map = pd.read_csv(pop_file, sep="\t")

# Merge to assign populations
merged = samples.merge(pop_map, on="Sample_ID", how="left")

# Fill missing populations with "-"
merged["Pop"].fillna("-", inplace=True)

# Output the .pop file (one pop per line)
merged["Pop"].to_csv(output_file, index=False, header=False)
print(f".pop file written to {output_file}")

