#!/usr/bin/env python3
import sys
import pandas as pd

in_maf = sys.argv[1]
out_maf = sys.argv[2]

# Read MAF, skip metadata lines
df = pd.read_csv(
    in_maf,
    sep="\t",
    comment="#",
    dtype=str
)

# Convert numeric fields
for col in ["t_depth", "t_alt_count"]:
    df[col] = pd.to_numeric(df.get(col), errors="coerce")

# Detect gnomAD column
af_col = None
for candidate in ["gnomADe_AF", "gnomAD_AF"]:
    if candidate in df.columns:
        af_col = candidate
        df[af_col] = pd.to_numeric(df[af_col], errors="coerce")
        break

# Functional consequence filter
functional_keep = df["Variant_Classification"].isin([
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Ins",
    "Frame_Shift_Del",
    "In_Frame_Ins",
    "In_Frame_Del",
    "Splice_Site",
    "Start_Codon_Del",
])

# Depth & allelic support
coverage_keep = (
    (df["t_depth"] >= 10) & (df["t_alt_count"] >= 5)
)

# Rare variant filter
if af_col:
    freq_keep = df[af_col].isna() | (df[af_col] < 0.01)
else:
    freq_keep = True  # If no AF column, keep all

filtered = df[functional_keep & coverage_keep & freq_keep]

filtered.to_csv(out_maf, sep="\t", index=False)

