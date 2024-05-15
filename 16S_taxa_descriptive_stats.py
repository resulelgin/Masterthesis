### Repository contains workflow for 16S rRNA sequencing taxonomical abundance descriptive statistics ###
### Input is obtained from Qiime2View ###

import pandas as pd
import openpyxl

data_df = pd.read_csv("level-3.csv") 
metadata_df = pd.read_csv("metadata.txt", sep="\t")

# Merge the datasets on the sample identifier
merged_df = pd.merge(data_df, metadata_df, how="left", left_on="index", right_on="Samples")

# Calculate relative abundances
numeric_cols = merged_df.select_dtypes(include="number").columns.tolist()
merged_df[numeric_cols] = merged_df[numeric_cols].div(merged_df[numeric_cols].sum(axis=1), axis=0)

# Aggregate data by Group for descriptives
grouped = merged_df.groupby("Group")[numeric_cols].agg(["mean", "std", "min", "max"])

# The `grouped` DataFrame now contains the mean, standard deviation, minimum, and maximum relative abundance for each taxon, grouped by WTSD, KOSD, WTKD, KOKD.
with pd.ExcelWriter("descriptive_stats_level3.xlsx", engine="openpyxl") as writer:
    grouped.to_excel(writer, sheet_name="Descriptive Stats")


