### Repository contains eggnog enzyme ID analysis of mouse gut metagenomic sequencing outputs ###
### See Xiao et al. 2015, A catalog of the mouse gut metagenome for datasets ###

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statannot import add_stat_annotation
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu
import scipy.stats as stats

# Define a custom conversion function to handle scientific notation
def scientific_notation(value):
    try:
        return float(value)
    except:
        return value

# Read the data from the text file, applying the custom conversion function
data = []
with open("Choweggnog_nonindexed.txt", "r") as file:
    for line in file:
        row = line.strip().split('\t')
        row = [scientific_notation(value) for value in row]
        data.append(row)


# Create a DataFrame
dfchow = pd.DataFrame(data, columns=["C_ID"] + [f"Value{i}" for i in range(1, len(data[0]))])
dfchow_long = pd.melt(dfchow, id_vars=["C_ID"], var_name="Variable", value_name="Value")
# Display the resulting DataFrame
print(dfchow_long)

# Define a custom conversion function to handle scientific notation
def scientific_notation(value):
    try:
        return float(value)
    except:
        return value

# Read the data from the text file, applying the custom conversion function
datahf = []
with open("Hfeggnog_nonindexed.txt", "r") as filehf:
    for line in filehf:
        row = line.strip().split('\t')
        row = [scientific_notation(value) for value in row]
        datahf.append(row)


# Create a DataFrame
dfhf = pd.DataFrame(datahf, columns=["C_ID"] + [f"Value{i}" for i in range(1, len(datahf[0]))])
dfhf_long = pd.melt(dfhf, id_vars=["C_ID"], var_name="Variable", value_name="Value")
# Display the resulting DataFrame
print(dfhf_long)

subsethf = dfhf_long[dfhf_long["C_ID"].isin(["COG3425", "COG0183", "COG0119", "KOG1610"])]
subsetchow = dfchow_long[dfchow_long["C_ID"].isin(["COG3425", "COG0183", "COG0119", "KOG1610"])]

subsetchow["Dataset"] = "Standard Diet"
subsethf["Dataset"] = "High Fat Diet"
combined_df = pd.concat([subsetchow, subsethf])
combined_df

plt.figure(figsize=(10, 6))
ax = sns.boxplot(x='C_ID', y='Value', hue='Dataset', data=combined_df, dodge=True)


plt.title('Swarm Plot of Two Datasets (Long Format)')
plt.legend(title='Dataset')
plt.show()
