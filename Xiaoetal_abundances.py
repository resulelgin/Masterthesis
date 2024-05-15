### Repository contains workflow for analysing Xiao et al (DOI: 10.1038/nbt.3353) database KEGG Orthology (KO) abundance ###
### Emphasizing an efficient workflow to filter/visualize the dataframes filtered for the genes of interest ###

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
with open("ChowKO_nonindexed.txt", "r") as file:
    for line in file:
        row = line.strip().split('\t')
        row = [scientific_notation(value) for value in row]
        data.append(row)


# Create a DataFrame
dfchow = pd.DataFrame(data, columns=["K_ID"] + [f"Value{i}" for i in range(1, len(data[0]))])
dfchow_long = pd.melt(dfchow, id_vars=["K_ID"], var_name="Variable", value_name="Value")
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
with open("HFKO_nonindexed.txt", "r") as filehf:
    for line in filehf:
        row = line.strip().split('\t')
        row = [scientific_notation(value) for value in row]
        datahf.append(row)


# Create a DataFrame
dfhf = pd.DataFrame(datahf, columns=["K_ID"] + [f"Value{i}" for i in range(1, len(datahf[0]))])
dfhf_long = pd.melt(dfhf, id_vars=["K_ID"], var_name="Variable", value_name="Value")
# Display the resulting DataFrame
print(dfhf_long)

subsethf = dfhf_long[dfhf_long["K_ID"].isin(["K01034", "K01641", "K01640", "K00019", "K01027", "K00626"])]
subsetchow = dfchow_long[dfchow_long["K_ID"].isin(["K01034", "K01641", "K01640", "K00019", "K01027", "K00626"])]

plt.figure(figsize=(12, 6))
sns.boxplot(x="Value", y="K_ID", data=subsethf, showfliers=True)
#plt.xscale("log")
#plt.xlim(-0.001, 0.001)
plt.xlabel("Relative Abundances")
plt.ylabel("Kegg Orthology Identifier")
plt.title("Relative Abundances of KO Identifiers in High Fat Diet")
plt.show()

plt.figure(figsize=(12, 6))
sns.boxplot(x="Value", y="K_ID", data=subsetchow, showfliers=True)
#plt.xscale("log")
#plt.xlim(-0.001, 0.001)
plt.xlabel("Relative Abundances")
plt.ylabel("Kegg Orthology Identifier")
plt.title("Relative Abundances of KO Identifiers in Standard Diet")
plt.show()

subsetchow["Dataset"] = "Standard Diet"
subsethf["Dataset"] = "High Fat Diet"
combined_df = pd.concat([subsetchow, subsethf])
combined_df

plt.figure(figsize=(10, 6))
ax = sns.boxplot(x='K_ID', y='Value', hue='Dataset', data=combined_df, dodge=True)


plt.title('Swarm Plot of Two Datasets (Long Format)')
plt.legend(title='Dataset')
plt.show()

