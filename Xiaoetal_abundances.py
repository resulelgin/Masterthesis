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

stat1, p_value1 = shapiro(subsetchow[subsetchow["K_ID"] == "K01641"]["Value"])
print(stat1, p_value1)
stat2, p_value2 = shapiro(subsethf[subsethf["K_ID"] == "K01641"]["Value"])
print(stat2, p_value2)
stat3, p_value3 = shapiro(subsetchow[subsetchow["K_ID"] == "K01640"]["Value"])
print(stat3, p_value3)
stat4, p_value4 = shapiro(subsethf[subsethf["K_ID"] == "K01640"]["Value"])
print(stat4, p_value4)
stat5, p_value5 = shapiro(subsetchow[subsetchow["K_ID"] == "K00019"]["Value"])
print(stat5, p_value5)
stat6, p_value6 = shapiro(subsethf[subsethf["K_ID"] == "K00019"]["Value"])
print(stat6, p_value6)
stat7, p_value7 = shapiro(subsetchow[subsetchow["K_ID"] == "K01034"]["Value"])
print(stat7, p_value7)
stat8, p_value8 = shapiro(subsethf[subsethf["K_ID"] == "K01034"]["Value"])
print(stat8, p_value8)
stat9, p_value9 = shapiro(subsetchow[subsetchow["K_ID"] == "K00626"]["Value"])
print(stat9, p_value9)
stat10, p_value10 = shapiro(subsethf[subsethf["K_ID"] == "K00626"]["Value"])
print(stat10, p_value10)

U1, p1 = mannwhitneyu(subsetchow[subsetchow["K_ID"] == "K01641"]["Value"], subsethf[subsethf["K_ID"] == "K01641"]["Value"], method="exact")
print(U1, p1)
U2, p2 = mannwhitneyu(subsetchow[subsetchow["K_ID"] == "K01640"]["Value"], subsethf[subsethf["K_ID"] == "K01640"]["Value"], method="exact")
print(U2, p2)
U3, p3 = mannwhitneyu(subsetchow[subsetchow["K_ID"] == "K00019"]["Value"], subsethf[subsethf["K_ID"] == "K00019"]["Value"], method="exact")
print(U3, p3)
U4, p4 = mannwhitneyu(subsetchow[subsetchow["K_ID"] == "K01034"]["Value"], subsethf[subsethf["K_ID"] == "K01034"]["Value"], method="exact")
print(U4, p4)
U5, p5 = mannwhitneyu(subsetchow[subsetchow["K_ID"] == "K00626"]["Value"], subsethf[subsethf["K_ID"] == "K00626"]["Value"], method="exact")
print(U5, p5)


