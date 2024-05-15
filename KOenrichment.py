### Repository contains workflow for analysing KEGG Orthology (KO) abundance analysis ###
### KO abundances were inferred from PICRUSt2 ###
### PICRUSt2 input was obtained from Qiime2 ###

import pandas as pd
import numpy as np
import requests
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

files = ['wtsd_filtered.txt', 'wtkd_filtered.txt', 'kosd_filtered.txt', 'kokd_filtered.txt']
conditions = ['WTSD', 'WTKD', 'KOSD', 'KOKD']
names = ['wtsd', 'wtkd', 'kosd', 'kokd']

# Read the data and create DataFrames with dynamically generated variable names
for name, file, condition in zip(names, files, conditions):
    globals()[name] = pd.read_csv(file, sep=' ', header=1)
    globals()[name] = globals()[name]

wtsd = wtsd.rename(columns = {'sample-10': 'WTSD_3', 'sample-11': 'WTSD_4', 'sample-12': 'WTSD_5', 'sample-8': 'WTSD_1', 'sample-9':'WTSD_2'})
wtkd = wtkd.rename(columns = {'sample-18': 'WTKD_1', 'sample-19': 'WTKD_2', 'sample-20': 'WTKD_3', 'sample-21': 'WTKD_4', 'sample-22':'WTKD_5'})
kosd = kosd.rename(columns = {'sample-3': 'KOSD_1', 'sample-4': 'KOSD_2', 'sample-5': 'KOSD_3', 'sample-6': 'KOSD_4', 'sample-7':'KOSD_5'})
kokd = kokd.rename(columns = {'sample-13': 'KOKD_1', 'sample-14': 'KOKD_2', 'sample-15': 'KOKD_3', 'sample-16': 'KOKD_4', 'sample-17':'KOKD_5'})

wtsd.index = wtsd['#OTUID']
wtkd.index = wtkd['#OTUID']
kosd.index = kosd['#OTUID']
kokd.index = kokd['#OTUID']

wtsd = wtsd.drop(columns='#OTUID')
wtkd = wtkd.drop(columns='#OTUID')
kosd = kosd.drop(columns='#OTUID')
kokd = kokd.drop(columns='#OTUID')

conditionwtsd = wtsd.sum(axis=1) < 100
conditionwtkd = wtkd.sum(axis=1) < 100
conditionkosd = kosd.sum(axis=1) < 100
conditionkokd = kokd.sum(axis=1) < 100


# Use the drop method to drop columns based on the condition
wtsd_filtered100 = wtsd.drop(index=wtsd.index[conditionwtsd == True])
wtkd_filtered100 = wtkd.drop(index=wtkd.index[conditionwtkd == True])
kosd_filtered100 = kosd.drop(index=kosd.index[conditionkosd == True])
kokd_filtered100 = kokd.drop(index=kokd.index[conditionkokd == True])
wtsd_filtered100

wtsd_filtered100_rel = wtsd_filtered100.div(wtsd_filtered100.sum(axis=0), axis=1)
wtkd_filtered100_rel = wtkd_filtered100.div(wtkd_filtered100.sum(axis=0), axis=1)
kosd_filtered100_rel = kosd_filtered100.div(kosd_filtered100.sum(axis=0), axis=1)
kokd_filtered100_rel = kokd_filtered100.div(kokd_filtered100.sum(axis=0), axis=1)
wtsd_filtered100_rel

wtkd_filtered100_rel['#OTUID'] = wtkd_filtered100_rel.index
wtkd_filtered100_rel['Condition'] = 'WTKD'
wtsd_filtered100_rel['#OTUID'] = wtsd_filtered100_rel.index
wtsd_filtered100_rel['Condition'] = 'WTSD'
kosd_filtered100_rel['#OTUID'] = kosd_filtered100_rel.index
kosd_filtered100_rel['Condition'] = 'KOSD'
kokd_filtered100_rel['#OTUID'] = kokd_filtered100_rel.index
kokd_filtered100_rel['Condition'] = 'KOKD'

wtkd_long = pd.melt(wtkd_filtered100_rel, id_vars=["#OTUID", "Condition"])
wtsd_long = pd.melt(wtsd_filtered100_rel, id_vars=["#OTUID", "Condition"])
kosd_long = pd.melt(kosd_filtered100_rel, id_vars=["#OTUID", "Condition"])
kokd_long = pd.melt(kokd_filtered100_rel, id_vars=["#OTUID", "Condition"])

df_combined = pd.concat([wtsd_long, wtkd_long, kosd_long, kokd_long])

df_combined_wovariable = df_combined.drop(columns=['variable'])

# Perform ANOVA to see if there are any overall significant differences
model = ols('value ~ C(Condition)', data=df_combined_wovariable).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
anova_table

# Perform Tukey's HSD test
tukey_results = pairwise_tukeyhsd(endog=df_combined_wovariable['value'], groups=df_combined_wovariable['Condition'], alpha=0.05)

# Convert the results to a DataFrame for better readability
tukey_df = pd.DataFrame(data=tukey_results._results_table.data[1:], columns=tukey_results._results_table.data[0])
significant_comparisons = tukey_df[tukey_df['reject']]
tukey_df  # Display the first 10 rows for an overview


ko_data = pd.read_csv("ko.biom.tsv", sep='\t')
ko_data

WTSD_samples = ['sample-10', 'sample-11', 'sample-12', 'sample-8', 'sample-9']
WTKD_samples = ['sample-18', 'sample-19', 'sample-20', 'sample-21', 'sample-22']
KOSD_samples = ['sample-3', 'sample-4', 'sample-5', 'sample-6', 'sample-7']
KOKD_samples = ['sample-13', 'sample-14', 'sample-15', 'sample-16', 'sample-17']

# Extract the relevant columns for each sample group
WTSD_data = ko_data.set_index('#OTU ID')[WTSD_samples]
WTKD_data = ko_data.set_index('#OTU ID')[WTKD_samples]
KOSD_data = ko_data.set_index('#OTU ID')[KOSD_samples]
KOKD_data = ko_data.set_index('#OTU ID')[KOKD_samples]

# Display a snippet of one of the groups to verify
WTSD_data

WTSD_data[WTSD_data[WTSD_data == 0].count(axis=1) >= 3] = np.nan
WTKD_data[WTKD_data[WTKD_data == 0].count(axis=1) >= 3] = np.nan
KOSD_data[KOSD_data[KOSD_data == 0].count(axis=1) >= 3] = np.nan
KOKD_data[KOKD_data[KOKD_data == 0].count(axis=1) >= 3] = np.nan

wtsd_filtered_rel = WTSD_data.div(WTSD_data.sum(axis=0), axis=1)
wtkd_filtered_rel = WTKD_data.div(WTKD_data.sum(axis=0), axis=1)
kosd_filtered_rel = KOSD_data.div(KOSD_data.sum(axis=0), axis=1)
kokd_filtered_rel = KOKD_data.div(KOKD_data.sum(axis=0), axis=1)

# Function to perform Mann-Whitney U test and interpret the results
def perform_mannwhitneyu_test(group1, group2):
    """
    Perform Mann-Whitney U test between two groups and return the p-value and direction of change.
    """
    stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
    direction = "increased" if group1.mean() > group2.mean() else "decreased"
    return p_value, direction

# Initialize a dictionary to store the results
significant_changes = {}

# Set the significance level
alpha = 0.05

# Iterate through each KO number and perform the tests
for ko_number in ko_data["#OTU ID"]:
    wt_sd_values = wtsd_filtered_rel.loc[ko_number]
    wt_kd_values = wtkd_filtered_rel.loc[ko_number]
    ko_sd_values = kosd_filtered_rel.loc[ko_number]
    ko_kd_values = kokd_filtered_rel.loc[ko_number]

    # WT vs KO (SD)
    p_value_wt_ko_sd, direction_wt_ko_sd = perform_mannwhitneyu_test(wt_sd_values, ko_sd_values)

    # WT vs KO (KD)
    p_value_wt_ko_kd, direction_wt_ko_kd = perform_mannwhitneyu_test(wt_kd_values, ko_kd_values)

    # SD vs KD (WT)
    p_value_wt_sd_kd, direction_wt_sd_kd = perform_mannwhitneyu_test(wt_sd_values, wt_kd_values)

    # SD vs KD (KO)
    p_value_ko_sd_kd, direction_ko_sd_kd = perform_mannwhitneyu_test(ko_sd_values, ko_kd_values)

    # Check if the changes are significant and store the results
    if p_value_wt_ko_sd < alpha:
        significant_changes[(ko_number, "WT vs KO (SD)")] = direction_wt_ko_sd
    
    if p_value_wt_ko_kd < alpha:
        significant_changes[(ko_number, "WT vs KO (KD)")] = direction_wt_ko_kd

    if p_value_wt_ko_sd < alpha:
        significant_changes[(ko_number, "SD vs KD (WT)")] = direction_wt_sd_kd
    
    if p_value_wt_ko_kd < alpha:
        significant_changes[(ko_number, "SD vs KD (KO)")] = direction_ko_sd_kd

# Display the significant changes
significant_changes

changes_wt_ko_sd = {key: value for key, value in significant_changes.items() if 'WT vs KO (SD)' in key}
changes_wt_ko_kd = {key: value for key, value in significant_changes.items() if 'WT vs KO (KD)' in key}
changes_wt_sd_kd = {key: value for key, value in significant_changes.items() if 'SD vs KD (WT)' in key}
changes_ko_sd_kd = {key: value for key, value in significant_changes.items() if 'SD vs KD (KO)' in key}
changes_wt_ko_sd

# Retrieve KO numbers from significant dictionaries.
def retrieve_ko_from_dicts(dict):
    ko_list = []
    for i in dict.keys():
        ko, cond = i
        ko_list.append(ko)
    return ko_list
ko_changes_wt_ko_sd = retrieve_ko_from_dicts(changes_wt_ko_sd)
ko_changes_wt_ko_sd

# Extract the relevant columns for each sample group and calculate their means
# ko_data.set_index('#OTU ID', inplace=True)
WTSD_mean = wtsd_filtered_rel.mean(axis=1)
WTKD_mean = wtkd_filtered_rel.mean(axis=1)
KOSD_mean = kosd_filtered_rel.mean(axis=1)
KOKD_mean = kokd_filtered_rel.mean(axis=1)

# Calculate log fold changes
# Adding a small constant to avoid division by zero or log of zero
epsilon = 1e-6
logFC_WT_KO_SD = np.log2((WTSD_mean + epsilon) / (KOSD_mean + epsilon))
logFC_WT_KO_KD = np.log2((WTKD_mean + epsilon) / (KOKD_mean + epsilon))
logFC_WT_SD_KD = np.log2((WTSD_mean + epsilon) / (WTKD_mean + epsilon))
logFC_KO_SD_KD = np.log2((KOSD_mean + epsilon) / (KOKD_mean + epsilon))

# Combine results into a DataFrame for easier visualization
logFC_data = pd.DataFrame({
    'KO Number': WTSD_data.index,
    'LogFC WT vs KO (SD)': logFC_WT_KO_SD,
    'LogFC WT vs KO (KD)': logFC_WT_KO_KD,
    'LogFC SD vs KD (WT)': logFC_WT_SD_KD,
    'LogFC SD vs KD (KO)': logFC_KO_SD_KD
}).reset_index(drop=True)

logFC_data

# Creating the heatmap for all KO numbers
plt.figure(figsize=(10, 15))
heatmap_data = logFC_data.loc[logFC_data['KO Number'].isin(["K00022", "K00626", "K01641", "K01640", "K00019", "K01027", "K01034", "0K1035",
                                             "K01843", "K01844", "K18012", "K18013", "K19709", "K23756", "K01896", "K01913", 
                                             "K00248", "K00209", "K17829", "K01907"])].set_index('KO Number').dropna()
# Creating the heatmap
sns.heatmap(heatmap_data, cmap="viridis", annot=True)
plt.title('Log2 Fold Changes of KO Enzymes')
plt.xlabel('Comparisons')
plt.ylabel('KO Numbers')
plt.show()

logFC_data.loc[logFC_data['KO Number'].isin(["K00022", "K00626", "K01641", "K01640", "K00019", "K01027", "K01034", "0K1035",
                                             "K01843", "K01844", "K18012", "K18013", "K19709", "K23756", "K01896", "K01913", 
                                             "K00248", "K00209", "K17829", "K01907", "K01782", "K01825"])]

