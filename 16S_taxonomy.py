### Workflow for an in-depth analysing 16S rRNA sequencing taxonomical abundance analysis across study groups ###
### Input is obtained from Qiime2View ###

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from skbio.stats.composition import ancom
from statsmodels.formula.api import ols
from skbio.stats.composition import multiplicative_replacement
import scipy
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from biom import Table
from biom import load_table
import h5py
import mygene

# Load the dataset
file_path = 'level7_taxa.csv'
data = pd.read_csv(file_path)

# Display the first few rows of the dataset to understand its structure
data

level6 = pd.read_csv('level6_taxa.csv')
#data = data[data["index"] != "sample-1"]
#data = data[data["index"] != "sample-2"]
#data = data[data["index"] != "sample-45"]
data.set_index("index", inplace=True)
data
new_sample_labels = []
genotype_diet_counts = {}

for index in data.index:
    if index == 'sample-1':
        new_label = 'SDC'
    elif index == 'sample-2':
        new_label = 'KDC'
    elif index == 'sample-45':
        new_label = 'Zymo'
    else:
        # Accessing other column values using the current index
        genotype = data.at[index, 'Genotype']
        diet = data.at[index, 'Diet']
        genotype_diet_key = f"{genotype}_{diet}"

        genotype_diet_counts[genotype_diet_key] = genotype_diet_counts.get(genotype_diet_key, 0) + 1
        new_label = f"{genotype_diet_key}_{genotype_diet_counts[genotype_diet_key]}"

    new_sample_labels.append(new_label)

# Update the DataFrame index
data.index = new_sample_labels
level6['d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia']

data_dropped = data.drop(["forward-absolute-filepath", "reverse-absolute-filepath", "Diet", "Genotype"], axis=1)
# species_16S = []
# for column in data.columns:
#    species_16S.append(column)
# with open("species_16S.txt", 'w') as file:
#        for item in species_16S:
#            file.write(f"{item}\n")

# Selecting only the abundance columns for the heatmap
abundance_data = data.drop(['forward-absolute-filepath', 'reverse-absolute-filepath', 'Diet', 'Genotype'], axis=1)
abundance_data
abundance_data['d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Eubacterium_sp.']
abundance_data = abundance_data.drop(['sample-1', 'sample-2', 'sample-45'])
abundance_data
taxa_abundances = abundance_data.sum(axis=0).sort_values(ascending=False)
taxa_abundances
condition = abundance_data.sum(axis=0) < 100

# Use the drop method to drop columns based on the condition
abundance_data_filtered_zeros = abundance_data.drop(columns=abundance_data.columns[condition == True])
abundance_data_filtered_zeros
relative_abundance_data_filtered_zeros = abundance_data_filtered_zeros.div(abundance_data_filtered_zeros.sum(axis=1), axis=0)
relative_abundance_data_filtered_zeros
relative_abundance_data_filtered_zeros['d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Lachnospiraceae_FCS020_group;s__mouse_gut']
wtko_sd_filtered = abundance_data_filtered_zeros.loc[['WT_SD_1', 'WT_SD_2', 'WT_SD_3', 'WT_SD_4', 'WT_SD_5',
                                   'KO_SD_1', 'KO_SD_2', 'KO_SD_3', 'KO_SD_4', 'KO_SD_5']]
wtko_sd_fitered_perc = wtko_sd_filtered.div(wtko_sd_filtered.sum(axis=1), axis=0)
wtko_sd_fitered_perc['d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridia_vadinBB60_group;f__Clostridia_vadinBB60_group;g__Clostridia_vadinBB60_group;s__uncultured_Clostridia']
wtko_sd_filtered.to_csv('wtko_sd_filtered_level7.csv')
relative_abundance_data_filtered_zeros = abundance_data_filtered_zeros.div(abundance_data_filtered_zeros.sum(axis=1), axis=0)
relative_abundance_data_filtered_zeros_transposed = relative_abundance_data_filtered_zeros.transpose()
# relative_abundance_data_filtered_zeros_transposed.to_csv('16s_taxa_relative_abundance_filtered100_transposed.csv')

wtko_sd_clusterdf = relative_abundance_data_filtered_zeros.loc[['WT_SD_1', 'WT_SD_2', 'WT_SD_3', 'WT_SD_4', 'WT_SD_5', 'KO_SD_1', 'KO_SD_2',
                                                  'KO_SD_3', 'KO_SD_4', 'KO_SD_5']]
wtko_sd_clusterdf['Ground'] = ['WT', 'WT', 'WT', 'WT', 'WT', 'KO', 'KO',
                                                  'KO', 'KO', 'KO']
X = wtko_sd_clusterdf.iloc[:, :-1]
y_true = wtko_sd_clusterdf.iloc[:, -1]
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
k = 2  # for example, change this based on your dataset

# Perform K-Means clustering
kmeans = KMeans(n_clusters=k, random_state=42)
y_pred = kmeans.fit_predict(X_scaled)

# Compare the clustering labels with the real labels
ari_score = adjusted_rand_score(y_true, y_pred)

print(f"Adjusted Rand Index: {ari_score:.4f}")
pca = PCA(n_components=2)  # Reduce to 2 dimensions for plotting
X_pca = pca.fit_transform(X_scaled)

# Plotting the clusters
fig, ax = plt.subplots(1, 2, figsize=(16, 6))

# Plot predicted clusters
for i in range(k):  # Assuming 'k' is the number of clusters used in KMeans
    ax[0].scatter(X_pca[y_pred == i, 0], X_pca[y_pred == i, 1], label=f'Cluster {i}')
ax[0].set_title('Predicted Clusters')
ax[0].set_xlabel('PCA 1')
ax[0].set_ylabel('PCA 2')
ax[0].legend()
if not pd.api.types.is_numeric_dtype(y_true):
    y_true_numeric = pd.Categorical(y_true).codes
else:
    y_true_numeric = y_true

unique_labels = set(y_true_numeric)
for i in unique_labels:
    ax[1].scatter(X_pca[y_true_numeric == i, 0], X_pca[y_true_numeric == i, 1], label=f'Label {i}')
ax[1].set_title('Actual Clusters')
ax[1].set_xlabel('PCA 1')
ax[1].set_ylabel('PCA 2')
ax[1].legend()

plt.show()

def biplot(score, coeff, labels=None):
    """
    Plot a biplot showing the score and the coefficients of the principal components.
    
    :param score: PCA scores, the transformed coordinates of the samples in the principal component space.
    :param coeff: Loadings of the original variables on the principal components.
    :param labels: Labels of the original features.
    """
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot scores
    ax.scatter(xs, ys, s=5)
    
    # Plot arrows and labels for coefficients
    for i in range(n):
        ax.arrow(0, 0, coeff[i,0]*max(xs), coeff[i,1]*max(ys), color='r', alpha=0.5)
        if labels is None:
            ax.text(coeff[i,0]*max(xs)*1.15, coeff[i,1]*max(ys)*1.15, "Var"+str(i+1), color='g', ha='center', va='center')
        else:
            ax.text(coeff[i,0]*max(xs)*1.15, coeff[i,1]*max(ys)*1.15, labels[i], color='g', ha='center', va='center')
    
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA Biplot")
    ax.grid(True)

# Assuming X_scaled is your scaled dataset
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Call the biplot function
biplot(X_pca, pca.components_.T, labels=wtko_sd_clusterdf.columns[:-1])  # Assuming df.columns[:-1] are your feature names
plt.show()

# Perform PCA without reducing the number of components to retain all variance information
pca_full = PCA()
X_pca_full = pca_full.fit_transform(X_scaled)

# Calculate the percentage of variance explained by each of the selected components
variance_explained = pca_full.explained_variance_ratio_

# Calculate cumulative variance explained
cumulative_variance = np.cumsum(variance_explained)

# Plotting
fig, ax1 = plt.subplots(figsize=(10, 7))

# Individual explained variance
ax1.bar(range(1, len(variance_explained) + 1), variance_explained, alpha=0.6, label='Individual Explained Variance')
ax1.set_xlabel('Principal Component')
ax1.set_ylabel('Explained Variance Ratio')
ax1.set_title('PCA Explained Variance')

# Make a second y-axis for the cumulative variance
ax2 = ax1.twinx()
ax2.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, 'r-o', label='Cumulative Explained Variance')
ax2.set_ylabel('Cumulative Explained Variance Ratio')

# Adding legends
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

plt.show()

def biplot(score, coeff, labels=None, arrow_scale=1, text_size=8):
    """
    Generate a biplot.
    
    :param score: PCA scores.
    :param coeff: Loadings of the original variables.
    :param labels: Labels of the original features.
    :param arrow_scale: Factor to scale the arrows (default is 1).
    :param text_size: Text size for labels.
    """
    xs = score[:, 0]
    ys = score[:, 1]
    n = coeff.shape[0]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.scatter(xs, ys, s=5)  # Plot scores as scatter
    
    # Plot arrows and labels
    for i in range(n):
        ax.arrow(0, 0, coeff[i, 0] * max(xs) * arrow_scale, coeff[i, 1] * max(ys) * arrow_scale, 
                 color='r', alpha=0.5)
        if labels is None:
            ax.text(coeff[i, 0] * max(xs) * arrow_scale * 1.15, 
                    coeff[i, 1] * max(ys) * arrow_scale * 1.15, 
                    "Var" + str(i + 1), color='g', ha='center', va='center', size=text_size)
        else:
            ax.text(coeff[i, 0] * max(xs) * arrow_scale * 1.15, 
                    coeff[i, 1] * max(ys) * arrow_scale * 1.15, 
                    labels[i], color='g', ha='center', va='center', size=text_size)
    
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA Biplot")
    ax.grid(True)

# Assuming X_scaled and pca.components_ are defined as before
biplot(X_pca, pca.components_.T, labels=wtko_sd_clusterdf.columns[:-1], arrow_scale=1, text_size=8)
plt.show()

# Calculate relative abundances
relative_abundances = abundance_data.div(abundance_data.sum(axis=1), axis=0)
perc_relative_abundances = abundance_data.div(abundance_data.sum(axis=1), axis=0) * 100

# Identify top 20 taxonomic groups based on relative abundances
top_taxa = perc_relative_abundances.sum(axis=0).sort_values(ascending=False).head(20).index

# Selecting only the top 20 taxonomic groups for the heatmap
top_taxa_data = perc_relative_abundances[top_taxa]
top_taxa_data
# relative_abundances[relative_abundances[relative_abundances == 0] == 0] = np.nan
perc_relative_abundances
relative_abundances_transposed = relative_abundances.transpose()
relative_abundances_transposed
# relative_abundances_transposed.to_csv('16S_relative_abundance_transposed.csv')
plt.figure(figsize=(20, 10))
sns.heatmap(top_taxa_data.T, cmap="YlGnBu", annot=True, fmt=".2f", xticklabels=new_sample_labels)
plt.title('Top 20 Taxonomic Groups by Relative Abundance')
plt.xlabel('Samples')
plt.xticks(rotation=45)  # Rotating the x-axis labels for better readability
plt.show()
ancom_data = relative_abundances
data_WTKO_SD = data_dropped.loc[["WT_SD_1", "WT_SD_2", "WT_SD_3", "WT_SD_4", "WT_SD_5",
                                   "KO_SD_1", "KO_SD_2", "KO_SD_3", "KO_SD_4", "KO_SD_5"]]
grouping_WTKO_SD = grouping = pd.Series(['wt_sd', 'wt_sd', 'wt_sd', 'wt_sd', 'wt_sd',
                                           'ko_sd', 'ko_sd', 'ko_sd', 'ko_sd', 'ko_sd'],
                     index=["WT_SD_1", "WT_SD_2", "WT_SD_3", "WT_SD_4", "WT_SD_5",
                            "KO_SD_1", "KO_SD_2", "KO_SD_3", "KO_SD_4", "KO_SD_5"])
data_WTKO_SD
# Epsilon addition approach.
epsilon = 0.00001
data_WTKO_SD_e = data_WTKO_SD + epsilon
data_WTKO_SD_e
# constant_WTKO_SD = data_WTKO_SD_e.columns[data_WTKO_SD_e.nunique() == 1]
# data_WTKO_SD_e.drop(columns=constant_WTKO_SD, inplace=True)
# data_WTKO_SD_e
ancom_WTKO_SD, percentile_WTKO_SD = ancom(data_WTKO_SD_e, grouping_WTKO_SD)
ancom_WTKO_SD.sort_values(by='W', ascending=False)
ancom_WTKO_SD[ancom_WTKO_SD["Reject null hypothesis"] == True]
percentile_WTKO_SD[50.0].loc['d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Alistipes;s__Alistipes_finegoldii']
data_WT_SDKD =  data_dropped.loc[["WT_SD_1", "WT_SD_2", "WT_SD_3", "WT_SD_4", "WT_SD_5",
                                   "WT_KD_1", "WT_KD_2", "WT_KD_3", "WT_KD_4", "WT_KD_5"]]
grouping_WT_SDKD = pd.Series(['wt_sd', 'wt_sd', 'wt_sd', 'wt_sd', 'wt_sd',
                                           'wt_kd', 'wt_kd', 'wt_kd', 'wt_kd', 'wt_kd'],
                     index=["WT_SD_1", "WT_SD_2", "WT_SD_3", "WT_SD_4", "WT_SD_5",
                            "WT_KD_1", "WT_KD_2", "WT_KD_3", "WT_KD_4", "WT_KD_5"])
data_WT_SDKD
data_WT_SDKD_e = data_WT_SDKD + epsilon
data_WT_SDKD_e
ancom_WT_SDKD, percentile_WT_SDKD = ancom(data_WT_SDKD_e, grouping_WT_SDKD)
ancom_WT_SDKD[ancom_WT_SDKD["Reject null hypothesis"] == True]
data_KO_SDKD =  data_dropped.loc[['KO_SD_1', 'KO_SD_2', 'KO_SD_3', 'KO_SD_4', 'KO_SD_5',
                                   'KO_KD_1', 'KO_KD_2', 'KO_KD_3', 'KO_KD_4', 'KO_KD_5']]
grouping_KO_SDKD = pd.Series(['ko_sd', 'ko_sd', 'ko_sd', 'ko_sd', 'ko_sd',
                              'ko_kd', 'ko_kd', 'ko_kd', 'ko_kd', 'ko_kd'],
                     index=['KO_SD_1', 'KO_SD_2', 'KO_SD_3', 'KO_SD_4', 'KO_SD_5',
                            'KO_KD_1', 'KO_KD_2', 'KO_KD_3', 'KO_KD_4', 'KO_KD_5'])
data_KO_SDKD
data_KO_SDKD_e = data_KO_SDKD + epsilon
data_KO_SDKD_e
ancom_KO_SDKD, percentile_KO_SDKD = ancom(data_KO_SDKD_e, grouping_KO_SDKD)
ancom_KO_SDKD[ancom_KO_SDKD['Reject null hypothesis'] == True]
percentile_KO_SDKD[50.0].loc['d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Sellimonas;s__uncultured_bacterium']
data_WTKO_KD =  data_dropped.loc[['WT_KD_1', 'WT_KD_2', 'WT_KD_3', 'WT_KD_4', 'WT_KD_5',
                                   'KO_KD_1', 'KO_KD_2', 'KO_KD_3', 'KO_KD_4', 'KO_KD_5']]
grouping_WTKO_KD = pd.Series(['wt_kd', 'wt_kd', 'wt_kd', 'wt_kd', 'wt_kd',
                              'ko_kd', 'ko_kd', 'ko_kd', 'ko_kd', 'ko_kd'],
                     index=['WT_KD_1', 'WT_KD_2', 'WT_KD_3', 'WT_KD_4', 'WT_KD_5',
                            'KO_KD_1', 'KO_KD_2', 'KO_KD_3', 'KO_KD_4', 'KO_KD_5'])
data_WTKO_KD
data_WTKO_KD_e = data_WTKO_KD + epsilon
data_WTKO_KD_e
ancom_WTKO_KD, percentile_WTKO_KD = ancom(data_WTKO_KD_e, grouping_WTKO_KD)
ancom_WTKO_KD[ancom_WTKO_KD['Reject null hypothesis'] == True]
percentile_WTKO_KD[50.0].loc['d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Angelakisella;s__Ruminococcus_sp.']
# grouping = pd.Series(['chow_control', 'wt_sd', 'wt_sd', 'wt_sd', 'ko_kd', 'ko_kd',
#                      'ko_kd', 'ko_kd', 'ko_kd', 'wt_kd', 'wt_kd', 'kd_control',
#                      'wt_kd', 'wt_kd', 'wt_kd', 'ko_sd', 'ko_sd', 'zymo', 'ko_sd',
#                      'ko_sd', 'ko_sd', 'wt_sd', 'wt_sd'],
#                     index=['SDC_SDC_1', 'WT_SD_1', 'WT_SD_2', 'WT_SD_3', 'KO_KD_1', 'KO_KD_2',
#       'KO_KD_3', 'KO_KD_4', 'KO_KD_5', 'WT_KD_1', 'WT_KD_2', 'KDC_KDC_1',
#       'WT_KD_3', 'WT_KD_4', 'WT_KD_5', 'KO_SD_1', 'KO_SD_2', 'Zymo_Zymo_1',
#       'KO_SD_3', 'KO_SD_4', 'KO_SD_5', 'WT_SD_4', 'WT_SD_5'])
wtko_sd = abundance_data.loc[['WT_SD_1', 'WT_SD_2', 'WT_SD_3', 'WT_SD_4', 'WT_SD_5', 'KO_SD_1', 'KO_SD_2',
                                                  'KO_SD_3', 'KO_SD_4', 'KO_SD_5']]
wtkosd_ra = wtko_sd.div(wtko_sd.sum(axis=1), axis=0)
wtkosd_ra
wtko_sd['d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Eubacterium_plexicaudatum']
wtkosd_ra['d__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Roseburia;s__Eubacterium_plexicaudatum']
