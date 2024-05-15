### Repository contains a workflow for taxonomical barplot visual output for 16S rRNA sequencing ###
### Qiime2View barplot csv file was used as an input ###

import pandas as pd

# Load the CSV data file
data = pd.read_csv(data)

# Load the metadata file
metadata_path = '/mnt/data/metadata.txt'
metadata = pd.read_csv(metadata_path, sep='\t')

# Merge the abundance data with metadata to map the sample names to their groups
merged_data = data.merge(metadata, left_on='index', right_on='Samples')

# Drop unnecessary columns and rename the 'Group' to something more informative
merged_data.drop(columns=['Samples', 'forward-absolute-filepath', 'reverse-absolute-filepath', 'Diet', 'Genotype'], inplace=True)
merged_data.rename(columns={'Group': 'Sample_Group'}, inplace=True)

# For each sample, add a distinctive name based on its group and a number indicating its order within the group
# For example, KO_SD1, KO_SD2, etc.
# First, create a dictionary to keep track of the counts
group_counts = {}

# Define a function to create distinctive names
def create_distinctive_name(row):
    group = row['Sample_Group']
    if group not in group_counts:
        group_counts[group] = 1
    else:
        group_counts[group] += 1
    return f"{group}{group_counts[group]}"

# Apply the function to each row
merged_data['Distinctive_Name'] = merged_data.apply(create_distinctive_name, axis=1)

# Set the new distinctive names as the index
merged_data.set_index('Distinctive_Name', inplace=True)

# Drop the original index and the Sample_Group columns
merged_data.drop(columns=['index', 'Sample_Group'], inplace=True)

# Calculate relative frequencies
relative_data = merged_data.div(merged_data.sum(axis=1), axis=0)

# Define a colormap
colormap = plt.cm.get_cmap('tab20', len(relative_data.columns))

# Assign a color for each taxon
taxa_colors = {taxon: colormap(i) for i, taxon in enumerate(relative_data.columns)}

# Create the stacked bar plot again with the new colors
ax = relative_data.plot(kind='bar', stacked=True, figsize=(20, 10), color=[taxa_colors[taxon] for taxon in relative_data.columns])

# Rename the axis labels
ax.set_xlabel("Sample Name")
ax.set_ylabel("Relative Frequency")

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Move the legend out of the plot
ax.legend(title="Taxa", bbox_to_anchor=(1.04, 0.5), loc="center left")

# Save the figure
plt.savefig(ax, bbox_inches='tight')

ax

