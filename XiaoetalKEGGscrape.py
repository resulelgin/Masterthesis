### Repository contains workflow for analysing KEGG Orthology database enzymes with filtering of KEGG species with host association ###
### Host association information of species/strains were obtained from: DOI: 10.1038/nbt.3353 ###

import requests
from bs4 import BeautifulSoup
import webbrowser
import re
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from matplotlib_venn import venn3
import venn
%matplotlib inline

def get_species_from_ko(ko_number):
    url = f'http://rest.kegg.jp/get/ko:{ko_number}'
    response = requests.get(url)
    species_list = set()

    if response.status_code != 200:
        print(f"Error fetching data: {response.status_code}")
        return []

    parsing = False
    for line in response.text.split('\n'):
        if line.startswith('GENES'):
            parsing = True
            continue

        if parsing:
            if line.startswith(' ') or line.startswith('\t'):
                species_code = line.split()[0].split(':')[0].lower()
                species_list.add(species_code)
            else:
                break

    return list(species_list)

def get_full_species_name_from_web(org_code):
    url = 'https://rest.kegg.jp/list/organism'
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Error fetching data for {org_code}: {response.status_code}")
        return org_code
    result = []
    for line in response.text.split('\n'):
        if '\t' in line:
            abv = line.split('\t')[1]
            full = line.split('\t')[2]
            family = line.split('\t')[3]
            for code in org_code:
                if code == abv:
                    result.append(f"{full}, {family}")
    return result
    
def filter_for_bacteria(species_list):
    return [species for species in species_list if 'Bacteria' in species]

def export_list_to_txt(filename, your_list):
    with open(filename, 'w') as file:
        for item in your_list:
            file.write(f"{item}\n")

def analyze_ko_number(ko_number):
    species_codes = get_species_from_ko(ko_number)
    full_species_names = get_full_species_name_from_web(species_codes)
    bacteria_species = filter_for_bacteria(full_species_names)
    export_list_to_txt(f"{ko_number}_bacteria.txt", bacteria_species)
    return bacteria_species

ko_number = "K07518"
result = analyze_ko_number(ko_number)
print(result)

bdh1 = pd.read_table("K00019_bacteria.txt", header=None)
acat = pd.read_table("K00626_bacteria.txt", header=None)
oxct = pd.read_table("K01027_bacteria.txt", header=None)
ato = pd.read_table("K01034_bacteria.txt", header=None)
hmgcl = pd.read_table("K01640_bacteria.txt", header=None)
hmgcs2 = pd.read_table("K01641_bacteria.txt", header=None)
kama = pd.read_table("K01843_bacteria.txt", header=None)
kamd = pd.read_table("K01844_bacteria.txt", header=None)
kdd = pd.read_table("K18012_bacteria.txt", header=None)
kce = pd.read_table("K18013_bacteria.txt", header=None)
fadj = pd.read_table("K01782_bacteria.txt", header=None)
fadb = pd.read_table("K01825_bacteria.txt", header=None)
phac = pd.read_table("K03821_bacteria.txt", header=None)
phae = pd.read_table("K22881_bacteria.txt", header=None)
phaz = pd.read_table("K05973_bacteria.txt", header=None)
hdep = pd.read_table("K07518_bacteria.txt", header=None)

set_bdh1 = set(bdh1[0])
set_acat = set(acat[0])
set_oxct = set(oxct[0])
set_ato = set(ato[0])
set_hmgcl = set(hmgcl[0])
set_hmgcs2 = set(hmgcs2[0])
set_kama = set(kama[0])
set_kamd = set(kamd[0])
set_kdd = set(kdd[0])
set_kce = set(kce[0])
set_fadj = set(fadj[0])
set_fadb = set(fadb[0])
set_phac = set(phac[0])
set_phae = set(phae[0])
set_phaz = set(phaz[0])
set_hdep = set(hdep[0])

# Define your set of interest.
set_kce & set_kama & set_kamd & set_kdd

mgs = pd.read_excel("MmMGS_summary.xlsx", header = 1)
mgs_species = mgs["Main species gene annotation"]
mgs_species

table184 = pd.read_csv("184sample.NR-tax.relativeAbun.table.species", sep='\t')
species_table = table184["Unnamed: 0"]
species_table

xiao_species = pd.concat([mgs_species, species_table], ignore_index=True)
xiao_species

with open("xiao_metagenome_species.txt", 'w') as file:
    for item in xiao_species:
        file.write(f"{item}\n")

def filter_ko_with_mgs(ko, mgs=xiao_species):
    ko_series = ko[0]
    result = []
    for sp in mgs:
        for element in ko_series:
            if sp in element:
                result.append(sp)
    return result

def encapsulate_filter_ko_with_mgs(gene):
    filtered_gene = filter_ko_with_mgs(gene)
    return filtered_gene

filtered_acat = encapsulate_filter_ko_with_mgs(acat)
filtered_bdh1 = encapsulate_filter_ko_with_mgs(bdh1)
filtered_oxct = encapsulate_filter_ko_with_mgs(oxct)
filtered_ato = encapsulate_filter_ko_with_mgs(ato)
filtered_hmgcl = encapsulate_filter_ko_with_mgs(hmgcl)
filtered_hmgcs2 = encapsulate_filter_ko_with_mgs(hmgcs2)

set_filtered_bdh1 = set(filtered_bdh1)
set_filtered_acat = set(filtered_acat)
set_filtered_oxct = set(filtered_oxct)
set_filtered_ato = set(filtered_ato)
set_filtered_hmgcl = set(filtered_hmgcl)
set_filtered_hmgcs2 = set(filtered_hmgcs2)

set1 = set_filtered_bdh1
set2 = set_filtered_acat
set3 = set_filtered_oxct
set4 = set_filtered_ato
set5 = set_filtered_hmgcl
set6 = set_filtered_hmgcs2

# Create a figure with subplots
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

# Plot Venn diagrams for the first three sets
venn3([set1, set2, set6], set_labels=('BDH1', 'ACAT', 'HMGCS2'), ax=axes[0, 0])
venn3([set1, set2, set5], set_labels=('BDH1', 'ACAT', 'HMGCL'), ax=axes[0, 1])
venn3([set1, set5, set6], set_labels=('BDH1', 'HMGCL', 'HMGCS2'), ax=axes[1, 0])

# Create an empty axis for the fourth set
axes[1, 1].axis('off')

# Show the plot
plt.suptitle("Xiao-filtered KEGG Genes")
plt.show()

labels = venn.get_labels([set1, set2, set5, set6], fill=['number'])
fig, ax = venn.venn4(labels, names=['BDH1', 'ACAT', 'HMGCL', 'HMGCS2'])
plt.title("Xiao-filtered KEGG Genes")
fig.show()

labels = venn.get_labels([set1, set2, set4, set5, set6], fill=['number'])
fig, ax = venn.venn5(labels, names=['BDH1', 'ACAT', 'ATO', 'HMGCL', 'HMGCS2'])
plt.title("Xiao-filtered KEGG Genes")
fig.show()

table_16s = pd.read_csv("level7_taxa.csv", sep=",", header=0)
table_16s

species_16s = []
for names in table_16s.columns:
    species_16s.append(names)

xiao_formatted = []
for element in (xiao_species):
    element_split = element.split(" ")[0] + "_" + element.split(" ")[1]
    xiao_formatted.append(element_split)
xiao_formatted

overlaps = []
for i in xiao_formatted:
    for j in species_16s:
        if i in j:
            overlaps.append(j)
set(overlaps)

for sp in xiao_species:
    if "Lactobacillus" in sp:
        print(sp)

for sp in species_16s:
    for element in (set_filtered_oxct):
        element_split = element.split(" ")[0] + "_" + element.split(" ")[1]
        if element_split in sp:
            print(sp)

      
