### Simple workflow for converting Ensembl ids into entrez ids ###
### Used for bulk RNA-seq analysis ###

import pandas as pd
import mygene

ensembl = pd.read_excel("RawData_Saherforreassignment.xlsx")
ensembl["gene_id"][1:101]

def convert_ensembl_to_entrez(ensembl_ids, species='mouse'):
    mg = mygene.MyGeneInfo()
    results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='entrezgene', species=species)

    entrez_id_dict = {}
    for result in results:
        ensembl_id = result['query']
        if 'entrezgene' in result:
            if ensembl_id in entrez_id_dict:
                entrez_id_dict[ensembl_id].add(result['entrezgene'])
            else:
                entrez_id_dict[ensembl_id] = {result['entrezgene']}
        else:
            entrez_id_dict[ensembl_id] = {None} # Placeholder for no-matches.

    return entrez_id_dict
    
ensembl_ids = ensembl["gene_id"]

# Convert Ensembl IDs to Entrez IDs
entrez_ids = convert_ensembl_to_entrez(ensembl_ids)

for ensembl_id, entrez_id_set in entrez_ids.items():
  print(f"{ensembl_id}: {entrez_id_set}")


def write_mappings_to_file(mapping_dict, filename):
    with open(filename, 'w') as file:
        for ensembl_id, entrez_id_set in mapping_dict.items():
            # Converting set to comma-separated string for entrez IDs
            entrez_ids_str = ', '.join(map(str, entrez_id_set))
            file.write(f"{ensembl_id}: {entrez_ids_str}\n")

# Assuming your dictionary is named `entrez_ids`
# and you want to write it to 'gene_mappings.txt'
write_mappings_to_file(entrez_ids, 'gene_mappings.txt')
def write_values_to_file(mapping_dict, filename):
    with open(filename, 'w') as file:
        for ensembl_id, entrez_id_set in mapping_dict.items():
            # Converting set to comma-separated string for entrez IDs
            entrez_ids_str = ', '.join(map(str, entrez_id_set))
            file.write(f"{entrez_ids_str}\n")

# Assuming your dictionary is named `entrez_ids`
# and you want to write it to 'gene_mappings.txt'
write_mappings_to_file(entrez_ids, 'entrez_ids.txt')
