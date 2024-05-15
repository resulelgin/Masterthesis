### Experimental workflow for obtaining homologous enzyme sequences ###
### Depends on KEGG Orthology API ###
### Knowledgebase-free inferences should be subjected to further bioinformatics confirmations along with experimental verification ###

import requests
import subprocess

def retrieve_species_codes(KO):
    print(f"Species retrieval for KO:{KO} started.")
    url = f'http://rest.kegg.jp/get/ko:{KO}'
    response = requests.get(url)
    species_codes = []
    parsing = False
    for line in response.text.split('\n'):
        if line.startswith('GENES'):
            parsing = True
            # Extract species from the 'GENES' line as well
            line_content = line[len('GENES'):].strip()
            if line_content:  # Check if there's anything to parse on the GENES line
                species_code = line_content.split()[0].split(':')[0].lower()
                species_codes.append(species_code)
        elif parsing:
            if line.startswith(' ') or line.startswith('\t'):
                species_code = line.strip().split()[0].split(':')[0].lower()
                species_codes.append(species_code)
            else:
                break
    print(f"Species retrieval for KO:{KO} ended. {len(species_codes)} species found.")
    return species_codes

def retrieve_gene_list(KO):
    print(f"Gene list retrieval for KO:{KO} started.")
    url = f'http://rest.kegg.jp/get/ko:{KO}'
    response = requests.get(url)
    genes_list = []
    parsing = False
    for line in response.text.split('\n'):
        if line.startswith('GENES'):
            parsing = True
            # Extract species from the 'GENES' line as well
            line_content = line[len('GENES'):].strip()
            if line_content:  # Check if there's anything to parse on the GENES line
                gene_code = line_content.split(':')[1].strip().lower().split(' ')[1]
                genes_list.append(gene_code.strip())
        elif parsing:
            if line.startswith(' ') or line.startswith('\t'):
                gene_code = line.strip().split(':')[1].strip().lower().split(' ')[0]
                genes_list.append(gene_code.strip())
            else:
                break
    print(f"Gene list retrieval for KO:{KO} ended. {len(genes_list)} genes found.")
    return genes_list

def combine_species_genes(species_list, genes_list):
    print(f"Combining species codes and genes started.")
    combined_list = [f'{prefix}:{suffix}' for prefix, suffix in zip(species_list, genes_list)]
    print(f"Combining species codes and genes ended. {len(combined_list)} combined entries created.")
    return combined_list

def fetch_and_write_sequences(ids, filename="sequences.fasta"):
    print(f"Sequence fetching and writing started for {len(ids)} entries.")
    base_url = "https://rest.kegg.jp/get/"
    
    with open(filename, "w") as file:
        for id in ids:
            # Extracting species and gene ID
            species_gene, gene_name = id.split(':')
            gene_id = gene_name.split('(')[0]  # Extracting gene ID
            
            # Constructing the request URL
            request_url = f"{base_url}{species_gene}:{gene_id}"
            response = requests.get(request_url)
            
            if response.status_code != 200:
                print(f"Failed to fetch data for {id}")
                continue
            
            # Parsing the response to extract the amino acid sequence
            parsing = False
            sequence_parts = []
            for line in response.text.split('\n'):
                if line.startswith('AASEQ'):
                    parsing = True
                elif line.startswith('NTSEQ') or line.startswith('///'):
                    parsing = False
                elif parsing and line.startswith('            '):
                    sequence_part = line.strip()
                    sequence_parts.append(sequence_part)
            
            # Joining the sequence parts and creating a FASTA entry
            if sequence_parts:
                seq = ''.join(sequence_parts)
                fasta_header = f">{species_gene.replace(':', '_')}_{gene_id}"
                fasta_sequence = f"{fasta_header}\n{seq}\n"
                file.write(fasta_sequence)
    print(f"Sequence fetching and writing ended. Data written to {filename}.")
                
def ProHKo(KO):
    print(f"Workflow for KO:{KO} started.")
    species_codes = retrieve_species_codes(KO)
    genes_list = retrieve_gene_list(KO)
    combined_list = combine_species_genes(species_codes, genes_list)
    fetch_and_write_sequences(combined_list, filename=f"{KO}_sequences.fasta")
    
    # Run Clustal Omega
    print(f"Running Clustal Omega for MSA...")
    clustalo_command = f"clustalo -i {KO}_sequences.fasta -o {KO}_aligned.aln"
    subprocess.run(clustalo_command, shell=True)
    
    # Build HMM with HMMER
    print(f"Building HMM with hmmbuild...")
    hmmbuild_command = f"hmmbuild {KO}_hmmsequences {KO}_aligned.aln"
    subprocess.run(hmmbuild_command, shell=True)
