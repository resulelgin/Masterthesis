import requests

### This workflow shows an efficient way to scrape KEGG Orthology (KO) API for a given KO number of an enzyme. 
### Workflow also includes a function to filter for taxa of interest. 
### Will be fully integrated to ProHKo workflow in the future.

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

# Example usage.
analyze_ko_number("K00626")
