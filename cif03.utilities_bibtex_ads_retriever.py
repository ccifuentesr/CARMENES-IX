import requests
import re
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.customization import convert_to_unicode

def fetch_bibtex(bibcode, api_key):
    url = f'https://api.adsabs.harvard.edu/v1/export/bibtex'
    headers = {
        'Authorization': f'Bearer {api_key}',
        'Content-Type': 'application/json'
    }
    payload = {
        'bibcode': [bibcode]
    }
    response = requests.post(url, headers=headers, json=payload)
    
    if response.status_code == 200:
        return response.json()['export']
    else:
        print(f"Error fetching bibtex for {bibcode}: {response.text}")
        return None

def customize_bibtex(bibtex_str):
    parser = BibTexParser()
    parser.customization = convert_to_unicode
    bib_database = bibtexparser.loads(bibtex_str, parser=parser)

    equivalence_dict = {}

    for entry in bib_database.entries:
        # Check if 'year' is present in the entry
        if 'year' in entry:
            # Determine the key based on available fields
            if 'author' in entry:
                author = entry['author'].split(' and ')[0].strip('{}')
                surname = author.split(',')[0]
                identifier = f"{surname[:3]}{entry['year'][-2:]}"
            elif 'title' in entry:
                identifier = entry['title'].split()[0].lower()[:5] + entry['year'][-2:]
            else:
                identifier = entry['ID']  # Use the original ID if no suitable fields are found

            # Update equivalence_dict
            equivalence_dict[entry['ID']] = identifier
            entry['ID'] = identifier
        else:
            print(f"Warning: Entry {entry['ID']} does not have 'year' field. Skipping.")

    return bib_database, equivalence_dict




def main(input_file, output_file, equivalence_file, api_key):
    with open(input_file, 'r') as infile:
        bibcodes = [line.strip() for line in infile if line.strip()]
    
    bibtex_entries = []
    equivalence_dict = {}

    for bibcode in bibcodes:
        bibtex = fetch_bibtex(bibcode, api_key)
        if bibtex:
            customized_bib, eq_dict = customize_bibtex(bibtex)
            bibtex_entries.append((customized_bib, eq_dict))

    # Filter out entries without 'author' fields
    bibtex_entries = [(bib, eq) for bib, eq in bibtex_entries if bib.entries and 'author' in bib.entries[0]]

    # Sort bibtex_entries by author surname
    bibtex_entries.sort(key=lambda x: x[0].entries[0]['author'].split(' and ')[0].split(',')[0])

    with open(output_file, 'w') as outfile:
        for bib_db, _ in bibtex_entries:
            for entry in bib_db.entries:
                bibtexparser.dump(bib_db, outfile)
                outfile.write('\n\n')

    with open(equivalence_file, 'w') as eqfile:
        for _, eq_dict in bibtex_entries:
            for original, new in eq_dict.items():
                eqfile.write(f"{original}, {new}\n")

if __name__ == "__main__":
    input_file = 'Data/bibcodes.txt'
    output_file = 'Output/bibtex_entries.txt'
    equivalence_file = 'Output/bibtex_abbreviations.txt'
    api_key = 'MKlR5MLEtF7ytw7DxdhGdNyl62w24w7gwUosVAHI'  # Replace with your ADS API key
    main(input_file, output_file, equivalence_file, api_key)


