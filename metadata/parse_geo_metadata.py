# parse_metadata.py
import re

def parse_geo_metadata(filepath):
    samples = {}
    characteristics = []
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!Sample_geo_accession'):
                ids = re.findall(r'"(GSM\d+)"', line)
            elif line.startswith('!Sample_characteristics_ch1'):
                values = re.findall(r'"([^"]+)"', line)
                characteristics.append(values)
    
    return ids, characteristics

# Parse each dataset
for dataset in ['GSE126669_series_matrix.txt.tar', 
                'GSE109761-GPL18573_series_matrix.txt.tar',
                'GSE180097_series_matrix.txt.tar']:
    ids, chars = parse_geo_metadata(dataset)
    print(f"\n{dataset}: {len(ids)} samples, {len(chars)} characteristic rows")
    for i, char_row in enumerate(chars[:3]):
        print(f"  Row {i}: {char_row[0]} ... {char_row[-1]}")