# build_metadata.py
import re
import pandas as pd

def parse_geo_metadata(filepath):
    ids = []
    characteristics = []
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!Sample_geo_accession'):
                ids = re.findall(r'"(GSM\d+)"', line)
            elif line.startswith('!Sample_characteristics_ch1'):
                key_match = re.search(r'"([^"]+?):', line)
                if key_match:
                    key = key_match.group(1).strip()
                    values = re.findall(r':\s*([^"]+)"', line)
                    values = [v.strip() for v in values]
                    characteristics.append((key, values))
    
    return ids, characteristics

# Build dataframes for each dataset
datasets = {
    'GSE126669': 'GSE126669_series_matrix.txt.tar',
    'GSE109761': 'GSE109761-GPL18573_series_matrix.txt.tar',
    'GSE180097': 'GSE180097_series_matrix.txt.tar'
}

all_dfs = []

for gse, filepath in datasets.items():
    ids, characteristics = parse_geo_metadata(filepath)
    df = pd.DataFrame({'sample_id': ids})
    df['dataset'] = gse
    for key, values in characteristics:
        if len(values) == len(ids):
            df[key] = values
    all_dfs.append(df)
    print(f"\n{gse} columns: {list(df.columns)}")
    print(df.head(3).to_string())

# Save each dataset separately
all_dfs[0].to_csv('metadata_GSE126669.csv', index=False)
all_dfs[1].to_csv('metadata_GSE109761.csv', index=False)
all_dfs[2].to_csv('metadata_GSE180097.csv', index=False)

print("\nMetadata tables saved.")