# filter_patient_samples.py
import pandas as pd

# Load all three
df_126669 = pd.read_csv('metadata_GSE126669.csv')
df_109761 = pd.read_csv('metadata_GSE109761.csv')
df_180097 = pd.read_csv('metadata_GSE180097.csv')

# GSE126669 - filter patient (Br61 only)
patient_126669 = df_126669[df_126669['model'] == 'Br61']
print(f"GSE126669 patient samples: {len(patient_126669)}")
print(patient_126669[['sample_id','model','hypoxia']].to_string())

# GSE109761 - filter patient CTCs only (no WBC, no REF)
patient_109761 = df_109761[
    (df_109761['origin'] == 'patient') & 
    (df_109761['cell type'] == 'CTC')
]
print(f"\nGSE109761 patient CTC samples: {len(patient_109761)}")
print(patient_109761[['sample_id','donor','sample type']].value_counts('sample type'))

# GSE180097 - filter patient only
patient_180097 = df_180097[df_180097['origin'] == 'patient']
print(f"\nGSE180097 patient samples: {len(patient_180097)}")
print(patient_180097[['sample_id','sample type','timepoint']].value_counts('timepoint'))

# Save filtered lists
patient_126669.to_csv('patient_samples_GSE126669.csv', index=False)
patient_109761.to_csv('patient_samples_GSE109761.csv', index=False)
patient_180097.to_csv('patient_samples_GSE180097.csv', index=False)

print("\nFiltered patient sample lists saved.")