import pandas as pd
# Load GSE126669
df1 = pd.read_csv('GSE126669_tpm_length_scaled_matrix.txt.gz', 
                   sep='\t', index_col=0, compression='gzip')
print(f"GSE126669: {df1.shape} (genes x samples)")
print(df1.index[:5])

# Load GSE109761
df2 = pd.read_csv('GSE109761_processed_normalized_matrix_hs.txt.gz',
                   sep='\t', index_col=0, compression='gzip')
print(f"GSE109761: {df2.shape} (genes x samples)")
print(df2.index[:5])