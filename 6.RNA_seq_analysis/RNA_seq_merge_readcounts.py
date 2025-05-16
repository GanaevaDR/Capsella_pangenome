import pandas as pd
import os

folder_path = './counts_per_gene'
output_file = 'RNA_seq_all_readcounts.csv'

# Initialize an empty DataFrame 
merged_df = pd.DataFrame()

# Loop through all files in the specified folder
for filename in sorted(os.listdir(folder_path)):
    if filename.endswith('.tab'):  
        file_path = os.path.join(folder_path, filename)
        
        # Read the file and extract the first two columns
        df = pd.read_csv(file_path, sep='\t', usecols=[0, 1], header=None)
        
        # Rename columns for clarity
        df.columns = ['gene_id', filename[:-4]]  
        
        # Merge based on 'gene_id'
        if merged_df.empty:
            merged_df = df  
        else:
            merged_df = merged_df.merge(df, on='gene_id', how='outer')

# Sort columns alphabetically, keeping 'gene_id' first
sorted_columns = ['gene_id'] + sorted(merged_df.columns[1:].tolist())
merged_df = merged_df[sorted_columns]

# Save the merged DataFrame to a new tab-delimited file
merged_df.to_csv(output_file, sep='\t', index=False)

print(f'Merged file saved as {output_file}')