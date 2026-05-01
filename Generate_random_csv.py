import pandas as pd
import random
import os

# Load the gene list from the CSV file
input_file = 'humans_all.csv'
df = pd.read_csv(input_file)

# Extract gene list from 'ENSEMBL_ID' column
gene_list = df['ENSEMBL_ID'].dropna().tolist()

# Number of random sets and genes per set
num_files = 20
genes_per_file = 247

# Check that genes_per_file isn't larger than total genes
if genes_per_file > len(gene_list):
    raise ValueError("Not enough genes in the input list to sample 247 genes.")

# Define output folder
output_folder = 'Random Files Attempt 2'

# Make sure the folder exists
os.makedirs(output_folder, exist_ok=True)

# Generate and save 20 random CSV files into the specified folder
for i in range(1, num_files + 1):
    random_genes = random.sample(gene_list, genes_per_file)
    random_df = pd.DataFrame(random_genes, columns=['ENSEMBL_ID'])
    output_file = os.path.join(output_folder, f'random_v2_{i}.csv')
    random_df.to_csv(output_file, index=False)

print("Done! 20 random CSV files created in the 'Random Files Attempt 2' folder.")
