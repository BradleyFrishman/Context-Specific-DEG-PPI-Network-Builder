import pandas as pd
import os
import sys

long_list_file = 'long_list_32.csv'  # adjust path if needed

# Load long list of ENSEMBL IDs
long_list_df = pd.read_csv(long_list_file, dtype=str)
long_list_ids = set(long_list_df['ENSEMBL_ID'].dropna().str.strip())

def count_genes(cell):
    genes = [g.strip() for g in cell.split(';') if g.strip()]
    return len(genes)

def filter_long_list(cell):
    genes = [g.strip() for g in cell.split(';') if g.strip()]
    matching_genes = [g for g in genes if g in long_list_ids]
    return ';'.join(matching_genes)

def interaction_set(cell):
    return set(cell.split(';')) if cell.strip() else set()

# Fixed file path - adjust if needed
file_path = 'protein_names_processes_and_interactions_with_connections.csv'

if not os.path.isfile(file_path):
    print(f"Error: File '{file_path}' not found.")
    sys.exit(1)

print(f"▶️ Processing {file_path}...")

# Load CSV
df = pd.read_csv(file_path, dtype=str).fillna('')

# Add counts for KNOWN_INTERACTIONS and REFINED_KNOWN_INTERACTIONS if present
if 'KNOWN_INTERACTIONS' in df.columns:
    df['KNOWN_INTERACTIONS_COUNT'] = df['KNOWN_INTERACTIONS'].apply(count_genes)
else:
    print("Warning: 'KNOWN_INTERACTIONS' column not found.")

if 'REFINED_KNOWN_INTERACTIONS' in df.columns:
    df['REFINED_KNOWN_INTERACTIONS_COUNT'] = df['REFINED_KNOWN_INTERACTIONS'].apply(count_genes)
else:
    print("Warning: 'REFINED_KNOWN_INTERACTIONS' column not found.")

# Create KNOWN_INTERACTIONS_W_LONG_LIST and count, if possible
if 'KNOWN_INTERACTIONS' in df.columns:
    df['KNOWN_INTERACTIONS_W_LONG_LIST'] = df['KNOWN_INTERACTIONS'].apply(filter_long_list)
    df['KNOWN_INTERACTIONS_W_LONG_LIST_COUNT'] = df['KNOWN_INTERACTIONS_W_LONG_LIST'].apply(count_genes)
else:
    print("Warning: 'KNOWN_INTERACTIONS' column not found; cannot create 'KNOWN_INTERACTIONS_W_LONG_LIST'.")
    df['KNOWN_INTERACTIONS_W_LONG_LIST'] = ''
    df['KNOWN_INTERACTIONS_W_LONG_LIST_COUNT'] = 0

# Compute CONNECTED_GENES_W_LONG_LIST and count
df['CONNECTED_GENES_W_LONG_LIST'] = ''

total_genes = len(df)

for idx_i in range(total_genes):
    gene_i = df.loc[idx_i, 'ENSEMBL_ID']
    interactions_i = interaction_set(df.loc[idx_i, 'KNOWN_INTERACTIONS_W_LONG_LIST'])

    connected_genes = []

    for idx_j in range(idx_i + 1, total_genes):
        gene_j = df.loc[idx_j, 'ENSEMBL_ID']
        interactions_j = interaction_set(df.loc[idx_j, 'KNOWN_INTERACTIONS_W_LONG_LIST'])

        if interactions_i & interactions_j:
            connected_genes.append(gene_j)

    df.loc[idx_i, 'CONNECTED_GENES_W_LONG_LIST'] = ';'.join(connected_genes)

    if (idx_i + 1) % 10 == 0 or idx_i == total_genes - 1:
        print(f'  Processed {idx_i + 1}/{total_genes} genes ({(idx_i + 1)/total_genes:.1%})')

df['CONNECTED_GENES_W_LONG_LIST_COUNT'] = df['CONNECTED_GENES_W_LONG_LIST'].apply(count_genes)

# Save updated CSV overwriting original
df.to_csv(file_path, index=False)
print(f'✅ Saved {file_path}.')
print('\n🎉 File processed and updated successfully.')
