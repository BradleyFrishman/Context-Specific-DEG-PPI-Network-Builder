import pandas as pd

# Load CSV file
df = pd.read_csv('long_list_32.csv')

# Ensure 'KNOWN_INTERACTIONS' is treated as string and fill NaNs with empty strings
df['KNOWN_INTERACTIONS'] = df['KNOWN_INTERACTIONS'].fillna('').astype(str)

# Initialize new column for connected genes
df['CONNECTED_GENES'] = ''

total_genes = len(df)

# Iterate through each gene
for i in range(total_genes):
    gene_i = df.loc[i, 'ENSEMBL_ID']
    interactions_i = set(df.loc[i, 'KNOWN_INTERACTIONS'].split(';')) if df.loc[i, 'KNOWN_INTERACTIONS'] else set()

    connected_genes = []

    for j in range(i + 1, total_genes):
        gene_j = df.loc[j, 'ENSEMBL_ID']
        interactions_j = set(df.loc[j, 'KNOWN_INTERACTIONS'].split(';')) if df.loc[j, 'KNOWN_INTERACTIONS'] else set()

        # Check if interacting sets intersect
        if interactions_i & interactions_j:
            connected_genes.append(gene_j)

    # Save connected genes (semicolon-separated string) to the new column
    df.loc[i, 'CONNECTED_GENES'] = ';'.join(connected_genes)

    # Progress printout every 10 rows
    if (i + 1) % 10 == 0 or i == total_genes - 1:
        print(f"Processed {i + 1} of {total_genes} genes ({(i + 1) / total_genes:.1%})")

# Overwrite original CSV file with new column
df.to_csv('long_list_32.csv', index=False)

print("Done. Updated 'long_list_32.csv' with 'CONNECTED_GENES' column.")
