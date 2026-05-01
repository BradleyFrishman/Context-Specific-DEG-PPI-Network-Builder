import pandas as pd

# Load CSV
df = pd.read_csv('new human genes 32 scraped.csv', encoding='utf-8-sig')

# Clean columns (strip and uppercase)
df.columns = [col.encode('ascii', errors='ignore').decode().strip().upper() for col in df.columns]
print(df.columns.tolist())  # confirm columns

# Build set of valid ENSEMBL_IDs from column A
valid_ids = set(df['ENSEMBL_ID'].dropna().astype(str))

# Function to refine known interactions
def refine_interactions(interactions):
    if pd.isna(interactions):
        return ''
    ids = [gene_id.strip() for gene_id in interactions.split(';')]
    filtered_ids = [gene_id for gene_id in ids if gene_id in valid_ids]
    return ';'.join(filtered_ids)

# Get index of KNOWN_INTERACTIONS column (B)
known_interactions_idx = df.columns.get_loc('KNOWN_INTERACTIONS')

# Create new column with refined interactions
df.insert(known_interactions_idx + 1, 'REFINED_KNOWN_INTERACTIONS', df['KNOWN_INTERACTIONS'].apply(refine_interactions))

# Save over the existing file
df.to_csv('new human genes 32 scraped.csv', index=False, encoding='utf-8-sig')

print("✅ Refinement complete — file updated with new 'REFINED_KNOWN_INTERACTIONS' column.")
