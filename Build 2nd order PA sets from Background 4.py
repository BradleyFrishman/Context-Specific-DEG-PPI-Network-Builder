import os
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
import networkx as nx
from collections import defaultdict, Counter
import numpy as np
from tqdm import tqdm
import random

# ===============================
# 🎛️ USER PARAMETERS
# ===============================

background_file = 'background_ensg.txt'                  
observed_file   = 'humans_G3_41_v1.csv'       
long_list_file  = 'humans_all_DEG_41.csv'                
output_folder   = 'July 19 Random Files for G3 41º (Removing 0 KI, NOT Filtering out observed from Background)'                 
num_random_files = 20                                     
max_threads = 20                                           
num_bins = 20

# ===============================
# 📦 PREPARATION: Load Input Data
# ===============================

os.makedirs(output_folder, exist_ok=True)

print("\n📖 Loading input files...")
with open(background_file) as f:
    background_genes = [line.strip() for line in f if line.strip()]

observed_df = pd.read_csv(observed_file, dtype=str).fillna('')
observed_genes = observed_df['ENSEMBL_ID'].tolist()

long_list_df = pd.read_csv(long_list_file, dtype=str)
long_list_ids = set(long_list_df['ENSEMBL_ID'].dropna().str.strip())

print(f"  • Background genes loaded: {len(background_genes)}")
print(f"  • Observed genes loaded: {len(observed_genes)}")
print(f"  • Long list interactors loaded: {len(long_list_ids)}")

# ===============================
# 🌐 SCRAPING PROTEINATLAS INTERACTIONS (WITH CACHE)
# ===============================

def scrape_proteinatlas(ensembl_id):
    url = f"https://www.proteinatlas.org/{ensembl_id}/structure+interaction"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code != 200:
            return set()
        html = r.text

        functions_content, start_pos = [], 0
        while True:
            func_index = html.find("function()", start_pos)
            if func_index == -1:
                break
            brace_start = html.find("{", func_index)
            open_braces, i = 1, brace_start + 1
            while i < len(html) and open_braces > 0:
                if html[i] == '{': open_braces += 1
                elif html[i] == '}': open_braces -= 1
                i += 1
            functions_content.append(html[brace_start + 1: i - 1])
            start_pos = i

        combined_str, gene_ids, pos = " ".join(functions_content), [], 0
        while True:
            pos = combined_str.find('"id":"int', pos)
            if pos == -1: break
            start = pos + len('"id":"int')
            end_pos = start
            while end_pos < len(combined_str) and combined_str[end_pos] not in ['"', ',']:
                end_pos += 1
            gene_ids.append(combined_str[start:end_pos])
            pos = end_pos

        interactors = set(gene_ids)
        interactors.discard(ensembl_id)
        return interactors

    except:
        return set()

cache_file = 'background_proteinatlas_cache.csv'

if os.path.exists(cache_file):
    print("📄 Loading cached interaction data...")
    cache_df = pd.read_csv(cache_file, dtype=str).fillna('')
    proteinatlas_map = {
        row['ENSEMBL_ID']: set(row['KNOWN_INTERACTIONS'].split(';')) if row['KNOWN_INTERACTIONS'] else set()
        for _, row in cache_df.iterrows()
    }
else:
    print("🔍 Scraping ProteinAtlas for background genes...")
    proteinatlas_map = {}
    def scrape_worker(gene):
        return gene, scrape_proteinatlas(gene)

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = {executor.submit(scrape_worker, g): g for g in background_genes}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Scraping", unit="gene"):
            gene, interactions = future.result()
            proteinatlas_map[gene] = interactions

    cache_df = pd.DataFrame([
        {'ENSEMBL_ID': gene, 'KNOWN_INTERACTIONS': ';'.join(sorted(list(interactors)))}
        for gene, interactors in proteinatlas_map.items()
    ])
    cache_df.to_csv(cache_file, index=False)
    print(f"💾 Cache saved: {cache_file}")

# ===============================
# 📊 OBSERVED GENE FILTERING
# ===============================

found_genes = [g for g in observed_genes if g in proteinatlas_map]
print(f"\n✅ Genes from observed file FOUND in Protein Atlas: {len(found_genes)} / {len(observed_genes)}")

# 👇 Optional: exclude genes with 0 known interactions too
found_genes_with_interactions = [g for g in found_genes if len(proteinatlas_map[g]) > 0]
print(f"🔍 Genes with ≥1 known Protein Atlas interactions: {len(found_genes_with_interactions)} / {len(observed_genes)}")

missing_genes = [g for g in observed_genes if g not in proteinatlas_map]
print(f"❌ Genes from observed file MISSING in Protein Atlas: {len(missing_genes)} / {len(observed_genes)}")

# 👇 Set whether you want to use only genes with ≥1 interaction or all found genes
use_zero_interactors = False  # set to True if you want to keep genes with 0 interactions

if use_zero_interactors:
    final_genes_to_use = found_genes
else:
    final_genes_to_use = found_genes_with_interactions

observed_df = observed_df[observed_df['ENSEMBL_ID'].isin(final_genes_to_use)]
observed_genes = observed_df['ENSEMBL_ID'].tolist()
random_sample_size = len(observed_genes)

print(f"🎛️ Final random sample size per file: {random_sample_size}\n")


observed_df = observed_df[observed_df['ENSEMBL_ID'].isin(found_genes)]
observed_genes = observed_df['ENSEMBL_ID'].tolist()
random_sample_size = len(observed_genes)

print(f"🎛️ Random sample size per file (matching found observed genes): {random_sample_size}\n")

# Remove observed genes from background before binning
#background_genes = [g for g in background_genes if g not in observed_genes]

# ===============================
# 📊 DEGREE BINNING
# ===============================

# Degree mapping for background and observed separately
degree_map_background = {gene: len(proteinatlas_map.get(gene, set())) for gene in background_genes}
degree_map_obs = {gene: len(proteinatlas_map.get(gene, set())) for gene in observed_genes}

obs_degrees = [degree_map_obs.get(g, 0) for g in observed_genes]

# Bin edges based on observed degrees
bin_edges = np.quantile(obs_degrees, np.linspace(0, 1, num_bins + 1))

def assign_bin(degree):
    for i in range(num_bins):
        if degree <= bin_edges[i+1]:
            return i
    return num_bins - 1

# Bin observed genes
observed_bins = [assign_bin(d) for d in obs_degrees]
bin_counts = Counter(observed_bins)

print(f"🎛️ Observed degree bin counts (using {num_bins} bins): {dict(bin_counts)}\n")

# Bin background genes
background_bins = defaultdict(list)
for gene in background_genes:
    b = assign_bin(degree_map_background.get(gene, 0))
    background_bins[b].append(gene)

# ===============================
# 🔍 OBSERVED NETWORK BUILDING
# ===============================

print("🔍 Building observed gene network...")

observed_df['KNOWN_INTERACTIONS_W_LONG_LIST'] = observed_df['ENSEMBL_ID'].apply(
    lambda g: ';'.join(sorted(proteinatlas_map.get(g, set()) & long_list_ids))
)

known_int_list = observed_df['KNOWN_INTERACTIONS_W_LONG_LIST'].tolist()
total_genes = len(observed_df)
connected_lists = [''] * total_genes

for idx_i in range(total_genes):
    int_i = set(known_int_list[idx_i].split(';')) if known_int_list[idx_i].strip() else set()
    connected = []
    for idx_j in range(idx_i + 1, total_genes):
        int_j = set(known_int_list[idx_j].split(';')) if known_int_list[idx_j].strip() else set()
        if int_i & int_j:
            connected.append(observed_df.iloc[idx_j]['ENSEMBL_ID'])  # <-- FIXED HERE
    connected_lists[idx_i] = ';'.join(connected)

observed_df['CONNECTED_GENES_W_LONG_LIST'] = connected_lists
# Build observed network graph
G_obs = nx.Graph()
G_obs.add_nodes_from(observed_df['ENSEMBL_ID'])
for idx, row in observed_df.iterrows():
    gene = row['ENSEMBL_ID']
    for cgene in row['CONNECTED_GENES_W_LONG_LIST'].split(';'):
        if cgene.strip():
            G_obs.add_edge(gene, cgene)

# Get connected components info
components = list(nx.connected_components(G_obs))
num_components_2plus = sum(1 for c in components if len(c) >= 2)
component_sizes = sorted([len(c) for c in components if len(c) >= 2], reverse=True)

print(f"🔗 Observed network: {num_components_2plus} connected components with ≥2 genes")
print(f"📏 Component sizes: {component_sizes}")

# Largest connected component
if components:
    largest_cc_obs = max(components, key=len)
    subgraph_obs = G_obs.subgraph(largest_cc_obs)
    print(f"✅ Observed network: Largest connected component size = {subgraph_obs.number_of_nodes()}, edges = {subgraph_obs.number_of_edges()}\n")
    observed_lcc_size = subgraph_obs.number_of_nodes()
    observed_edge_count = subgraph_obs.number_of_edges()

else:
    print("⚠️ No connected components found in observed network.\n")

# ===============================
# 🌟 RANDOM FILE GENERATION + NETWORK BUILDING
# ===============================
random_lcc_sizes = []
random_edge_counts = []

def sample_genes_by_bins():
    sampled = []
    for b in range(num_bins):
        count = bin_counts.get(b, 0)
        pool = background_bins[b]
        if len(pool) >= count:
            sampled.extend(random.sample(pool, count))
        else:
            sampled.extend(random.choices(pool, k=count))
    return sampled

for i in range(1, num_random_files + 1):
    sampled_genes = sample_genes_by_bins()
    df_random = pd.DataFrame({'ENSEMBL_ID': sampled_genes}).reset_index(drop=True)

    df_random['KNOWN_INTERACTIONS_W_LONG_LIST'] = df_random['ENSEMBL_ID'].apply(
        lambda g: ';'.join(sorted(proteinatlas_map.get(g, set()) & long_list_ids))
    )

    known_int_list, total_genes = df_random['KNOWN_INTERACTIONS_W_LONG_LIST'].tolist(), len(df_random)
    connected_lists = [''] * total_genes

    for idx_i in range(total_genes):
        int_i = set(known_int_list[idx_i].split(';')) if known_int_list[idx_i].strip() else set()
        connected = []
        for idx_j in range(idx_i + 1, total_genes):
            int_j = set(known_int_list[idx_j].split(';')) if known_int_list[idx_j].strip() else set()
            if int_i & int_j:
                connected.append(df_random.iloc[idx_j]['ENSEMBL_ID'])  # <-- FIXED HERE
        connected_lists[idx_i] = ';'.join(connected)

    df_random['CONNECTED_GENES_W_LONG_LIST'] = connected_lists

    G = nx.Graph()
    G.add_nodes_from(df_random['ENSEMBL_ID'])
    for idx, row in df_random.iterrows():
        gene = row['ENSEMBL_ID']
        for cgene in row['CONNECTED_GENES_W_LONG_LIST'].split(';'):
            if cgene.strip():
                G.add_edge(gene, cgene)

    components = list(nx.connected_components(G))
    num_components_2plus = sum(1 for c in components if len(c) >= 2)
    component_sizes = sorted([len(c) for c in components if len(c) >= 2], reverse=True)

    print(f"🔗 Random file {i:2}: {num_components_2plus} connected components with ≥2 genes")
    print(f"📏 Component sizes: {component_sizes}")

    largest_cc = max(components, key=len) if len(components) > 0 else set()
    subgraph = G.subgraph(largest_cc)
    # Store LCC stats
    random_lcc_sizes.append(subgraph.number_of_nodes())
    random_edge_counts.append(subgraph.number_of_edges())

    print(f"Random file {i:2}: LCC size = {subgraph.number_of_nodes()}, edges = {subgraph.number_of_edges()}")
    
    save_path = os.path.join(output_folder, f'random_g3_j16__20bins_41_{i}.csv')
    df_random[['ENSEMBL_ID', 'KNOWN_INTERACTIONS_W_LONG_LIST', 'CONNECTED_GENES_W_LONG_LIST']].to_csv(save_path, index=False)

print("\n✅ All random files generated and networks built!")

from scipy import stats

print("\n📊 1-Sample T-Tests: Random vs Observed Network Statistics")

# LCC size
t_stat_lcc, p_val_lcc = stats.ttest_1samp(random_lcc_sizes, observed_lcc_size)
print("\n🔹 LCC Size:")
print(f"  • Mean of random LCC sizes: {np.mean(random_lcc_sizes):.2f}")
print(f"  • Observed LCC size: {observed_lcc_size}")
print(f"  • T-statistic: {t_stat_lcc:.3f}")
print(f"  • P-value: {p_val_lcc:.5f}")

# Edge count
t_stat_edges, p_val_edges = stats.ttest_1samp(random_edge_counts, observed_edge_count)
print("\n🔹 LCC Edge Count:")
print(f"  • Mean of random LCC edge counts: {np.mean(random_edge_counts):.2f}")
print(f"  • Observed LCC edge count: {observed_edge_count}")
print(f"  • T-statistic: {t_stat_edges:.3f}")
print(f"  • P-value: {p_val_edges:.5f}")
print(random_lcc_sizes)
print(observed_lcc_size)

