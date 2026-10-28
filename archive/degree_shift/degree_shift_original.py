import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import random

# =====================
# USER PARAMETERS
# =====================
observed_file = 'humans_g1g2g4_32_v1.csv'
background_file = 'background_proteinatlas_cache.csv'
long_list_file = 'humans_all_DEG_32.csv'
num_iterations = 20
threads = 8

# =====================
# LOAD DATA
# =====================
def load_gene_interactions(file):
    df = pd.read_csv(file, dtype=str).fillna('')
    return {
        row['ENSEMBL_ID']: set(row['KNOWN_INTERACTIONS'].split(';')) if row['KNOWN_INTERACTIONS'] else set()
        for _, row in df.iterrows()
    }

def get_degree_bins(genes, interactions):
    return Counter([len(interactions.get(g, set())) for g in genes])

def shift_distribution(bins, direction='up', aggressive=False):
    shifted = Counter()
    for deg, count in bins.items():
        new_deg = deg + 1 if (direction == 'up' and aggressive) else \
                  max(deg - 1, 0) if (direction == 'down' and aggressive) else deg
        shifted[new_deg] += count

    if not aggressive:
        nonzero = sorted(k for k in shifted if shifted[k] > 0)
        if nonzero:
            d = nonzero[0]
            target = d + 1 if direction == 'up' else max(d - 1, 0)
            shifted[target] += shifted[d]
            del shifted[d]
    return shifted

def sample_genes_from_bins(bins, bg_degrees):
    sampled = []
    pool = defaultdict(list)
    fallback_summary = defaultdict(int)

    for gene, deg in bg_degrees.items():
        pool[deg].append(gene)

    available_degrees = sorted(pool.keys(), reverse=True)

    for target_deg, count in bins.items():
        remaining = count
        current_deg = target_deg

        while remaining > 0:
            options = pool.get(current_deg, [])
            if not options:
                fallback_summary[target_deg] += remaining
                current_deg = max([d for d in available_degrees if d < current_deg], default=None)
                if current_deg is None:
                    print(f"❌ No lower-degree genes available to fulfill bin for degree {target_deg}. {remaining} genes unsampled.")
                    break
                continue

            if len(options) <= remaining:
                sampled += options
                remaining -= len(options)
                pool[current_deg] = []
                current_deg = max([d for d in available_degrees if d < current_deg], default=None)
            else:
                chosen = random.sample(options, remaining)
                sampled += chosen
                pool[current_deg] = list(set(options) - set(chosen))
                remaining = 0

    return sampled, fallback_summary

def build_network_second_order_only(genes, interactions, long_list):
    G = nx.Graph()
    long_set = set(long_list)
    for i, g in enumerate(genes):
        ki_g = interactions.get(g, set())
        for h in genes[i+1:]:
            ki_h = interactions.get(h, set())
            if ki_g & ki_h & long_set:
                G.add_edge(g, h)
    return G

def largest_component_size(G):
    return len(max(nx.connected_components(G), key=len)) if G.number_of_nodes() > 0 else 0

# =====================
# MAIN LOGIC
# =====================
observed_df = pd.read_csv(observed_file, dtype=str).fillna('')
background_df = pd.read_csv(background_file, dtype=str).fillna('')
long_list = pd.read_csv(long_list_file, header=None)[0].tolist()

observed_genes = observed_df['ENSEMBL_ID'].tolist()
background_genes = background_df['ENSEMBL_ID'].tolist()

observed_interactions = load_gene_interactions(observed_file)
background_interactions = load_gene_interactions(background_file)

bg_degrees = {g: len(background_interactions.get(g, set())) for g in background_genes}
observed_bins = get_degree_bins(observed_genes, observed_interactions)

print("📊 Observed degree bin distribution:", dict(observed_bins))

obs_net = build_network_second_order_only(observed_genes, observed_interactions, long_list)
obs_lcc = largest_component_size(obs_net)
print(f"✅ Observed LCC size: {obs_lcc}\n")

scenarios = {
    "Conservative Up": shift_distribution(observed_bins, 'up', aggressive=False),
    "Aggressive Up": shift_distribution(observed_bins, 'up', aggressive=True),
    "Conservative Down": shift_distribution(observed_bins, 'down', aggressive=False),
    "Aggressive Down": shift_distribution(observed_bins, 'down', aggressive=True)
}

results = {}

for label, shifted_bins in scenarios.items():
    print(f"🔄 Running: {label} ...")
    lccs = []
    fallback_counter = defaultdict(int)

    def worker():
        sample, fallback = sample_genes_from_bins(shifted_bins, bg_degrees)
        for bin_deg, fallback_count in fallback.items():
            fallback_counter[bin_deg] += fallback_count
        return largest_component_size(build_network_second_order_only(sample, background_interactions, long_list))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(worker) for _ in range(num_iterations)]
        for f in tqdm(futures):
            lccs.append(f.result())

    results[label] = {
        'mean_lcc': np.mean(lccs),
        'fallback_total': sum(fallback_counter.values()),
        'fallback_bins': len([b for b in fallback_counter if fallback_counter[b] > 0]),
        'mean_fallback_per_bin': np.mean(list(fallback_counter.values())) if fallback_counter else 0
    }
    print(f"✅ Mean LCC: {results[label]['mean_lcc']:.2f}\n")

# =====================
# FINAL OUTPUT
# =====================
print("📊 Summary of Results:")
print(f"Observed LCC size: {obs_lcc}")
for label, stats in results.items():
    print(f"\n📁 {label}:")
    print(f"Mean LCC                  : {stats['mean_lcc']:.2f}")
    # print(f"Total fallback genes      : {stats['fallback_total']}")
    # print(f"Bins with fallback        : {stats['fallback_bins']}")
    # print(f"Mean fallback per bin     : {stats['mean_fallback_per_bin']:.2f}")
