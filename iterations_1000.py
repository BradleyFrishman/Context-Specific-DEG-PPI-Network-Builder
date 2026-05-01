import pandas as pd
import requests
import random
import concurrent.futures
from collections import defaultdict
import time

# CONFIG
humans_file = 'humans_all.csv'
long_list_file = 'long_list_32.csv'
sample_size = 247
iterations = 1000
num_workers = 16  # Adjust based on your CPU
observed_value = 183  # Your observed largest network size

# Load gene lists
humans_df = pd.read_csv(humans_file, dtype=str)
humans_genes = humans_df['ENSEMBL_ID'].dropna().str.strip().tolist()
long_list_ids = set(pd.read_csv(long_list_file, dtype=str)['ENSEMBL_ID'].dropna().str.strip())

# Scraper cache - global, filled once upfront
scrape_cache = {}

def scrape_interacting_genes(ensembl_id):
    if ensembl_id in scrape_cache:
        return scrape_cache[ensembl_id]
    try:
        url = f"https://www.proteinatlas.org/{ensembl_id}/structure+interaction"
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            scrape_cache[ensembl_id] = []
            return []
        html = response.text
        contents, start = [], 0
        while True:
            idx = html.find('"id":"int', start)
            if idx == -1:
                break
            start = idx + 9
            end = start
            while end < len(html) and html[end] not in ['"', ',']:
                end += 1
            gene = html[start:end]
            contents.append(gene)
            start = end
        interactions = list(set(contents) - {ensembl_id})
        scrape_cache[ensembl_id] = interactions
        return interactions
    except:
        scrape_cache[ensembl_id] = []
        return []

def build_neighbors(sample_genes):
    known_map = {gene: scrape_cache.get(gene, []) for gene in sample_genes}

    refined_map = {
        gene: [g for g in interactors if g in sample_genes]
        for gene, interactors in known_map.items()
    }

    long_list_map = {
        gene: [g for g in interactors if g in long_list_ids]
        for gene, interactors in known_map.items()
    }

    neighbor_map = defaultdict(set)

    for gene, first_neighbors in refined_map.items():
        neighbor_map[gene].update(first_neighbors)

    for i, gene_i in enumerate(sample_genes):
        for gene_j in sample_genes[i+1:]:
            if set(long_list_map[gene_i]) & set(long_list_map[gene_j]):
                neighbor_map[gene_i].add(gene_j)
                neighbor_map[gene_j].add(gene_i)

    return neighbor_map

def expand_network(seed, neighbor_map, available):
    network, frontier = {seed}, {seed}
    while frontier:
        new_genes = set()
        for gene in frontier:
            new_neighbors = neighbor_map[gene] & available - network
            new_genes.update(new_neighbors)
        if not new_genes:
            break
        network.update(new_genes)
        frontier = new_genes
    return network

def build_networks(sample_genes, neighbor_map):
    available = set(sample_genes)
    sizes = []
    while available:
        seed = next(iter(available))
        network = expand_network(seed, neighbor_map, available)
        sizes.append(len(network))
        available -= network
    return sizes

def single_iteration(seed):
    random.seed(seed)
    sample_genes = random.sample(humans_genes, sample_size)
    neighbor_map = build_neighbors(sample_genes)
    sizes = build_networks(sample_genes, neighbor_map)
    return max(sizes) if sizes else 0

def main():
    print("▶️ Pre-scraping all genes. This may take a while...")
    start_scrape = time.time()
    total_genes = len(humans_genes)
    batch_size = max(1, total_genes // 20)  # print progress every 5%
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i, gene in enumerate(humans_genes, 1):
            futures.append(executor.submit(scrape_interacting_genes, gene))
            if i % batch_size == 0 or i == total_genes:
                elapsed = time.time() - start_scrape
                percent = i / total_genes * 100
                est_total = elapsed / (i / total_genes)
                est_remain = est_total - elapsed
                print(f"  Pre-scraped {i}/{total_genes} genes ({percent:.1f}%), "
                      f"Elapsed: {elapsed:.1f}s, Est. remaining: {est_remain:.1f}s")
        # Wait for all to finish
        concurrent.futures.wait(futures)
    print(f"✅ Pre-scraping done in {time.time()-start_scrape:.1f} seconds.\n")

    print(f"▶️ Running {iterations} iterations on {num_workers} threads...")
    start_iter = time.time()
    batch_size = max(1, iterations // 20)  # print progress every 5%

    seeds = [42 + i for i in range(iterations)]
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for i, seed in enumerate(seeds, 1):
            futures.append(executor.submit(single_iteration, seed))
            if i % batch_size == 0 or i == iterations:
                elapsed = time.time() - start_iter
                percent = i / iterations * 100
                est_total = elapsed / (i / iterations)
                est_remain = est_total - elapsed
                print(f"  Completed {i}/{iterations} iterations ({percent:.1f}%), "
                      f"Elapsed: {elapsed:.1f}s, Est. remaining: {est_remain:.1f}s")
        concurrent.futures.wait(futures)
        results = [f.result() for f in futures]

    count_extreme = sum(1 for x in results if x >= observed_value)
    p_value = count_extreme / iterations
    mean_size = sum(results) / iterations

    print(f"\n✅ Empirical p-value: {p_value:.4f}")
    print(f"📊 Mean largest network size: {mean_size:.2f}")
    print(f"⏱️ Total runtime: {time.time() - start_scrape:.1f} seconds")

if __name__ == "__main__":
    main()
