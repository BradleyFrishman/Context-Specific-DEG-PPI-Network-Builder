import pandas as pd

def load_neighbors(csv_path):
    df = pd.read_csv(csv_path, dtype=str).fillna('')
    neighbor_map = {}

    for _, row in df.iterrows():
        gene = row['ENSEMBL_ID'].strip()
        raw = row['REFINED_KNOWN_INTERACTIONS'] + ';' + row['CONNECTED_GENES']
        neighs = {g.strip() for g in raw.split(';') if g.strip()}

        if gene not in neighbor_map:
            neighbor_map[gene] = set()

        neighbor_map[gene].update(neighs)

        # Symmetrize: for each neighbor, add the gene back
        for n in neighs:
            if n not in neighbor_map:
                neighbor_map[n] = set()
            neighbor_map[n].add(gene)

    return df['ENSEMBL_ID'].tolist(), neighbor_map


def expand_network(seed_gene, neighbor_map, available_genes):
    network_genes = set([seed_gene])
    current_frontier = set([seed_gene])

    round_num = 1
    while current_frontier:
        new_genes = set()
        for gene in current_frontier:
            neighbors = neighbor_map.get(gene, set())
            new_neighbors = neighbors.intersection(available_genes) - network_genes
            new_genes.update(new_neighbors)

        print(f"  Round {round_num}: {len(new_genes)} new genes")
        if not new_genes:
            break

        network_genes.update(new_genes)
        current_frontier = new_genes
        round_num += 1

    return network_genes


def main():
    csv_file = 'protein_names_processes_and_interactions_with_connections.csv'
    ordered_genes_list, neighbor_map = load_neighbors(csv_file)
    unassigned_genes = set(ordered_genes_list)

    disjoint_sets = []
    network_num = 1

    for gene in ordered_genes_list:
        if gene in unassigned_genes:
            print(f"\n=== Building network {network_num} starting from {gene} ===")
            closed_set = expand_network(gene, neighbor_map, unassigned_genes)

            print(f"  Final size of network {network_num}: {len(closed_set)} genes")
            disjoint_sets.append(closed_set)

            unassigned_genes -= closed_set
            network_num += 1

    total_genes_in_networks = sum(len(s) for s in disjoint_sets)
    print(f"\nTotal unique genes across all networks: {total_genes_in_networks}")
    print(f"Number of disjoint networks built: {len(disjoint_sets)}")

    # Optional: print executive summary
    print("\nExecutive Summary:")
    sizes = [len(s) for s in disjoint_sets]
    multi_gene_networks = [s for s in sizes if s > 1]
    print(f"  {len(disjoint_sets)} total disjoint networks")
    print(f"  {len(multi_gene_networks)} have 2 or more genes")
    for n in sorted(set(multi_gene_networks)):
        count = multi_gene_networks.count(n)
        print(f"    {count} networks have {n} genes")

if __name__ == "__main__":
    main()
