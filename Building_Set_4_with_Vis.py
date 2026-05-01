import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

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

        # Symmetrize
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

def build_graph(neighbor_map, genes_to_include):
    G = nx.Graph()
    for gene in genes_to_include:
        for neigh in neighbor_map[gene]:
            if neigh in genes_to_include:
                G.add_edge(gene, neigh)
    return G

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

    # Find the largest network
    largest_network = max(disjoint_sets, key=len)
    print(f"\nLargest network has {len(largest_network)} genes. Visualizing now...")

    # Build and draw graph for largest network
    G = build_graph(neighbor_map, largest_network)

    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, k=0.1, iterations=50)  # force-directed layout, tweak if needed
    nx.draw_networkx_nodes(G, pos, node_size=20, node_color='blue', alpha=0.6)
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    plt.title(f'Gene Interaction Network (Largest Cluster, {len(largest_network)} genes)')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('largest_gene_network.png')
    plt.show()

    print("Network visualization saved as largest_gene_network.png")

if __name__ == "__main__":
    main()
