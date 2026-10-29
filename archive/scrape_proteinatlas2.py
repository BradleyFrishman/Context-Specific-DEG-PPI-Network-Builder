import requests

def scrape_interacting_genes_simple(ensembl_id):
    """
    Fetch the ProteinAtlas interaction page for a gene (given by ensembl_id),
    extract all JavaScript function() blocks, 
    then find and return the list of interacting gene IDs (without the 'int' prefix).
    """
    # Step 1: Download the page HTML
    url = f"https://www.proteinatlas.org/{ensembl_id}/structure+interaction"
    response = requests.get(url)
    if response.status_code != 200:
        print("Failed to fetch the page.")
        return []

    html = response.text

    # Step 2: Extract all text inside {} after every 'function()' found
    functions_content = []
    start_pos = 0
    while True:
        func_index = html.find("function()", start_pos)
        if func_index == -1:  # No more function() found
            break
        brace_start = html.find("{", func_index)
        if brace_start == -1:
            break

        # Find the matching closing brace for the opening brace
        open_braces = 1
        i = brace_start + 1
        while i < len(html) and open_braces > 0:
            if html[i] == '{':
                open_braces += 1
            elif html[i] == '}':
                open_braces -= 1
            i += 1
        brace_end = i

        # Extract function content inside braces (excluding braces)
        content = html[brace_start + 1: brace_end - 1]
        functions_content.append(content)

        start_pos = brace_end  # Move on for next function()

    # Step 3: Combine all extracted contents into one big string
    combined_str = " ".join(functions_content)

    # Step 4: Search for all occurrences of '"id":"int' and extract gene IDs after that
    gene_ids = []
    search_str = '"id":"int'
    pos = 0
    while True:
        pos = combined_str.find(search_str, pos)
        if pos == -1:
            break
        start = pos + len(search_str)

        # Extract gene id characters until next '"' or ','
        end_pos = start
        while end_pos < len(combined_str) and combined_str[end_pos] not in ['"', ',']:
            end_pos += 1
        gene_id = combined_str[start:end_pos]  # This excludes the 'int' prefix
        gene_ids.append(gene_id)

        pos = end_pos  # Continue searching

    return gene_ids

# ----------- Main script part ------------
if __name__ == "__main__":
    # Ask user to enter the Ensembl gene ID
    gene_to_search = input("Enter Ensembl gene ID (e.g., ENSG00000163344): ").strip()

    genes = scrape_interacting_genes_simple(gene_to_search)

    print(f"\nFound {len(genes)} interacting genes (including duplicates):\n")
    for gene in genes:
        print(gene)
