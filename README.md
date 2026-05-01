# Manuscript Title: " A Context-Specific, Literature-Supported Framework for Validating Stress Response Models in Mammals"
# Context-Specific-DEG-PPI-Network-Builder

This repository contains the Python source code used to construct and analyze context-specific protein-protein interaction (PPI) networks for validating computational models of differential gene expression (DEG) under temperature stress.

The goal of this project is to determine whether model-selected genes form biologically meaningful and literature-supported mechanisms rather than reflecting statistical noise.

This work was developed as part of undergraduate research in computational biology and temperature-response genomics.

---

## Biological Question

When mammals experience environmental change (such as extreme heat or cold), computational models attempt to predict which genes drive the physiological response.

However, identifying differentially expressed genes (DEGs) is not enough—we must determine whether those genes actually participate in a coherent biological mechanism & exclude noisy genes.

This project tests whether a model-selected subset of genes:

- form significantly interconnected protein-protein interaction networks
        - as measured by their largest connected subnetworks (Largest Connected Components, LCCs)
- outperform degree-matched random gene sets
- show stronger biological specificity than traditional enrichment tools like STRING-db

---

## Core Method

### First-order connection

A direct protein-protein interaction between two model-selected genes.

### Second-order connection

Two model-selected genes are connected if they both interact with the same DEG.

Importantly:

**Second-order interactions are only allowed when the bridge/intermediate gene is itself also a DEG.**

This DEG-restricted approach improves biological specificity and prevents inflated connectivity caused by generic hub proteins.

---

## Primary Validation Metric

### Largest Connected Component (LCC)

The LCC is the largest subset of genes that are mutually connected through known protein-protein interactions.

We use LCC instead of raw edge count because:

- edge count may not not guarantee mechanistic coherence
- LCC may better reflect coordinated biological function
- stress-response pathways are expected to behave as cohesive functional subnetworks

---

## Null Model

### Degree-Preserving Randomization

To control for bias from highly studied genes, random gene sets are generated using degree-matched sampling.

Genes are binned based on known interaction counts (20 bins), and random sets are sampled to preserve the same degree distribution (or as similar distribution as possible) as the observed gene set.

This prevents inflated significance caused by genes with naturally high connectivity.

---

## Main Scripts

### `Building_Set_4.py`

Constructs context-specific first- and second-order DEG-restricted PPI networks.

### `Building_Set_4_with_Vis.py`

Generates network visualizations and outputs the largest connected component.

### `Build 2nd order PA sets from Background 4.py`

Builds second-order networks using Human Protein Atlas interactions and analyzes network fragmentation.

### `Full Pipeline Single File 2.py`

Runs network analysis for a single input file and reports network statistics.

### `Generate_random_csv.py`

Generates randomized degree-matched null gene sets.

### `iterations_1000.py`

Performs large-scale randomization iterations for null distribution generation.

### `Refined_Known_Interactions.py`

Processes known interactions from Protein Atlas and prepares DEG-restricted interaction mappings.

### `Connections_2.py`

Builds second-order connection columns for gene interaction datasets.

---

## Data Sources

### Human Protein Atlas

Used to identify literature-supported known protein-protein interactions. This website compiles protein-protein interactions from BioGRID, IntAct, BioPlex, and OpenCell. 

https://www.proteinatlas.org/

### STRING-db

Used for comparison with standard enrichment-based PPI validation.

https://string-db.org/

### HRT Atlas

Used for housekeeping gene validation.

---

## Key Findings

Across both 32°C and 41°C temperature stress conditions:

- >75% of Principal Response genes from the Gonzalez et al., 2026 paper (https://doi.org/10.1186/s12864-026-12823-7) (our test-case) joined the Largest Connected Component
- observed LCCs were significantly larger than degree-matched random expectations
- Support Group genes showed strong housekeeping gene enrichment
- STRING confirmed enrichment, but produced less stable and less context-specific results. (randomly produced gene sets tended to have a significant PPI-enrichment P-value)

These findings support the biological validity of the tested stress-response model.

---

## Citation

If using this code, please cite:

Frishman, B., Gonzalez, J., & Forbes, V.  
*A Context-Specific, Literature-Supported Framework for Validating Stress Response Models in Mammals*  
(Working manuscript)

---

## Author

Bradley Frishman  
University of Florida  
Undergraduate Researcher — Computational Biology / Stress Genomics

---
