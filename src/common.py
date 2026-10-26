"""Portable input/output helpers for the manuscript analyses."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = ('principal_32', 'principal_41', 'support_32', 'support_41', 'limma_32', 'limma_41')


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def write(path, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    records = list(records)
    if not records:
        raise ValueError(f'No records for {path}')
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def ids(path):
    values = [r['ENSEMBL_ID'].strip() for r in rows(path)]
    if not all(values) or len(values) != len(set(values)):
        raise ValueError(f'Blank or duplicated IDs in {path}')
    return values


def inputs():
    data = ROOT / 'data'
    background = ids(data / 'background_interactions.csv.gz')
    partners = {
        r['ENSEMBL_ID']: set(filter(None, r['KNOWN_INTERACTIONS'].split(';'))) - {r['ENSEMBL_ID']}
        for r in rows(data / 'background_interactions.csv.gz')
    }
    selected = {case: ids(data / 'gene_sets' / f'{case}.csv') for case in CASES}
    bridges = {temp: set(ids(data / 'deg_bridges' / f'{temp}.csv')) for temp in (32, 41)}
    housekeeping = set(ids(data / 'housekeeping_genes.csv'))
    if any(not set(genes) <= set(background) for genes in selected.values()):
        raise ValueError('Selected genes must occur in the cached background')
    return background, partners, selected, bridges, housekeeping


def verify_inputs():
    manifest = json.loads((ROOT / 'data/manifest.json').read_text())
    for entry in manifest:
        content = (ROOT / entry['path']).read_bytes()
        if hashlib.sha256(content).hexdigest() != entry['sha256']:
            raise ValueError(f'Input checksum mismatch: {entry["path"]}')
    return manifest
