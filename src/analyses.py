"""Regenerate the reported network, housekeeping, and degree results."""
import random
from collections import Counter
import numpy as np
from scipy.stats import fisher_exact
from .common import CASES, ROOT, inputs, rows, write
from .network import networks, sampling_setup

SEED = 20261011
N_MAIN = 1000
CHECKPOINTS = (100, 500, 1000, 1500)


def network_analysis(out):
    background, partners, selected, bridges, _ = inputs()
    universe = sorted(set().union(*partners.values()))
    index = {gene: i for i, gene in enumerate(universe)}
    full_masks = {g: sum(1 << index[p] for p in ps) for g, ps in partners.items()}
    summaries, paired, sensitivity, membership, edges = [], [], [], [], []
    for case in CASES:
        genes = [g for g in selected[case] if partners[g]]
        temp = int(case.rsplit('_', 1)[1])
        degs = bridges[temp]
        masks = {g: sum(1 << index[p] for p in ps & degs) for g, ps in partners.items()}
        setup, _, _ = sampling_setup(genes, background, partners)
        actual_r = networks(genes, partners, masks)
        actual_u = networks(genes, partners, full_masks)
        for gene in genes:
            membership.append(dict(analysis=case, ENSEMBL_ID=gene, degree=len(partners[gene])))
        for i, a in enumerate(genes):
            for b in genes[i + 1:]:
                first = b in partners[a] or a in partners[b]
                second = bool(masks[a] & masks[b])
                if first or second:
                    edges.append(dict(analysis=case, gene_a=a, gene_b=b,
                                      first_order=int(first), second_order=int(second)))
        restricted = {mode: [] for mode in actual_r}
        unrestricted = {mode: [] for mode in actual_r}
        for i in range(max(CHECKPOINTS)):
            rng = random.Random(f'{SEED}:{case}:{i}')
            sample = [g for pool, count in setup for g in rng.sample(pool, count)]
            result_r = networks(sample, partners, masks)
            result_u = networks(sample, partners, full_masks) if i < N_MAIN else None
            for mode in actual_r:
                restricted[mode].append(result_r[mode])
                if result_u is not None:
                    unrestricted[mode].append(result_u[mode])
                    for j, metric in enumerate(('lcc_size', 'lcc_edges')):
                        paired.append(dict(analysis=case, iteration=i + 1, network=mode, metric=metric,
                                           restricted=result_r[mode][j], unrestricted=result_u[mode][j],
                                           loss=result_u[mode][j] - result_r[mode][j]))
        for mode in actual_r:
            for j, metric in enumerate(('lcc_size', 'lcc_edges')):
                r_all = np.array(restricted[mode])[:, j]
                r = r_all[:N_MAIN]
                u = np.array(unrestricted[mode])[:, j]
                loss = u - r
                vr, vu = actual_r[mode][j], actual_u[mode][j]
                vl = vu - vr
                summaries.append(dict(
                    analysis=case, network=mode, metric=metric, n_genes=len(genes),
                    observed_restricted=vr, observed_unrestricted=vu, observed_loss=vl,
                    restricted_null_mean=float(r.mean()), unrestricted_null_mean=float(u.mean()),
                    null_mean_loss=float(loss.mean()),
                    p_restricted=(1 + int(np.sum(r >= vr))) / (N_MAIN + 1),
                    p_restricted_lower=(1 + int(np.sum(r <= vr))) / (N_MAIN + 1),
                    p_unrestricted=(1 + int(np.sum(u >= vu))) / (N_MAIN + 1),
                    p_smaller_loss=(1 + int(np.sum(loss <= vl))) / (N_MAIN + 1),
                    p_larger_loss=(1 + int(np.sum(loss >= vl))) / (N_MAIN + 1)))
                for tail in ('upper', 'lower'):
                    record = dict(gene_set=case, network=mode, metric=metric, tail=tail, observed=vr)
                    for n in CHECKPOINTS:
                        count = int(np.sum(r_all[:n] >= vr if tail == 'upper' else r_all[:n] <= vr))
                        record[f'p_{n}'] = (count + 1) / (n + 1)
                    record['same_classification_all_n'] = len({record[f'p_{n}'] < .05 for n in CHECKPOINTS}) == 1
                    sensitivity.append(record)
        print(f'{case}: regenerated paired networks and iteration sensitivity', flush=True)
    write(out / 'network_summary.csv', summaries)
    write(out / 'paired_randomizations.csv', paired)
    write(out / 'iteration_sensitivity.csv', sensitivity)
    write(out / 'observed_network_nodes.csv', membership)
    write(out / 'observed_network_edges.csv', edges)


def housekeeping_and_degree(out):
    background, partners, selected, bridges, hk = inputs()
    bg = set(background)
    degree = {g: len(ps) for g, ps in partners.items()}
    hk_summary, degree_summary, random_medians, gene_degrees, descriptive = [], [], [], [], []
    for temp in (32, 41):
        g3 = set(selected[f'support_{temp}'])
        principal = set(selected[f'principal_{temp}'])
        assert not g3 & principal
        h = len(g3 & hk)
        for comparison, other in [('other_model_groups', principal), ('background_outside_g3', bg - g3)]:
            k = len(other & hk)
            table = [[h, len(g3) - h], [k, len(other) - k]]
            hk_summary.append(dict(temperature=temp, comparison=comparison,
                                   g3_hk=h, g3_non_hk=len(g3) - h,
                                   comparison_hk=k, comparison_non_hk=len(other) - k,
                                   g3_hk_percent=100 * h / len(g3),
                                   comparison_hk_percent=100 * k / len(other),
                                   full_background_hk=len(bg & hk), full_background_n=len(bg),
                                   full_background_hk_percent=100 * len(bg & hk) / len(bg),
                                   fold_enrichment=(h / len(g3)) / (len(bg & hk) / len(bg)),
                                   g3_share_selected_percent=100 * len(g3) / len(g3 | principal),
                                   g3_share_selected_hk_percent=100 * h / len((g3 | principal) & hk),
                                   fisher_p=float(fisher_exact(table, alternative='greater').pvalue)))
        for label in ('All genes', 'Non-housekeeping only'):
            genes = g3 if label == 'All genes' else g3 - hk
            universe = bg if label == 'All genes' else bg - hk
            pool = np.array([degree[g] for g in sorted(universe)], dtype=float)
            observed = float(np.median([degree[g] for g in genes]))
            rng = np.random.default_rng(20261021 + temp + (1000 if label != 'All genes' else 0))
            sims = np.array([np.median(rng.choice(pool, len(genes), replace=False)) for _ in range(N_MAIN)])
            count = int(np.sum(sims >= observed))
            degree_summary.append(dict(temperature=temp, analysis=label, n_genes=len(genes),
                                       n_background=len(pool), observed_median=observed,
                                       background_median=float(np.median(pool)), tail_count=count,
                                       iterations=N_MAIN, empirical_p=(count + 1) / (N_MAIN + 1)))
            random_medians.extend(dict(temperature=temp, analysis=label, iteration=i + 1, median_degree=v)
                                  for i, v in enumerate(sims))
        for name, genes in [('background', bg), ('all_degs', bridges[temp]),
                            ('principal', principal), ('support', g3)]:
            values = np.array([degree[g] for g in sorted(genes)], dtype=float)
            descriptive.append(dict(temperature=temp, gene_set=name, n_genes=len(values),
                                    median_degree=float(np.median(values)), mean_degree=float(values.mean()),
                                    sample_sd=float(values.std(ddof=1))))
        for g in sorted(g3):
            gene_degrees.append(dict(temperature=temp, gene=g, degree=degree[g], housekeeping=int(g in hk)))
    write(out / 'housekeeping_fisher.csv', hk_summary)
    write(out / 'g3_degree_summary.csv', degree_summary)
    write(out / 'g3_random_medians.csv', random_medians)
    write(out / 'g3_gene_degrees.csv', gene_degrees)
    write(out / 'degree_descriptive_statistics.csv', descriptive)
    print('Regenerated housekeeping tests and degree analyses', flush=True)


def verify_results(out):
    """Compare every numerical field with frozen reference files."""
    checks = [('network_summary.csv', 'network_summary.csv'),
              ('paired_randomizations.csv', 'paired_randomizations.csv.gz'),
              ('iteration_sensitivity.csv', 'iteration_sensitivity.csv'),
              ('g3_random_medians.csv', 'g3_random_medians.csv.gz'),
              ('g3_degree_summary.csv', 'g3_degree_summary.csv'),
              ('housekeeping_fisher.csv', 'housekeeping_fisher.csv'),
              ('degree_descriptive_statistics.csv', 'degree_descriptive_statistics.csv')]
    count = 0
    for actual_name, expected_name in checks:
        actual, expected = rows(out / actual_name), rows(ROOT / 'reference' / expected_name)
        if len(actual) != len(expected):
            raise AssertionError(f'{actual_name}: row count mismatch')
        for i, (a, e) in enumerate(zip(actual, expected), 2):
            for key in e:
                try:
                    equal = np.isclose(float(a[key]), float(e[key]), rtol=1e-12, atol=1e-12)
                except ValueError:
                    equal = a[key] == e[key]
                if not equal:
                    raise AssertionError(f'{actual_name}, row {i}, {key}: {a[key]} != {e[key]}')
                count += 1
    print(f'PASS: {count:,} reference values match', flush=True)
    return count
