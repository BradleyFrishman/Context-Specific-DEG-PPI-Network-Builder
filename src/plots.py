"""Regenerate presentation plots from the same numerical outputs and inputs."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from .common import ROOT, inputs, rows


def make_plots(out):
    directory = out / 'figures'
    directory.mkdir(exist_ok=True)
    background, partners, selected, bridges, _ = inputs()
    colors = {'Background': '#808080', 'All DEGs': '#e69f00',
              'Principal Response': '#0072b2', 'G3': '#cc79a7'}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    grid = np.linspace(0, 200, 600)
    for row, temp in enumerate((32, 41)):
        groups = {'Background': background, 'All DEGs': bridges[temp],
                  'Principal Response': selected[f'principal_{temp}'], 'G3': selected[f'support_{temp}']}
        for name, genes in groups.items():
            values = np.array([len(partners[g]) for g in genes])
            density = gaussian_kde(values)(grid)
            axes[row, 0].plot(grid, density, color=colors[name], label=name)
            axes[row, 1].hist(values, bins=np.arange(0, 1020, 20), density=True,
                              histtype='step', color=colors[name], label=name)
        for col in (0, 1):
            axes[row, col].set(xlim=(0, 200), xlabel='Number of recorded interaction partners',
                               ylabel='Density', title=f'{temp}°C')
            axes[row, col].legend(fontsize=8)
    fig.savefig(directory / 'degree_distributions.png', dpi=200)
    fig.savefig(directory / 'degree_distributions.pdf')
    plt.close(fig)
    records = [r for r in rows(out / 'network_summary.csv') if r['metric'] == 'lcc_size']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, mode, title in zip(axes, ('second_only', 'first_and_second'),
                               ('Second-order only', 'First + second order')):
        subset = [r for r in records if r['network'] == mode]
        x = np.arange(len(subset))
        ax.bar(x - .19, [float(r['observed_restricted']) for r in subset], .38,
               color='#0072b2', label='Observed')
        ax.bar(x + .19, [float(r['restricted_null_mean']) for r in subset], .38,
               color='#999999', label='Null mean')
        ax.set_xticks(x, [r['analysis'].replace('principal', 'Principal').replace('support', 'G3')
                         .replace('_', '\n') + '°C' for r in subset], fontsize=8)
        ax.set(title=title, ylabel='DEG-restricted LCC size')
        ax.legend(fontsize=8)
    fig.savefig(directory / 'restricted_lcc_sizes.png', dpi=200)
    fig.savefig(directory / 'restricted_lcc_sizes.pdf')
    plt.close(fig)
    records = rows(ROOT / 'archive/degree_shift/reported_means.csv')
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, temp in zip(axes, (32, 41)):
        subset = [r for r in records if int(r['temperature']) == temp]
        values = [float(r['lcc_size']) for r in subset]
        bars = ax.bar(range(5), values, color=['#d55e00', '#d55e00', '#009e73', '#009e73', '#0072b2'])
        ax.bar_label(bars, labels=[r['lcc_size'] for r in subset], padding=3, fontsize=8)
        ax.set_xticks(range(5), [r['scenario'].replace(' ', '\n') for r in subset], fontsize=8)
        ax.set(title=f'{temp}°C — reported degree-shift values', ylabel='LCC size (random mean or observed)',
               ylim=(0, 160))
    fig.savefig(directory / 'degree_shift_reported_means.png', dpi=200)
    fig.savefig(directory / 'degree_shift_reported_means.pdf')
    plt.close(fig)
