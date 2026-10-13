"""
Plot Mozambique demographics: population, pregnancies, and pregnancy by maternal age
"""

import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
from analyzers import PregnancyByAge

pal = dict(model='#2a6fb0', model2='#8fb8de', data='k')


def load_data():
    pop = pd.read_csv('data/mozambique_pop.csv')
    pop = pop[pop.Sex == 'Both sexes'].set_index('Time').Value
    asfr = pd.read_csv('data/mozambique_asfr.csv')
    asfr['group'] = (asfr.AgeStart // 5) * 5
    asfr = asfr.groupby(['Time', 'group']).Value.mean().unstack()  # Births per 1000 women, 5-year groups
    tfr = asfr.sum(axis=1) * 5 / 1000
    return pop, asfr, tfr


def plot_demographics(df, filename='figures/demographics.png', asfr_years=(1985, 2005, 2025)):
    pop, asfr, tfr = load_data()
    df = df.set_index('timevec').iloc[:-1]  # Drop final partial year
    yrs = df.index
    labels = PregnancyByAge().labels
    k = 'pregnancybyage'
    women = df[[f'{k}.women_{lb}' for lb in labels]].values
    preg = df[[f'{k}.new_pregnancies_{lb}' for lb in labels]].values
    births = df[[f'{k}.births_{lb}' for lb in labels]].values
    pregnant = df[[f'{k}.pregnant_{lb}' for lb in labels]].values
    agecols = plt.cm.Blues(np.linspace(0.35, 1, len(labels)))

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()

    # Population
    ax = axes[0]
    ax.plot(yrs, df.n_alive / 1e6, color=pal['model'], lw=2, label='Model')
    d = pop.loc[yrs.min():yrs.max()]
    ax.scatter(d.index, d / 1e6, color=pal['data'], s=12, label='UN WPP', zorder=3)
    ax.set_title('Population (millions)')
    ax.legend(frameon=False)

    # Pregnancies and births per year
    ax = axes[1]
    ax.plot(yrs, df['pregnancy.pregnancies'] / 1e6, color=pal['model'], lw=2, label='Pregnancies')
    ax.plot(yrs, df['pregnancy.births'] / 1e6, color=pal['model2'], lw=2, label='Births')
    ax.set_title('Pregnancies and births per year (millions)')
    ax.legend(frameon=False)

    # TFR
    ax = axes[2]
    ax.plot(yrs, df['pregnancy.tfr'], color=pal['model'], lw=2, label='Model')
    d = tfr.loc[yrs.min():yrs.max()]
    ax.scatter(d.index, d, color=pal['data'], s=20, label='UN WPP', zorder=3)
    ax.set_title('Total fertility rate')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)

    # Age-specific fertility vs WPP, selected years
    ax = axes[3]
    x = np.arange(len(labels))
    ycols = plt.cm.Greys(np.linspace(0.4, 1, len(asfr_years)))
    for yr, c in zip(asfr_years, ycols):
        i = yrs.get_loc(yr)
        ax.plot(x, births[i] / women[i] * 1000, color=c, lw=2, label=f'{yr} model')
        ax.scatter(x, asfr.loc[yr, [int(lb[:2]) for lb in labels]], color=c, marker='o', s=25, facecolor='none', label=f'{yr} WPP')
    ax.set_xticks(x, labels)
    ax.set_xlabel('Maternal age')
    ax.set_title('Births per 1000 women')
    ax.legend(frameon=False, ncol=2, fontsize=8)

    # Pregnancies by maternal age over time
    ax = axes[4]
    ax.stackplot(yrs, preg.T / 1e6, colors=agecols, labels=labels, edgecolor='white', lw=0.5)
    ax.set_title('Pregnancies per year by maternal age (millions)')
    ax.legend(frameon=False, fontsize=8, loc='upper left', reverse=True)

    # Share currently pregnant by age
    ax = axes[5]
    for j, lb in enumerate(labels):
        ax.plot(yrs, pregnant[:, j] / women[:, j] * 100, color=agecols[j], lw=2, label=lb)
    ax.set_title('Women currently pregnant (%)')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    sc.figlayout()
    sc.savefig(filename)
    return fig


if __name__ == '__main__':
    df = sc.loadobj('results/sim.df')
    plot_demographics(df)
