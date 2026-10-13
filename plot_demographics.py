"""
Plot Mozambique demographics: population, pregnancies, and pregnancy by maternal age
"""

import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as plt
from analyzers import PregnancyByAge, AgePyramid

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


def plot_parity(df, filename='figures/parity.png', start_plot=2000, min_age=15):
    """ Parity by age; dashed where the age group includes women older than min_age at sim start """
    df = df.set_index('timevec').iloc[:-1]  # Drop final partial year
    start = df.index.min()
    df = df.loc[start_plot:]
    yrs = df.index
    an = PregnancyByAge()
    labels, parities = an.labels, an.parities
    uppers = an.age_bins[1:]
    k = 'pregnancybyage'
    agecols = plt.cm.Blues(np.linspace(0.35, 1, len(labels)))
    parcols = plt.cm.Purples(np.linspace(0.15, 1, len(parities)))

    def plot_split(ax, y, upper, color, label):
        complete = (yrs - upper) >= (start - min_age)
        ax.plot(yrs, np.where(complete, np.nan, y), color=color, lw=2, ls='--')
        ax.plot(yrs, np.where(complete, y, np.nan), color=color, lw=2, label=label)
        if complete.any():  # Join the dashed and solid segments
            i = np.argmax(complete)
            ax.plot(yrs[max(i-1, 0):i+1], y[max(i-1, 0):i+1], color=color, lw=2, ls='--')

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    for lb, ub, c in zip(labels, uppers, agecols):
        plot_split(ax, df[f'{k}.mean_parity_{lb}'].values, ub, c, lb)
    ax.set_title('Mean parity by age')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    ax = axes[1]
    yr = yrs.max()
    dist = np.array([[df.loc[yr, f'{k}.parity{p}_{lb}'] for p in parities] for lb in labels])
    dist = dist / dist.sum(axis=1, keepdims=True) * 100
    bottom = np.zeros(len(labels))
    for j, (p, c) in enumerate(zip(parities, parcols)):
        ax.bar(labels, dist[:, j], bottom=bottom, color=c, edgecolor='white', lw=0.5, label=p)
        bottom += dist[:, j]
    ax.set_title(f'Parity distribution by age, {yr} (%)')
    ax.set_xlabel('Age')
    ax.legend(frameon=False, fontsize=8, title='Parity', title_fontsize=8, bbox_to_anchor=(1, 1), loc='upper left', reverse=True)

    ax = axes[2]
    for lb, ub, c in zip(labels, uppers, agecols):
        y = df[f'{k}.parity0_{lb}'] / df[f'{k}.women_{lb}'] * 100
        plot_split(ax, y.values, ub, c, lb)
    ax.set_title('Nulliparous women by age (%)')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    sc.figlayout()
    sc.savefig(filename)
    return fig


def plot_pyramids(df, filename='figures/age_pyramids.png', years=(1985, 2005, 2025)):
    """ Age pyramids by sex, model (bars) vs UN WPP (outline), as % of total population """
    df = df.set_index('timevec')
    an = AgePyramid()
    labels, bins = an.labels, an.age_bins
    data = pd.read_csv('data/mozambique_age_sex.csv')
    data['group'] = pd.cut(data.AgeStart, bins, right=False, labels=labels)
    data = data.groupby(['Time', 'Sex', 'group'], observed=False).Value.sum()
    y = np.arange(len(labels))
    sexes = dict(m=('Male', '#2a6fb0', -1), f=('Female', '#c2577a', 1))

    fig, axes = plt.subplots(1, len(years), figsize=(5*len(years), 5), sharey=True)
    for ax, yr in zip(axes, years):
        model = {s: df.loc[yr, [f'agepyramid.{s}_{lb}' for lb in labels]].values for s in sexes}
        mtot = sum(v.sum() for v in model.values())
        dtot = data.loc[yr].sum()
        for s, (name, color, sign) in sexes.items():
            ax.barh(y, sign * model[s] / mtot * 100, height=0.85, color=color, alpha=0.8, label=f'{name} (model)')
            d = data.loc[(yr, name)].values / dtot * 100
            ax.plot(sign * np.repeat(d, 2), np.column_stack([y - 0.5, y + 0.5]).ravel(), color='k', lw=1.2,
                    label='UN WPP' if s == 'f' else None)
        lim = ax.get_xlim()
        m = max(abs(lim[0]), abs(lim[1]))
        ax.set_xlim(-m, m)
        ticks = ax.get_xticks()
        ax.set_xticks(ticks, [f'{abs(t):g}' for t in ticks])
        ax.set_xlim(-m, m)
        ax.set_title(str(yr))
        ax.set_xlabel('% of population')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_yticks(y, labels)
    axes[0].legend(frameon=False, fontsize=8, loc='upper left')
    sc.figlayout()
    sc.savefig(filename)
    return fig


if __name__ == '__main__':
    df = sc.loadobj('results/sim.df')
    plot_demographics(df)
    plot_parity(df)
    plot_pyramids(df)
