"""
Create Mozambique demographic model
"""

import sciris as sc
import stisim as sti
from analyzers import PregnancyByAge


def make_sim(seed=1, start=1985, stop=2030, n_agents=10e3, verbose=1/12):
    sim = sti.Sim(
        n_agents=n_agents,
        start=start,
        stop=stop,
        datafolder='data/',
        demographics='mozambique',
        rand_seed=seed,
        analyzers=PregnancyByAge(),
        verbose=verbose,
    )
    return sim


if __name__ == '__main__':
    from plot_demographics import plot_demographics, plot_parity
    sim = make_sim()
    sim.run()
    df = sim.to_df(resample='year', use_years=True, sep='.')
    sc.saveobj('results/sim.df', df)
    plot_demographics(df)
    plot_parity(df)
