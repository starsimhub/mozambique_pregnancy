"""
Create Mozambique demographic model
"""

import sciris as sc
import stisim as sti


def make_sim(seed=1, start=1985, stop=2030, n_agents=10e3, verbose=1/12):
    sim = sti.Sim(
        n_agents=n_agents,
        start=start,
        stop=stop,
        datafolder='data/',
        demographics='mozambique',
        rand_seed=seed,
        verbose=verbose,
    )
    return sim


if __name__ == '__main__':
    sim = make_sim()
    sim.run()
    sim.plot()
    sc.savefig('figures/demographics.png')
