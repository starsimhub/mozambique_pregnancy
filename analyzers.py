"""
Analyzers for pregnancy by maternal age
"""

import numpy as np
import starsim as ss


class PregnancyByAge(ss.Analyzer):
    """
    Women, currently pregnant, new pregnancies, births, and parity by 5-year
    maternal age group. Parity is 0 for everyone at sim start, so values are only
    complete for cohorts aged <= min_age at the start.
    """

    def __init__(self, age_bins=np.arange(15, 55, 5), max_parity=6, **kwargs):
        super().__init__(**kwargs)
        self.age_bins = age_bins
        self.labels = [f'{a}-{b-1}' for a, b in zip(age_bins[:-1], age_bins[1:])]
        self.parities = [str(p) for p in range(max_parity)] + [f'{max_parity}+']
        return

    def init_results(self):
        super().init_results()
        results = []
        for lb in self.labels:
            results += [
                ss.Result(f'women_{lb}', summarize_by='mean', scale=True),
                ss.Result(f'pregnant_{lb}', summarize_by='mean', scale=True),
                ss.Result(f'new_pregnancies_{lb}', summarize_by='sum', scale=True),
                ss.Result(f'births_{lb}', summarize_by='sum', scale=True),
                ss.Result(f'mean_parity_{lb}', summarize_by='mean', scale=False),
            ]
            results += [ss.Result(f'parity{p}_{lb}', summarize_by='mean', scale=True) for p in self.parities]
        self.define_results(*results)
        return

    def step(self):
        ppl = self.sim.people
        preg = self.sim.demographics.pregnancy
        ti = self.ti
        f = ppl.female & (ppl.age >= 0)
        new_preg = f & (preg.ti_pregnant == ti)
        births = f & ~preg.pregnant & (preg.ti_delivery == ti)
        age = ppl.age
        top = len(self.parities) - 1
        for lb, a, b in zip(self.labels, self.age_bins[:-1], self.age_bins[1:]):
            in_age = (age >= a) & (age < b)
            self.results[f'women_{lb}'][ti] = (f & in_age).count()
            self.results[f'pregnant_{lb}'][ti] = (f & preg.pregnant & in_age).count()
            self.results[f'new_pregnancies_{lb}'][ti] = (new_preg & in_age).count()
            self.results[f'births_{lb}'][ti] = (births & in_age).count()
            par = preg.parity[(f & in_age).uids]
            self.results[f'mean_parity_{lb}'][ti] = par.mean() if len(par) else np.nan
            counts = np.bincount(np.minimum(par, top).astype(int), minlength=top + 1)
            for p, c in zip(self.parities, counts):
                self.results[f'parity{p}_{lb}'][ti] = c
        return
