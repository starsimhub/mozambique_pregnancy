"""
Analyzers for pregnancy by maternal age
"""

import numpy as np
import starsim as ss


class PregnancyByAge(ss.Analyzer):
    """ Women, currently pregnant, new pregnancies, and births by 5-year maternal age group """

    def __init__(self, age_bins=np.arange(15, 55, 5), **kwargs):
        super().__init__(**kwargs)
        self.age_bins = age_bins
        self.labels = [f'{a}-{b-1}' for a, b in zip(age_bins[:-1], age_bins[1:])]
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
            ]
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
        for lb, a, b in zip(self.labels, self.age_bins[:-1], self.age_bins[1:]):
            in_age = (age >= a) & (age < b)
            self.results[f'women_{lb}'][ti] = (f & in_age).count()
            self.results[f'pregnant_{lb}'][ti] = (f & preg.pregnant & in_age).count()
            self.results[f'new_pregnancies_{lb}'][ti] = (new_preg & in_age).count()
            self.results[f'births_{lb}'][ti] = (births & in_age).count()
        return
