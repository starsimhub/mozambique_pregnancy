# Malaria in pregnancy: design

## Goal

Test the thesis that protection against malaria in pregnancy rises with gravidity because women acquire antibodies against placental (VAR2CSA-expressing) parasites with each infected pregnancy. The model should let the gravidity effect *emerge* from exposure rather than be imposed, so that it can vary with transmission intensity and with the fertility transition.

Uses: (1) adverse-outcome burden by gravidity, (2) interventions in pregnancy, (3) infant protection via transferred antibodies, (4) effects of falling fertility on the share of first pregnancies.

## Stages

Each stage is its own spec → plan → implementation cycle. This spec details stage 1 and fixes the interfaces that stages 2–4 hook into.

1. **Core**: maternal malaria infection, placental infection, antibodies, burn-in/branch workflow
2. **Outcomes**: per-outcome risk from placental history (growth restriction, preterm, anaemia, stillbirth/miscarriage)
3. **Infants**: antibody transfer at delivery, decay, infant susceptibility; age-acquired immunity in children
4. **Interventions**: IPTp, ANC nets, screen-and-treat
5. **Fertility-transition scenarios**: analysis on 1–4

## Decisions

| Decision | Choice | Reason |
|---|---|---|
| Framework | starsim v4, `ss.Sim(solver='abm')`; drop `sti.Sim` | `ss.Pregnancy` and `FetalHealth` are agent-only in v4 (`E804`); stisim adds nothing for a non-STI and is untested on v4 |
| Vectors | Ross-Macdonald mosquito states as counts in the same `ss.Disease` | Native in v4 under `abm` (`graph.count_states`); keeps community-wide vector control possible |
| Immunity | Continuous per-woman antibody level | Gravidity effect emerges and weakens at low transmission; gives a natural quantity to transfer to the newborn |
| Setting | National demography, transmission set by one `intensity` par, run at low/moderate/high | Tests the transmission dependence of the gravidity effect without spatial structure |
| Timesteps | `Pregnancy`/`Deaths` monthly from 1985; `Malaria` daily or 2-day from its own later start | v4 interleaves modules with their own `dt` and `start` by absolute time |
| Gravidity | `Pregnancy.n_pregnancies`, not parity | Antibodies respond to pregnancies, including those not ending in birth |
| Naming | Terse: `Malaria` (module), `placental`, `outcomes` (connectors, lowercase as in `ssl.mnch.fetal_infection`), `Scenarios` (container) | Starsim conventions |

## Components

### Demography

`ss.Pregnancy(dt=ss.months(1))` and `ss.Deaths(dt=ss.months(1))` from the committed WPP files in `data/`, with the people's age distribution from `mozambique_age_1985.csv`. The existing demography figures must reproduce.

**Initial gravidity.** `n_pregnancies` starts at 0 for everyone, which makes older women look naive. At t = 0, draw each woman's `n_pregnancies` from her age: Poisson with mean equal to cumulative ASFR up to her age (1985 rates), multiplied by `1 + p_loss` (pregnancies not ending in birth). Lives in this repo; candidate to promote into `ss.Pregnancy`.

### `Malaria` (`ss.Disease`)

Humans S → E → I → R → S; mosquitoes `null` → Sv → Ev → Iv → `null` as counts.

- `infection: Sh >> Eh @ a·b·Iv.n/Nh`, `latency: Eh >> Ih @ ss.days(12)`, `recovery: Ih >> Rh @ ss.days(dur_inf)`, `waning: Rh >> Sh @ ss.days(dur_imm)`
- `emergence: null >> Sv @ intensity · base_emergence · seasonal(t)`, `biting: Sv >> Ev @ a·c·Ih.n/Nh`, `incubation: Ev >> Iv @ ss.days(10)`, `mortality: (Sv|Ev|Iv) >> null @ ss.days(10)`
- `seasonal`: `ss.timeseries(mode='cosine', period=ss.years(1))`, peak in February–March (rainy season November–April)
- Pregnancy raises susceptibility: `ss.modify(malaria.rel_sus, multiply=rr_preg, where=pregnancy.pregnant)`
- **Every rate carries a unit** (`ss.perday`, `ss.days`). A bare number takes the sim's unit; with a monthly sim this made transmission ~30× too low and malaria died out in the prototype.

Human infectiousness depends on peripheral infection only.

### `placental` (connector)

States per woman:

- `ab` (float, 0–1): VAR2CSA antibody level
- `plac` (int): 0 none, 1 acute (parasites, no pigment), 2 chronic (parasites + pigment), 3 past (pigment only)
- `ti_plac`: time placental infection began
- `dur_plac_inf`, `dur_plac_chronic` (float): accumulated time in acute+chronic and in chronic during this pregnancy (stage 2 uses these)

Each step, among pregnant women:

1. **New infection → placental?** Newly infected pregnant women with `plac ∈ {0, 3}` (none or past) become acute with probability `p0 / (1 + (ab/ab50)^h)`.
2. **Progression.** Acute → chronic after `t_chronic`. Acute/chronic → past at rate `1/dur_plac` (untreated clearance, months), **independent of peripheral recovery**: placental infection can persist with a negative peripheral smear.
3. **Boost.** On becoming placentally infected: `ab += boost · (1 − ab)`. On a new infection during pregnancy that does not become placental: `ab += f · boost · (1 − ab)`. `f = 0` gives the placental-only variant.
4. **Waning.** All women: `ab *= exp(−ln2 · dt / ab_halflife)`.
5. **Delivery** (Pregnancy delivery callback): record `plac` as the histology class at delivery, then reset `plac`, `ti_plac`, and durations. Stage 3 adds the newborn transfer here: `ab_newborn = ab · eff_transfer · (1 − imp_transfer · [plac ∈ {1, 2}])`.

**Initial antibody.** At malaria start, set each woman's `ab` by applying the boost/wane rules over her `n_pregnancies` past pregnancies, with the per-pregnancy probability of placental infection taken from a scenario-level par `p_plac_past` (and peripheral-only infection `p_inf_past`). These are calibrated per intensity level so that `ab` by age is stationary across the malaria burn-in.

Postpartum parasite release at delivery is out of scope for stage 1.

### `outcomes` (connector, stage 2 interface only)

Harm accrues with exposure and does not reverse when infection clears. For outcome *k*:

- Growth restriction: per-step penalty to `FetalHealth.apply_growth_restriction` while `plac ∈ {1, 2}`, larger in chronic (`gr_acute`, `gr_chronic`)
- Preterm: `FetalHealth.apply_timing_shift` with probability depending on placental history
- Stillbirth/miscarriage: via `ss.Pregnancy` (hook to be confirmed)
- Maternal anaemia: state on the mother, driven by peripheral and placental infection, flatter across gravidity

Shapes differ by outcome (growth restriction steepest first-pregnancy effect, preterm intermediate, anaemia near-flat), which comes from each outcome having its own dependence on placental history rather than one curve.

### `Scenarios` (container, prototype)

Burn in once, branch many times.

- **Burn-ins.** For each intensity level, K sims with different `rand_seed`, run with `sim.run(until=branch_date)` and saved: `burnins/<intensity>_<k>.sim`. A sim keeps one seed for its whole life; there is no reseeding at the branch point. Stochastic spread comes from the K seeds.
- **Branching.** `Scenarios(burnins, scenarios)` loads each burn-in, applies each scenario's changes, and runs to `stop`, in parallel. Scenarios compared within a burn-in share its random streams.
- **Scenario changes.** Parameter changes after load are confirmed to work. Interventions are declared (inactive) in the burn-in and activated at the branch point; adding modules after load is not relied on. How a scenario is specified (`ss.modify`, an intervention, a par dict) is open.
- `dt` cannot change at the branch point (v4 `Disease.dt` has no setter), so malaria runs at its final `dt` through its burn-in; the cost is paid once per burn-in.

Verified in v4 prototypes: continuation after save/load is bit-identical to an uninterrupted run; a saved 5,000-agent sim is 258 KB and loads in < 1 s.

## Parameters (priors, to be calibrated)

| Par | Prior | Notes |
|---|---|---|
| `a`, `b`, `c` | 0.3/day, 0.5, 0.5 | Ross-Macdonald |
| `dur_inf`, `dur_imm` | 100 d, 365 d | Untreated infection, temporary immunity |
| `intensity` | low / moderate / high | Calibrated to MAP *Pf*PR₂₋₁₀ ranges for Mozambique |
| `rr_preg` | 1.5–2 | Pregnant women more attractive to *Anopheles* |
| `p0`, `ab50`, `h` | 0.7, 0.3, 2 | Placental establishment vs antibody |
| `boost`, `f` | 0.5, 0.15 (0–0.2) | Graded boost; `f = 0` switch |
| `ab_halflife` | 3 y | Waning between pregnancies |
| `t_chronic` | 6 weeks | Acute → chronic (pigment) |
| `dur_plac` | 90 d | Untreated placental clearance |
| `p_loss` | 0.15 | Pregnancies not ending in birth, for initial gravidity |

## Testing

pytest in `tests/`:

- Boost/wane arithmetic; `p_plac` decreases monotonically with `ab`
- Unit consistency: `Malaria` gives the same results in a sim with a monthly unit and one with a daily unit
- Initial `n_pregnancies` matches cumulative ASFR × (1 + `p_loss`) by age
- Placental clearance is independent of peripheral recovery
- Branching: continued run equals uninterrupted; a parameter change after load takes effect
- `sim.check_dt(factor=2)`: daily vs 2-day malaria
- Speed benchmark: 10k agents, 1985 monthly demography → 10-y malaria burn-in → 5-y projection; decides 1- vs 2-day `dt`

**Stage 1 success:** placental infection prevalence at delivery is clearly higher in first than later pregnancies at high intensity; the gradient is smaller at low intensity; the projection step is fast enough for hundreds of scenarios.

**Calibration targets (later stages):** *Pf*PR₂₋₁₀ per intensity; peripheral parasitaemia at first ANC by gravidity (Mozambique ANC surveillance); placental histology class by gravidity (e.g. Manhiça cohorts; Desai et al. 2007). Age-specific prevalence in children will not fit until stage 3 adds age-acquired immunity.

## Open questions

- Shape of the `Scenarios` API: how a scenario is specified
- Stillbirth/miscarriage hook in v4 `ss.Pregnancy`
- Whether to model postpartum parasite release
- Final malaria `dt` (benchmark)
