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
- Pregnancy raises susceptibility: `ss.modify(malaria.rel_sus, multiply=rr_preg, where=pregnancy.pregnant)`. This is a reduced-form stand-in for higher peripheral parasitaemia in pregnancy (attractiveness to *Anopheles*, and seeding of the periphery from the placenta), not a mechanism; it does not raise onward transmission, and must not double-count the placental niche
- **Every rate carries a unit** (`ss.perday`, `ss.days`). A bare number takes the sim's unit; with a monthly sim this made transmission ~30× too low and malaria died out in the prototype.

Human infectiousness depends on peripheral infection only. Placental infection never feeds back into transmission, so no intervention can act through that route.

### `placental` (connector)

States per woman:

- `ab` (float, 0–1): VAR2CSA antibody level (latent; never observed directly)
- `parasites` (bool): placental parasites present
- `pigment` (bool): haemozoin in placental tissue; once set, never reverts within a pregnancy
- `ti_parasites`: time the current placental infection began
- `dur_parasites`, `dur_parasites_pigment` (float): accumulated time with parasites, and with parasites plus pigment, this pregnancy (stage 2 harm uses these, not the histology label)

Histology class is derived, not stored: none (neither), acute (parasites, no pigment), chronic (parasites and pigment), past (pigment only). A re-infection over old pigment is therefore chronic.

Each step, among pregnant women:

1. **New infection → placental?** Newly infected pregnant women with gestation ≥ `ga_plac` (~13 weeks; no placental surface before) and no current placental parasites establish placental infection with probability `p0 / (1 + (ab/ab50)^h)`.
2. **Pigment.** Set after `t_chronic` of continuous placental parasites.
3. **Clearance.** Placental parasites clear at rate `1/dur_plac` (untreated, months), **independent of peripheral recovery**: placental infection can persist with a negative peripheral smear.
4. **Boost.** On establishing placental infection: `ab += boost · (1 − ab)`. On a new infection during pregnancy that does not establish placentally: `ab += f · boost · (1 − ab)`. `f = 0` gives the placental-only variant.
5. **Waning.** All women: `ab *= exp(−ln2 · dt / ab_halflife)`.
6. **Delivery** (Pregnancy delivery callback): record the histology class at delivery, then reset the placental states. Stage 3 adds the newborn transfer here: `ab_newborn = ab · eff_transfer · (1 − imp_transfer · parasites)`.

`clear_placental(uids)` sets `parasites = False` and leaves `pigment` (acute → none, chronic → past). This is the hook for treatment and IPTp in stage 4; placental clearance by drugs is less reliable than peripheral, so stage 4 uses its own efficacy.

**Initial antibody.** At malaria start, replay each woman's `n_pregnancies` past pregnancies with a per-pregnancy infection probability `q_past` (one constant per intensity level) and the birth interval implied by her age and gravidity. Each infected past pregnancy goes through the same Hill function with her running `ab`, then boost and wane. Calibrate `q_past` per intensity so that `ab` by age at malaria start matches `ab` by age after the burn-in.

### `outcomes` (connector, stage 2 interface only)

Harm accrues with exposure and does not reverse when infection clears.

- **Growth restriction:** per-step penalty to `FetalHealth.apply_growth_restriction` while placental parasites are present, larger with pigment (`gr_acute`, `gr_chronic`). Steepest gravidity gradient.
- **Preterm birth:** mainly from peripheral (febrile) episodes, with a smaller placental term, via `FetalHealth.apply_timing_shift`. Stage 1's SEIRS has no clinical state, so peripheral infection is the proxy until one is added.
- **Miscarriage:** from peripheral infection only. It is mostly first-trimester, before placental infection can establish (`ga_plac`), and is driven by systemic illness. Hook into `ss.Pregnancy` to be confirmed.
- **Stillbirth:** from placental and peripheral infection.
- **Maternal anaemia:** state on the mother, driven by peripheral and placental infection; flattest across gravidity.

The flatter gradients for preterm birth and anaemia emerge from their dependence on peripheral infection, rather than being imposed per outcome.

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
| `rr_preg` | 1.5–2 | Reduced-form; see `Malaria` |
| `p0`, `ab50` | 0.7, 0.3 | Placental establishment vs antibody |
| `h` | 2 (fixed) | Fixed to limit the `boost`/`ab50` ridge |
| `boost` | 0.4 (0.3–0.5) | Per infection |
| `f` | 0.15 (fixed; 0 as switch) | Boost from non-placental infection in pregnancy |
| `ab_halflife` | 8 y (6–10) | Waning between pregnancies; strongest control on the gradient |
| `ga_plac` | 13 weeks | Earliest placental establishment |
| `t_chronic` | 6 weeks | Parasites → pigment |
| `dur_plac` | 90 d | Untreated placental clearance |
| `p_loss` | 0.15 | Pregnancies not ending in birth, for initial gravidity |

**Joint calibration.** `p0`, `boost`, `ab50`, `ab_halflife` and the force of infection act on one observable (the gravidity gradient) and cannot have independent priors. A per-woman Monte Carlo of the mechanism (one infection per pregnancy, 3-year birth interval) gives G2/G1 and G3/G1 of 0.71 and 0.67 with a 3-year half-life and `boost = 0.5`; reaching the observed ≈ 0.6 and 0.4 needs a half-life of 6–10 years with `boost` 0.3–0.5. At high intensity, multiple infections per pregnancy raise the effective per-pregnancy boost and push G1 towards saturation, so the per-infection values come out of calibration, not literature.

**Identifiability.** `ab` is latent and the targets are parasitaemia and histology. Expect a ridge between `boost` and `ab50`, broken only weakly by the `1 − ab` saturation. Fix `h` and `f`; calibrate `p0`, `boost`, `ab50`, `ab_halflife` jointly and check for the ridge.

## Testing

pytest in `tests/`:

- Boost/wane arithmetic; `p_plac` decreases monotonically with `ab`
- **Null gradient:** with `boost = 0` and `f = 0`, placental infection is flat across gravidity (the gradient must come from antibodies, not from age/parity structure)
- **Gradient check:** the mechanism alone (no transmission model) reproduces the Monte Carlo gradients above
- Longer `ab_halflife` gives a steeper gradient
- **Birth interval:** at fixed gravidity, longer intervals give lower `ab` (the mechanism behind the fertility-transition result)
- Pigment never reverts within a pregnancy; `clear_placental` maps acute → none and chronic → past
- No placental infection before `ga_plac`
- Placental clearance is independent of peripheral recovery
- **Initial antibody:** `ab` by age at malaria start ≈ `ab` by age after the burn-in
- Initial `n_pregnancies` matches cumulative ASFR × (1 + `p_loss`) by age
- Unit consistency: `Malaria` epi outputs agree within tolerance between a sim with a monthly unit and one with a daily unit
- Branching: continued run equals uninterrupted; a parameter change after load takes effect
- `sim.check_dt(factor=2)`: daily vs 2-day malaria
- Speed benchmark: 10k agents, 1985 monthly demography → 10-y malaria burn-in → 5-y projection; decides 1- vs 2-day `dt`

**Stage 1 success:** the 4-class histology distribution at delivery (none / acute / chronic / past) differs by gravidity in the observed direction at high intensity, within age bands, and the gradient is smaller at low intensity. Comparisons average over enough years that conception-month composition (seasonality × gestation) does not drive them. The projection step is fast enough for hundreds of scenarios.

**Calibration targets (later stages):** *Pf*PR₂₋₁₀ per intensity; peripheral parasitaemia at first ANC by gravidity (Mozambique ANC surveillance); placental histology class by gravidity and age (e.g. Manhiça cohorts; Desai et al. 2007), using HIV-negative strata where reported. Age-specific prevalence in children will not fit until stage 3 adds age-acquired immunity.

## Limitations

- **No HIV.** HIV attenuates gravidity-dependent protection (HIV-positive multigravidae resemble primigravidae), and Mozambique and the Manhiça cohorts have high HIV prevalence. Fitting without HIV biases `boost` low; prefer HIV-negative strata.
- Antibodies block placental establishment only; they do not shorten placental infections. Candidate addition if chronic-class prevalence by gravidity will not fit.
- Placental infection does not feed back into transmission.
- No postpartum parasite release.

## Open questions

- Shape of the `Scenarios` API: how a scenario is specified
- Miscarriage/stillbirth hook in v4 `ss.Pregnancy`
- Whether to add a link from `ab` to anti-VAR2CSA IgG titre so published titres by gravidity can be a target (would pin the scale of `ab` and resolve the `boost`/`ab50` ridge)
- Whether to model postpartum parasite release
- Final malaria `dt` (benchmark)
