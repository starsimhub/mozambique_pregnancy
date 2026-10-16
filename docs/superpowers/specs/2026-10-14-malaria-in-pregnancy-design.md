# Malaria in pregnancy: design

*Revised 2026-10-15 to the full project scope (concept memo, kept locally, not in the repo).*

## Goal

The model is Activity 2 of a project that couples an analysis of health-system barriers (Activity 1) with a malaria transmission and impact model. It should:

1. identify the population segments and settings where malaria in pregnancy causes the most harm;
2. take segment-specific effective coverage of integrated ANC-malaria care from Activity 1;
3. estimate the impact of targeted interventions on maternal and neonatal outcomes under realistic health-system constraints, across settings on the path to elimination.

Two dynamics have to emerge, not be imposed:

- **Gravidity:** protection against placental malaria rises with gravidity, because women acquire VAR2CSA antibodies with each infected pregnancy.
- **Immunity transition:** as transmission falls, women reach pregnancy with less general immunity, so clinical and severe malaria in pregnancy rises.

## Stages

Each stage is its own spec → plan → implementation cycle. This spec details stage 1 and fixes the interfaces that stages 2–6 hook into.

1. **Transmission and immunity:** `Malaria` with exposure-acquired immunity (`rel_sev`), `placental` (VAR2CSA antibodies, placental infection), initial gravidity and immunity, burn-in/branch workflow
2. **Outcomes:** maternal (severe malaria, severe anaemia, death) and birth/neonatal (growth restriction and low birth weight, preterm, miscarriage, stillbirth, neonatal death)
3. **Segments and ANC cascade:** segment attributes; ANC process (first contact, visits, IPTp-SP with resistance, RDT, iron-folic acid, deworming, referral); input format for Activity 1's cascade table
4. **Co-factors:** background anaemia, helminth co-infection, iron-malaria interaction
5. **Infants:** transfer of maternal immunity at delivery, impaired by placental infection
6. **Scenarios:** settings on the path to elimination, intervention packages, segment-level impact

## Decisions

| Decision | Choice | Reason |
|---|---|---|
| Framework | starsim v4, `ss.Sim(solver='abm')`; drop `sti.Sim` | `ss.Pregnancy` and `FetalHealth` are agent-only (`E804`); per-woman gravidity histories, segment-specific ANC pathways and dose counts need agents; stisim adds nothing for a non-STI |
| Vectors | Ross-Macdonald mosquito states as counts in the same `ss.Disease` | Native in v4 under `abm` (`graph.count_states`); needed for transmission-decline scenarios and vector control |
| General immunity | Per-agent `rel_sev`, lowered by each infection and waning | Malaria immunity mainly prevents disease, not infection: keeps the asymptomatic reservoir, separates general from placental immunity, and makes the immunity transition show up as more clinical disease |
| Placental immunity | Per-woman VAR2CSA antibody level `placental.ab` | Gravidity effect emerges and weakens at low transmission |
| Settings | Geography-agnostic bundles (below); study geographies not yet chosen | Activity 1 selects the 2–3 geographies |
| Timesteps | `Pregnancy`/`Deaths` monthly; `Malaria` daily or 2-day from its own later `start` | v4 interleaves modules with their own `dt` and `start` by absolute time |
| Gravidity | `Pregnancy.n_pregnancies`, not parity | Antibodies respond to pregnancies, including those not ending in birth |
| Naming | Terse: `Malaria` (module), `placental`, `outcomes`, `anc` (connectors/interventions, lowercase as in `ssl.mnch.fetal_infection`), `Scenarios` (container) | Starsim conventions |

### Settings

A setting is a bundle, one per study geography or scenario:

- `country`: demography source (default Mozambique, from the committed WPP files)
- `intensity` and its trajectory: a time series multiplying mosquito emergence (constant, or declining towards elimination)
- `sp_resistance`: absent or prevalent (stage 3)
- `cascade`: Activity 1's cascade table (stage 3)
- `past`: exposure history used to initialise immunity (`foi_past`, `q_past`), calibrated per setting

## Components

### Demography

`ss.Pregnancy(dt=ss.months(1))` and `ss.Deaths(dt=ss.months(1))` from the WPP files in `data/`, with the people's age distribution from `mozambique_age_1985.csv`. The existing demography figures must reproduce.

**Initial gravidity.** `n_pregnancies` starts at 0 for everyone, which makes older women look naive. At t = 0, draw each woman's `n_pregnancies` from her age: Poisson with mean equal to cumulative ASFR up to her age (1985 rates), multiplied by `1 + p_loss` (pregnancies not ending in birth). Lives in this repo; candidate to promote into `ss.Pregnancy`.

### `Malaria` (`ss.Disease`)

**Humans:**

```
S → E → C (clinical)      p_clin = p_clin0 · rel_sev   (branch on onset)
      → A (asymptomatic)  remainder
C → S  treated            p_treat, dur_treat
C → A  untreated          dur_clin
A → S  clearance          dur_asym
```

- **`rel_sev`** (per agent, base default 1, registered as an `ss.modify` target like `rel_sus`/`rel_trans`). Its base is the immunity memory: each new infection sets `rel_sev ← rel_sev · (1 − sev_boost)`, and between infections it relaxes towards 1 with half-life `sev_halflife`. Children build immunity over repeated infections; when transmission falls, immunity wanes and clinical disease rises at older ages.
- Later effects (transferred maternal immunity in infants, HIV, pregnancy) are `ss.modify(malaria.rel_sev, ...)` modifiers on top of the base, so they never touch the immunity memory.
- Severe disease (stage 2) branches from C with `p_sev = p_sev0 · rel_sev^k`, k > 1, since protection against severe disease is acquired after fewer infections than protection against clinical disease.
- Background case management of clinical episodes (`p_treat`) is part of stage 1, because it shapes prevalence; ANC-based care is stage 3.
- `rel_sus` is reserved for factors that change the chance of infection: `rr_preg` (below) and vector control.

**Mosquitoes** (counts): `emergence: null >> Sv @ intensity(t) · base_emergence · seasonal(t)`, `biting: Sv >> Ev @ a·c·(Ch.n + r_A·Ah.n)/Nh`, `incubation: Ev >> Iv @ ss.days(10)`, `mortality: (Sv|Ev|Iv) >> null @ ss.days(10)`. `seasonal` is `ss.timeseries(mode='cosine', period=ss.years(1))`, peaking in February–March for Mozambique. Asymptomatic infections are less infectious (`r_A < 1`) but last longer, and form the reservoir.

**Infection:** `S >> E @ a·b·Iv.n/Nh`, scaled per agent by `rel_sus`.

- Pregnancy: `ss.modify(malaria.rel_sus, multiply=rr_preg, where=pregnancy.pregnant)`. A reduced-form stand-in for higher peripheral parasitaemia in pregnancy (attractiveness to *Anopheles*; seeding from the placenta), not a mechanism; it does not raise onward transmission and must not double-count the placental niche.
- **Every rate carries a unit** (`ss.perday`, `ss.days`). A bare number takes the sim's unit; with a monthly sim this made transmission ~30× too low in a prototype.

**Initial immunity.** At malaria start, set each agent's `rel_sev` by replaying an expected number of past infections from age and the setting's `foi_past`, through the same boost/wane rules. Calibrate `foi_past` so that `rel_sev` by age at malaria start matches `rel_sev` by age after the burn-in.

**Detectability** (used by RDTs in stage 3): high in C, lower in A, zero for placental-only infection. This is the source of false negatives in primigravidae.

Human infectiousness depends on peripheral infection only. Placental infection never feeds back into transmission.

### `placental` (connector)

States per woman:

- `ab` (float, 0–1): VAR2CSA antibody level (latent)
- `parasites` (bool): placental parasites present
- `pigment` (bool): haemozoin in placental tissue; once set, never reverts within a pregnancy
- `ti_parasites`: time the current placental infection began
- `dur_parasites`, `dur_parasites_pigment` (float): time with parasites, and with parasites plus pigment, this pregnancy (stage 2 harm uses these, not the histology label)

Histology class is derived: none (neither), acute (parasites, no pigment), chronic (parasites and pigment), past (pigment only). A re-infection over old pigment is therefore chronic.

Each step, among pregnant women:

1. **New infection → placental?** Newly infected pregnant women with gestation ≥ `ga_plac` (~13 weeks) and no current placental parasites establish placental infection with probability `p0 / (1 + (ab/ab50)^h)`.
2. **Pigment.** Set after `t_chronic` of continuous placental parasites.
3. **Clearance.** Placental parasites clear at rate `1/dur_plac` (untreated, months), **independent of peripheral recovery**.
4. **Boost.** On establishing placental infection: `ab += boost · (1 − ab)`. On a new infection in pregnancy that does not establish placentally: `ab += f · boost · (1 − ab)`; `f = 0` gives the placental-only variant.
5. **Waning.** All women: `ab *= exp(−ln2 · dt / ab_halflife)`.
6. **Delivery** (Pregnancy delivery callback): record the histology class, then reset the placental states.

`clear_placental(uids)` sets `parasites = False` and leaves `pigment` (acute → none, chronic → past): the hook for treatment and IPTp in stage 3, with its own (lower) efficacy than peripheral clearance.

**Initial antibody.** At malaria start, replay each woman's `n_pregnancies` past pregnancies with a per-pregnancy infection probability `q_past` and the birth interval implied by her age and gravidity; each infected past pregnancy goes through the Hill function with her running `ab`, then boost and wane. Calibrate `q_past` per setting so that `ab` by age at malaria start matches `ab` by age after the burn-in.

### `outcomes` (connector, stage 2 interface)

Harm accrues with exposure and does not reverse when infection clears.

- **Growth restriction / low birth weight:** per-step penalty to `FetalHealth.apply_growth_restriction` while placental parasites are present, larger with pigment. Steepest gravidity gradient.
- **Preterm birth:** mainly from clinical episodes (C), with a smaller placental term, via `FetalHealth.apply_timing_shift`.
- **Miscarriage:** from clinical infection only; mostly first-trimester, before placental infection can establish.
- **Stillbirth:** from placental and clinical infection.
- **Maternal severe malaria and death:** branch from C with `p_sev`; case fatality depends on treatment.
- **Maternal anaemia:** state on the mother, from peripheral and placental infection (plus stage 4 co-factors); flattest across gravidity.
- **Neonatal death:** from low birth weight and preterm birth.

Hooks into v4 `ss.Pregnancy` for miscarriage, stillbirth and neonatal death to be confirmed.

### `anc` and segments (stage 3 interface)

- **Segments:** agent attributes `wealth` (quintile) and `rural`, plus age and gravidity. Assigned at initialisation from survey distributions; newborns inherit their mother's. Results are stratified by segment.
- **Cascade input:** a table from Activity 1 with one row per segment × step (× dose), and columns `value`, `low`, `high`. Steps: timing of first ANC contact, number of visits, IPTp dose *k* given a visit, RDT at a visit, iron-folic acid, deworming, referral for complications.
- **`anc` process** per pregnant woman: draw first-contact gestational age and a visit schedule; at each visit apply IPTp (SP efficacy by `sp_resistance`: clears peripheral infection, calls `clear_placental`, gives a prophylactic window), RDT (sensitivity by detectability) with treatment, iron-folic acid and deworming.
- Products (SP, RDT) are separate from delivery, as in starsim.
- Optional: availability terms that respond to transmission (stock-outs at peaks).

### Co-factors (stage 4 interface)

Background anaemia (haemoglobin state), helminth co-infection (cleared by deworming), and an iron-malaria interaction (iron supplementation acting on `rel_sus` or `rel_sev`).

### Infants (stage 5 interface)

At delivery, the newborn's protection comes from the mother's *general* immunity (the base of her `rel_sev`), not from `placental.ab`, which is placenta-specific. Transferred protection = `(1 − rel_sev_mother) · eff_transfer · (1 − imp_transfer · parasites)`, decaying over a few months, applied as `ss.modify(malaria.rel_sev, multiply=..., where=infant_protected)`.

### `Scenarios` (container, prototype)

Burn in once, branch many times.

- **Burn-ins.** For each setting, K sims with different `rand_seed`, run with `sim.run(until=branch_date)` and saved: `burnins/<setting>_<k>.sim`. A sim keeps one seed for its whole life; there is no reseeding at the branch point. Stochastic spread comes from the K seeds.
- **Branching.** `Scenarios(burnins, scenarios)` loads each burn-in, applies each scenario's changes, and runs to `stop`, in parallel. Scenarios compared within a burn-in share its random streams.
- **Scenario changes.** Parameter changes after load work. Interventions are declared (inactive) in the burn-in and activated at the branch point; adding modules after load is not relied on. How a scenario is specified (`ss.modify`, an intervention, a par dict) is open.
- `dt` cannot change at the branch point (`Disease.dt` has no setter), so malaria runs at its final `dt` through its burn-in; the cost is paid once per burn-in.

Verified on v4 (re-checked 2026-10-15): continuation after save/load is bit-identical; parameter changes after load take effect; setting `rand_seed` after load has no effect; a saved 5,000-agent sim is 258 KB and loads in < 1 s; per-agent branch probabilities from a custom agent quantity work.

## Parameters (priors, to be calibrated)

| Par | Prior | Notes |
|---|---|---|
| `a`, `b`, `c` | 0.3/day, 0.5, 0.5 | Ross-Macdonald |
| `r_A` | 0.3 | Infectiousness of asymptomatic relative to clinical |
| `dur_clin`, `dur_asym` | 14 d, 150 d | Untreated clinical phase; asymptomatic clearance |
| `p_treat`, `dur_treat` | 0.5, 3 d | Background case management |
| `p_clin0` | 0.8 | Clinical probability in the fully naive |
| `sev_boost`, `sev_halflife` | 0.2, 5 y | General immunity per infection; waning |
| `p_sev0`, `k` | 0.05, 2 | Stage 2 |
| `intensity` | per setting | Calibrated to MAP *Pf*PR₂₋₁₀ |
| `rr_preg` | 1.5–2 | Reduced-form; see `Malaria` |
| `p0`, `ab50` | 0.7, 0.3 | Placental establishment vs antibody |
| `h` | 2 (fixed) | Fixed to limit the `boost`/`ab50` ridge |
| `boost` | 0.4 (0.3–0.5) | Per infection |
| `f` | 0.15 (fixed; 0 as switch) | Boost from non-placental infection in pregnancy |
| `ab_halflife` | 8 y (6–10) | Strongest control on the gravidity gradient |
| `ga_plac` | 13 weeks | Earliest placental establishment |
| `t_chronic` | 6 weeks | Parasites → pigment |
| `dur_plac` | 90 d | Untreated placental clearance |
| `p_loss` | 0.15 | Pregnancies not ending in birth, for initial gravidity |

**Joint calibration.** `p0`, `boost`, `ab50`, `ab_halflife` and the force of infection act on one observable (the gravidity gradient). A per-woman Monte Carlo of the mechanism (one infection per pregnancy, 3-year birth interval) gives G2/G1 and G3/G1 of 0.71 and 0.67 with a 3-year half-life and `boost = 0.5`; reaching the observed ≈ 0.6 and 0.4 needs a half-life of 6–10 years with `boost` 0.3–0.5. At high intensity, multiple infections per pregnancy raise the effective per-pregnancy boost, so per-infection values come out of calibration. Likewise `sev_boost`, `sev_halflife` and `p_clin0` are calibrated jointly with intensity against age patterns of clinical disease.

**Identifiability.** `ab` is latent and the targets are parasitaemia and histology. Expect a ridge between `boost` and `ab50`, broken only weakly by the `1 − ab` saturation. Fix `h` and `f`; calibrate `p0`, `boost`, `ab50`, `ab_halflife` jointly and check for the ridge.

## Testing

pytest in `tests/`:

- `rel_sev`: boost/wane arithmetic; clinical share falls as `rel_sev` falls; modifiers on `rel_sev` do not change its base
- **Immunity transition:** after lowering `intensity`, the clinical share of infections in adults rises over time
- **Asymptomatic reservoir:** at high intensity, infection prevalence in adults stays well above zero while the clinical share is low
- Placental boost/wane arithmetic; `p_plac` decreases monotonically with `ab`
- **Null gradient:** with `boost = 0` and `f = 0`, placental infection is flat across gravidity
- **Gradient check:** the placental mechanism alone reproduces the Monte Carlo gradients above
- Longer `ab_halflife` gives a steeper gradient
- **Birth interval:** at fixed gravidity, longer intervals give lower `ab`
- Pigment never reverts within a pregnancy; `clear_placental` maps acute → none and chronic → past
- No placental infection before `ga_plac`; placental clearance is independent of peripheral recovery
- **Initial immunity:** `rel_sev` and `ab` by age at malaria start ≈ by age after the burn-in
- Initial `n_pregnancies` matches cumulative ASFR × (1 + `p_loss`) by age
- Unit consistency: `Malaria` epi outputs agree within tolerance between a sim with a monthly unit and one with a daily unit
- Branching: continued run equals uninterrupted; a parameter change after load takes effect
- `sim.check_dt(factor=2)`: daily vs 2-day malaria
- Speed benchmark: 10k agents, monthly demography → 10-y malaria burn-in → 5-y projection; decides 1- vs 2-day `dt`

**Stage 1 success:**

- Clinical share of infections falls with age at high intensity, and the age pattern shifts older at low intensity.
- The 4-class histology distribution at delivery differs by gravidity in the observed direction at high intensity, within age bands, and the gradient is smaller at low intensity. Comparisons average over enough years that conception-month composition does not drive them.
- The projection step is fast enough for hundreds of scenarios.

**Calibration targets (later stages):** *Pf*PR₂₋₁₀ and clinical incidence by age per setting; peripheral parasitaemia at first ANC by gravidity (ANC surveillance); placental histology class by gravidity and age (e.g. Manhiça cohorts; Desai et al. 2007), using HIV-negative strata where reported; cascade coverage by segment from Activity 1.

## Limitations

- **No HIV.** HIV attenuates gravidity-dependent protection, and Mozambique and the Manhiça cohorts have high HIV prevalence. Fitting without HIV biases `boost` low; prefer HIV-negative strata.
- Antibodies block placental establishment only; they do not shorten placental infections.
- Immunity reduces disease, not parasite density or infection duration.
- Placental infection does not feed back into transmission.
- No postpartum parasite release.
- No link from `ab` to anti-VAR2CSA IgG titre (out of scope, decided 2026-10-14).

## Open questions

- Shape of the `Scenarios` API: how a scenario is specified
- Miscarriage, stillbirth and neonatal-death hooks in v4 `ss.Pregnancy`
- Segment attributes beyond wealth and residence, once Activity 1 defines the segments
- Whether to model postpartum parasite release
- Final malaria `dt` (benchmark)
