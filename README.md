# mozambique_pregnancy

Agent-based model of malaria in pregnancy in Mozambique, built on [STIsim](https://github.com/starsimhub/stisim) / [Starsim](https://github.com/starsimhub/starsim).

The motivating question is how malaria immunity in pregnancy depends on reproductive history. Women acquire antibodies against placental malaria with successive pregnancies, so first-time mothers are most vulnerable and protection rises with each pregnancy. Capturing this requires a realistic distribution of pregnancies by age and over time, which this model provides as a first step.

**Status:** demographics only — births via `ss.Pregnancy`, background deaths, and UN WPP age structure, compared against WPP. Malaria and immunity are not yet included.

## Repo layout

- `model.py` — sim builder (`make_sim`); run as a script to simulate 1985–2030, save `results/sim.df`, and make figures
- `analyzers.py` — `PregnancyByAge`: women, currently pregnant, new pregnancies, births, and parity by 5-year maternal age group; `AgePyramid`: population by sex and 5-year age group
- `plot_demographics.py` — population/fertility vs UN WPP (`figures/demographics.png`), parity by age (`figures/parity.png`), and age pyramids by sex (`figures/age_pyramids.png`)
- `data/` — UN WPP 2024 files for Mozambique: age distribution (1985), age-specific fertility, age/sex-specific mortality, total population, population by age and sex (1985–2025, every 10 years)

## Usage

```bash
python model.py               # run sim + make figures
python plot_demographics.py   # replot from results/sim.df
```

Requires stisim >= 1.7.1 (starsim 3.6).

## Data

Files in `data/` were downloaded once from the [UN Data Portal API](https://population.un.org/dataportalapi/index.html) via stisim's downloader and are committed, so runs need no network access. To refresh them (e.g. for a different start year), run `sti.Sim(demographics='mozambique', start=...)` without `datafolder`; stisim downloads to `stisim/data/files/` (needs a UN token at `stisim/data/auth_key.txt`), then copy the files into `data/`.

## Known limitations

- Parity is 0 for everyone at sim start, so it is only complete for cohorts aged ≤ 15 in 1985 (dashed lines in `figures/parity.png`).
- Primary infertility defaults to 0, so almost no women aged 30+ are nulliparous.
- Births at ages 15–19 are below WPP: `min_age=15` applies at conception, so no births occur before ~15.75.
- No migration, so the 1990–93 population plateau in WPP is not reproduced.
