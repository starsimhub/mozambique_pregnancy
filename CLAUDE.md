# CLAUDE.md

Malaria-in-pregnancy ABM for Mozambique on starsim v4. Read `README.md` and the design spec first.

## State (paused 2026-10-14)

- **Done:** demographics-only sim (`model.py`, still on `sti.Sim`), `analyzers.py` (`PregnancyByAge`, `AgePyramid`), `plot_demographics.py`, WPP data in `data/`.
- **Design:** `docs/superpowers/specs/2026-10-14-malaria-in-pregnancy-design.md` — stage 1 spec, revised after external review. Awaiting the user's final approval.
- **Next:** confirm approval, then write the stage-1 implementation plan (superpowers:writing-plans). No malaria code exists yet.

## Environment

- conda env `starsim`: editable starsim on branch `v4` (4.0.0), editable stisim 1.7.1.
- Demography data is committed; runs need no network. Re-downloading needs a UN Data Portal token (see `README.md` § Data).

## v4 gotchas (verified by probe)

- `ss.Pregnancy` and `FetalHealth` are agent-only: use `solver='abm'`.
- **Bare-number rates take the sim's time unit** — always type rates (`ss.perday`, `ss.days`).
- Modules can have their own `dt` and `start`; `Disease.dt` cannot be changed after a run starts.
- `sim.run(until=...)` + `sim.save()` + `sc.load()` continues bit-identically; parameter changes after load take effect; setting `rand_seed` after load is ignored (by design: one seed per sim).

## Conventions

- Names terse, starsim-style: `Malaria` (module), `placental`, `outcomes` (connectors), `Scenarios` (container).
- Burn in once (K seeds per intensity), branch many scenarios from the saved burn-ins.
- Don't push without asking.
