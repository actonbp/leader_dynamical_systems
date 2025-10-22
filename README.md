# Dynamic Leadership Identity Dynamics Toolkit

This folder contains the shareable code and documentation that accompany Chapter 5 (Leadership Identity at a Crossroads). It is self-contained so it can be dropped directly into a public repository.

## Contents

### scripts/
- `dynamic_states_walkthrough.R` – Step 1–3 pipeline (interpolation → residuals → GOLD derivatives) with simulation fallback and optional CSV input.
- `dynamic_states_walkthrough.qmd` – Quarto narrative of the same workflow.
- `dynamic_states_sim_demo.qmd` – Simulation-based demo with visual diagnostics (phase plots, distributions, scatter).
- `simulate_dynamic_identity.R` – Helper to simulate damped identity trajectories with event inputs.
- `setup_packages.R` – Convenience function for installing missing dependencies (uses `pak` when available).

### resources/
- `dynamic_states_method_guide.md` – When to use derivative proxies vs. ODE/ctsem approaches; interpretation guidance.
- `data_requirements.md` – Expected data schema and CSV example if providing your own time series.

### data/
- `example_identity_timeseries.csv` – Example data file with 2 participants across 3 weeks (7 days each)
- Place your own data here as `identity_timeseries.csv` matching the schema in `resources/data_requirements.md`

## Quick Start
1. Optionally run `source("scripts/setup_packages.R"); ensure_packages(c("tidyverse","here","nlme","imputeTS","doremi","lmerTest","sessioninfo"))`.
2. If you have your own data, save it as `data/identity_timeseries.csv` matching `resources/data_requirements.md`.
3. Run `scripts/dynamic_states_walkthrough.R` or render `scripts/dynamic_states_walkthrough.qmd` for the core analysis.
4. Render `scripts/dynamic_states_sim_demo.qmd` for the simulated illustration and plots.

## License
This toolkit inherits the root repository's license. Update this README if you distribute it separately.
