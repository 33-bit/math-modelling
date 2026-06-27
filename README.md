# Monte Carlo SIR Simulator with Superspreaders

This project implements a spatial Monte Carlo SIR simulator for reimplementation
experiments on epidemic spreading with superspreaders. It supports the normal
infection model, strong-infectiousness superspreaders, and hub superspreaders, then
generates CSV summaries and figures for the full experiment suite.

The implementation covers all 15 numbered paper-figure workflows plus four
supplemental comparison figures. The generated outputs are reproducible from fixed
random seeds, but they should be described as a reimplementation rather than an
exact numerical copy of every published curve.

## Project Structure

- `epidemic_model/config.py`: shared constants, state labels, defaults, and model
  names.
- `epidemic_model/simulator.py`: core Monte Carlo SIR simulator and infection
  probability kernels.
- `epidemic_model/experiment_pipeline/`: scenario builders, metric collection,
  CSV writers, plot generation, and saved-output reconstruction.
- `experiments.py`: thin command-line entry point for the experiment pipeline.
- `config.py`, `simulator.py`, `experiment_pipeline/`: backward-compatible import
  wrappers for older code that imports from the project root.
- `demo_run.py`: small sanity-check script for the three model variants.
- `realtime_visualization.py`: interactive infection-route visualization.
- `results/`: committed generated CSV files and PNG plots.
- `tests/`: regression tests for formulas, boundary behavior, and experiment helpers.
- `requirements.txt`: Python dependencies.

## Installation

Use Python 3.12 or newer if available.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

If `.venv` already exists, reinstall dependencies only when needed:

```bash
.venv/bin/pip install -r requirements.txt
```

## Quick Demo

```bash
.venv/bin/python demo_run.py
```

The demo runs:

- `normal`
- `strong` with superspreader fraction `0.2`
- `hub` with superspreader fraction `0.2`

It prints total infected count, duration, early epidemic curve values, and the
largest secondary-infection counts.

## Realtime Visualization

Run the interactive route-of-infection viewer with:

```bash
.venv/bin/python realtime_visualization.py
```

The window animates the simulation step by step and keeps the infection tree,
active/recovered states, and epidemic curve updated in realtime. Use the controls
to change `N`, `lambda`, `model`, `L`, `r0`, `w0`, `gamma`, and playback speed, then
press `Apply` to start a new run from those settings.

You can also set initial values from the command line:

```bash
.venv/bin/python realtime_visualization.py --N 477 --model hub --lambda 0.2 --L 10 --seed 42
```

## Full Experiment Run

The committed reimplementation result set was generated with:

```bash
python3 experiments.py --seeds 1000 --seed-offset 100000 --max-steps 200 --jobs 1
```

This is a 284,000-simulation run: 284 experiment configurations times the clean
seed block `100000..100999`. Every replicate records its seed in
`results/summary_metrics.csv`.

For a quick all-configuration smoke run, use a smaller seed count and a temporary
output directory:

```bash
python3 experiments.py --seeds 1 --seed-offset 200000 --max-steps 20 --jobs 1 --output-dir /tmp/experiment-smoke
```

To regenerate plots from existing CSV files without rerunning the full simulation
batch:

```bash
python3 experiments.py --plot-only
```

`--plot-only` regenerates plots from CSV outputs and reruns only the three exact
infection-route configurations stored in `results/route_plot_selections.csv`.

Useful options:

```bash
python3 experiments.py --seeds 30
python3 experiments.py --seeds 1000 --seed-offset 100000
python3 experiments.py --seeds 30 --output-dir results/quick_check
python3 experiments.py --max-steps 300
python3 experiments.py --jobs 4 --chunksize 20
```

## Running One Experiment Group

List available groups with:

```bash
python3 experiments.py --list-experiments
```

Run one group with `--experiment`:

```bash
python3 experiments.py --experiment percolation --seeds 30
python3 experiments.py --experiment routes --seeds 100 --max-steps 200
python3 experiments.py --experiment sars --seeds 1000 --seed-offset 100000
```

You can combine groups by repeating the flag:

```bash
python3 experiments.py --experiment routes --experiment sars --seeds 100
```

When `--output-dir` is omitted, the full run writes to `results/`, while selected
runs write to group-specific folders such as `results/percolation/` or
`results/selected_routes_sars/`. This keeps a partial run from overwriting the
committed full result set by accident.

`--plot-only` also respects `--experiment`:

```bash
python3 experiments.py --experiment percolation --plot-only
```

Running `experiments.py` overwrites CSV, Markdown, and PNG files in the selected
output directory so that the directory contains only the current clean run. The
project intentionally does not keep a separate `results/README.md`; this file is
the single project documentation entry point.

## Generated Outputs

CSV files in `results/`:

- `summary_metrics.csv`: one row per simulation replicate.
- `baseline_summary.csv`: averaged fixed-density comparison.
- `percolation_probability.csv`: percolation probability by model and density.
- `critical_density.csv`: density where percolation probability crosses 0.5.
- `critical_density_reference_curves.csv`: analytical Eq. 3-5 reference values used
  for Fig. 5.
- `front_speed.csv`: front-speed summaries across experiment groups.
- `front_distance.csv`: mean infection front distance over time by superspreader
  fraction.
- `epidemic_curves.csv`: mean new, active, and cumulative infections.
- `secondary_distribution.csv`: secondary-infection distribution from simulations.
- `infection_probability_functions.csv`: model-definition values for Figs. 1-2.
- `sars_singapore_secondary_patients.csv`: reconstructed SARS Singapore secondary
  patient distribution for Fig. 14.
- `sars_singapore_epidemic_curve.csv`: approximate six-day SARS Singapore epidemic
  curve bins for Fig. 15.
- `sars_epidemic_model_curves.csv`: simulated model curves for the SARS comparison.
- `sensitivity_summary.csv`: superspreader-fraction sweep summaries for
  front-speed and case-fraction plots.
- `route_plot_selections.csv`: exact seeds and configurations selected for the
  infection-route figures.

PNG files in `results/plots/`:

- `fig01_infection_probability_strong.png`
- `fig02_infection_probability_hub.png`
- `fig03_percolation_probability_strong.png`
- `fig04_percolation_probability_hub.png`
- `fig05_critical_density.png`
- `fig06_front_distance_strong.png`
- `fig07_velocity_vs_superspreader_fraction.png`
- `fig08_epidemic_curves.png`
- `fig09_infection_route_strong.png`
- `fig10_infection_route_hub.png`
- `fig11_infection_route_normal.png`
- `fig12_secondary_distribution_normal.png`
- `fig13_secondary_distribution_superspreaders.png`
- `fig14_sars_secondary_patients.png`
- `fig15_sars_epidemic_curve_comparison.png`
- `fig_extra_front_distance_hub.png`
- `fig_extra_percolation_probability_comparison.png`
- `fig_extra_secondary_distribution_comparison.png`
- `fig_extra_sensitivity_superspreader_fraction_case_fraction.png`

## Model Details

The simulator keeps individuals in a square box. By default, the domain side
length is `10 r0`, `w0 = 1`, `gamma = 1`, and the initial infected individual is
placed at the bottom center of the domain.

Supported infection models:

- `normal`: standard distance-based infection probability
  `w(r) = w0 * (1 - r / r0)^2` for `r <= r0`.
- `strong`: superspreaders infect with probability `w0` inside the normal radius
  `r0`.
- `hub`: superspreaders use a larger radius `sqrt(6) * r0` with the same quadratic
  decay form.

For `strong` and `hub`, the superspreader fraction is the share of the population
assigned the superspreader infection kernel.

## Experiment Choices

Percolation is defined as infection reaching the top band of the spatial system.
The paper states periodic boundary conditions while also defining bottom-to-top
percolation from an initial case at the bottom. This implementation uses the
operational interpretation implied by the paper figures: horizontal wrapping with
an open vertical direction. Full vertical wrapping would make the top-reaching
criterion ambiguous.

The percolation sweep uses a domain side length of `10 r0` and population sizes
from `150..900`, extended with `50, 75, 100, 125`. These extra low-density points
allow the hub-model critical density near `rho * pi * r0^2 = 3.2` to be
interpolated rather than clipped to the first sampled point.

For Fig. 5, analytical reference rows use:

```text
R0 = [superspreader_fraction + (1 - superspreader_fraction) / 6] * rho * pi * r0^2
rho_c * pi * r0^2 = Rc / [superspreader_fraction + (1 - superspreader_fraction) / 6]
```

with `Rc = 4.5` for the strong model and `Rc = 3.2` for the hub model.

Propagation and epidemic-curve plots use population size `637` in the same box,
giving `rho * pi * r0^2 = 20.012`, with superspreader fraction `0.2` for the
superspreader models. Route and secondary-distribution plots use population size
`477`, giving `rho * pi * r0^2 = 14.985`, with superspreader fraction `0.2`.

The SARS comparison uses population size `477`, `rho * pi * r0^2 = 14.985`,
superspreader fraction `0.4`, and `1 timestep = 6 days`. Fig. 15 compares normal,
strong-superspreader, and hub-superspreader model curves.

## SARS Data Notes

The SARS observations are not raw official case-level data.

- Fig. 14 uses reconstructed/digitized CDC MMWR secondary-infection frequencies and
  accounts for all 201 probable cases.
- Fig. 15 uses approximate six-day case bins digitized from the published Singapore
  SARS epidemic curve.

The original digitized values are encoded in
`epidemic_model/experiment_pipeline/settings.py` as `SARS_SECONDARY_CASE_FREQUENCIES`
and `SARS_SIX_DAY_CASE_COUNTS`. The generated CSV copies are
`results/sars_singapore_secondary_patients.csv` and
`results/sars_singapore_epidemic_curve.csv`. Both the CSV source column and the
Fig. 15 plot label mark the epidemic-curve data as approximate/digitized.

## Tests

Run the regression suite with:

```bash
python3 -m unittest discover -s tests -v
```

The tests cover infection probability formulas, optional vertical wrapping, the
top-band percolation rule, critical-density interpolation, Eq. 3-5 reference rows,
SARS secondary-patient totals, route-selection restoration, legacy CSV header
compatibility, and output schemas.
