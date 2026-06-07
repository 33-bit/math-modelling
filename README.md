# Monte Carlo SIR Simulator with Superspreaders

Spatial Monte Carlo SIR simulator for studying the effects of superspreaders in an
epidemic. The code supports the normal model, strong-infectiousness superspreaders,
and hub superspreaders, then runs the paper-matched experiment suite.

## Project Structure

- `config.py`: constants, default parameters, and model names.
- `simulator.py`: core SIR simulator and superspreader infection probability models.
- `demo_run.py`: small sanity-check run for the three model types.
- `experiments.py`: Member 4 experiment pipeline.
- `results/`: generated clean CSV files, plots, report, and output README.
- `requirements.txt`: Python dependencies.

## Installation

Use Python 3.12 or newer if available.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

If `.venv` already exists, only reinstall dependencies when needed:

```bash
.venv/bin/pip install -r requirements.txt
```

## Run a Demo

```bash
.venv/bin/python demo_run.py
```

The demo runs:

- `normal`
- `strong` with `lambda_ss = 0.2`
- `hub` with `lambda_ss = 0.2`

It prints total infected count, duration, early epidemic curve values, and the largest
secondary-infection counts.

## Run Member 4 Experiments

```bash
.venv/bin/python experiments.py
```

By default, this runs 1000 random seeds per parameter setting to match the paper's
Monte Carlo averaging and writes clean outputs to `results/`. For quick coursework
checks, pass a smaller seed count explicitly.

Useful options:

```bash
.venv/bin/python experiments.py --seeds 30
.venv/bin/python experiments.py --seeds 1000 --seed-offset 100000
.venv/bin/python experiments.py --seeds 30 --output-dir results/quick_check
.venv/bin/python experiments.py --max-steps 300
.venv/bin/python experiments.py --output-dir results/member4_large
```

## Member 4 Outputs

Generated files:

- `results/README.md`: output-specific notes.
- `results/summary_metrics.csv`: one row per replicate.
- `results/baseline_summary.csv`: averaged paper fixed-density comparison.
- `results/percolation_probability.csv`: percolation probability by density.
- `results/critical_density.csv`: density where percolation probability crosses 0.5.
- `results/propagation_speed.csv`: speed summaries across experiment groups.
- `results/front_distance.csv`: infection front distance over time by lambda.
- `results/epidemic_curves.csv`: mean new, active, and cumulative infections.
- `results/secondary_distribution.csv`: secondary-infection distribution.
- `results/infection_probability_functions.csv`: model-definition values for paper
  Fig. 1 and Fig. 2.
- `results/sars_singapore_secondary_patients.csv`: empirical/reconstructed SARS
  Singapore secondary-patient distribution.
- `results/sars_singapore_epidemic_curve.csv`: approximate six-day SARS Singapore
  epidemic-curve bins for paper Fig. 15 comparison.
- `results/sars_epidemic_model_curves.csv`: model curves generated with the paper's
  SARS comparison settings.
- `results/sensitivity_summary.csv`: paper lambda-sweep summaries for velocity and
  attack-rate plots.
- `results/plots/*.png`: figures for analysis or presentation.

Percolation is defined as an outbreak reaching the top band of the spatial system.
Propagation speed is estimated from the slope of the infection front radius over time.
Percolation uses the paper's `L = 10 r0` setup and varies `N = 150..900`.
Propagation and paper epidemic-curve plots use `N = 637` in the same box, giving
`rho * pi * r0^2 = 20.012`, with `lambda = 0.2` for the superspreader models.
Route and secondary-distribution plots use `N = 477`, giving
`rho * pi * r0^2 = 14.985`, with `lambda = 0.2`. The SARS comparison uses
`N = 477`, giving `rho * pi * r0^2 = 14.985`, with `lambda = 0.4` and
`1 timestep = 6 days`.

## Basic Usage

```python
from simulator import MonteCarloSIRSimulator

sim = MonteCarloSIRSimulator(
    N=477,
    model="hub",
    lambda_ss=0.2,
    seed=42,
)

result = sim.run()

print(result.total_infected)
print(result.duration)
print(result.new_infections_per_step)
print(result.secondary_counts)
```

## Supported Models

- `model="normal"`: standard distance-based infection.
- `model="strong"`: superspreaders infect with stronger probability within the normal
  infection radius.
- `model="hub"`: superspreaders use a larger effective infection radius.

For `strong` and `hub`, `lambda_ss` is the fraction of superspreaders in the population.

## Reproducibility Notes

- All batch outputs record the random seed used by each replicate; use `--seed-offset`
  to generate a fresh reproducible seed block.
- Density is controlled as `N / L^2`; the experiment script computes `L` from the target
  density.
- Distances wrap horizontally, while the vertical axis remains open so bottom-to-top
  percolation and front-distance plots are meaningful.
- Existing generated CSV, Markdown, and plot files in the selected output directory are
  overwritten so the directory contains only the current clean run.
