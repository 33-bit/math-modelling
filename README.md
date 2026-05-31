# Monte Carlo SIR Simulator with Superspreaders

Spatial Monte Carlo SIR simulator for studying the effects of superspreaders in an
epidemic. The code supports the normal model, strong-infectiousness superspreaders,
and hub superspreaders, then runs the Member 4 experiment suite for robustness and
sensitivity analysis.

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

By default, this runs 10 random seeds per parameter setting and writes clean outputs to
`results/`.

Useful options:

```bash
.venv/bin/python experiments.py --seeds 30
.venv/bin/python experiments.py --max-steps 300
.venv/bin/python experiments.py --output-dir results/member4_large
```

## Member 4 Outputs

Generated files:

- `results/REPORT.md`: written experiment report for Member 4.
- `results/README.md`: output-specific notes.
- `results/summary_metrics.csv`: one row per replicate.
- `results/baseline_summary.csv`: averaged baseline comparison.
- `results/percolation_probability.csv`: percolation probability by density.
- `results/critical_density.csv`: density where percolation probability crosses 0.5.
- `results/propagation_speed.csv`: speed summaries across experiment groups.
- `results/front_distance.csv`: infection front distance over time by lambda.
- `results/epidemic_curves.csv`: mean new, active, and cumulative infections.
- `results/secondary_distribution.csv`: secondary-infection distribution.
- `results/sensitivity_summary.csv`: sensitivity over `N`, `lambda_ss`, density,
  and random seeds.
- `results/plots/*.png`: figures for report or presentation.

Percolation is defined as an outbreak reaching the top band of the spatial system.
Propagation speed is estimated from the slope of the infection front radius over time.

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

- All batch outputs record the random seed used by each replicate.
- Density is controlled as `N / L^2`; the experiment script computes `L` from the target
  density.
- Distances wrap horizontally, while the vertical axis remains open so bottom-to-top
  percolation and front-distance plots are meaningful.
- Existing generated CSV, Markdown, and plot files in the selected output directory are
  overwritten so the directory contains only the current clean run.
