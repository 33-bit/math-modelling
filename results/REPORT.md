# Member 4 Report: Experiments and Sensitivity Analysis

## Scope

This report covers Member 4 responsibilities: percolation probability experiments,
propagation speed, epidemic curves, secondary-infection distribution, and sensitivity
checks over population size, superspreader fraction, density, and random seeds.

## Reproducibility

- Simulator: `MonteCarloSIRSimulator`
- Random seeds per setting: `10`
- Baseline population: `N = 477`
- Baseline density: `4.77`
- Baseline superspreader fraction for strong and hub models: `lambda = 0.20`
- Percolation event definition: final attack rate at least `50%`
- Generated outputs: `summary_metrics.csv`, `percolation_probability.csv`,
  `propagation_speed.csv`, `epidemic_curves.csv`,
  `secondary_distribution.csv`, and `sensitivity_summary.csv`

## Baseline Results

| Model | Mean attack rate | Percolation probability | Mean duration | Mean propagation speed | Mean secondary infections |
| --- | ---: | ---: | ---: | ---: | ---: |
| normal | 0.097 | 0.000 | 10.800 | 0.160 | 0.858 |
| strong | 0.792 | 0.800 | 15.300 | 0.281 | 0.798 |
| hub | 0.794 | 0.800 | 10.600 | 0.385 | 0.798 |

The hub model produces the fastest baseline spread. The strong-infectiousness model
also increases final epidemic size relative to the normal model, but its propagation
speed is lower because superspreaders still infect within the shorter normal range.

## Percolation and Density Sweep

| Model | Density | Percolation probability | Mean attack rate | Mean propagation speed |
| --- | ---: | ---: | ---: | ---: |
| hub | 2.50 | 0.20 | 0.28 | 0.17 |
| hub | 3.50 | 0.70 | 0.65 | 0.30 |
| hub | 4.77 | 0.80 | 0.79 | 0.39 |
| hub | 6.00 | 0.90 | 0.90 | 0.51 |
| hub | 8.00 | 0.90 | 0.90 | 0.53 |
| normal | 2.50 | 0.00 | 0.01 | 0.09 |
| normal | 3.50 | 0.00 | 0.03 | 0.08 |
| normal | 4.77 | 0.00 | 0.10 | 0.16 |
| normal | 6.00 | 0.80 | 0.73 | 0.16 |
| normal | 8.00 | 1.00 | 0.98 | 0.29 |
| strong | 2.50 | 0.00 | 0.07 | 0.19 |
| strong | 3.50 | 0.50 | 0.50 | 0.15 |
| strong | 4.77 | 0.80 | 0.79 | 0.28 |
| strong | 6.00 | 0.90 | 0.90 | 0.39 |
| strong | 8.00 | 0.90 | 0.90 | 0.46 |

Percolation probability increases with density because each infectious individual has
more neighbors inside the infection radius. Hub superspreaders reach the high
percolation regime earlier than the strong model because their effective infection range
is larger.

## Sensitivity Findings

- Increasing `lambda` raises attack rate and speed for both superspreader models.
- At baseline `N` and density, the strongest hub speed in the sweep occurs at
  `lambda = 0.40` with mean speed
  `0.675`.
- At baseline `N` and density, the strongest strong-model speed in the sweep occurs at
  `lambda = 0.40` with mean speed
  `0.493`.
- Changing `N` while keeping density fixed gives similar attack-rate trends, so the
  model behavior is mainly controlled by density and superspreader fraction rather than
  raw population size alone.
- Lower density is the most important robustness stress test: it reduces outbreak
  probability and makes outcomes more seed-dependent.

## Figures

- `plots/percolation_probability.png`
- `plots/epidemic_curves.png`
- `plots/secondary_distribution.png`
- `plots/sensitivity_lambda_attack_rate.png`

## Notes for Final Group Report

The experiments support the qualitative conclusion that a small fraction of
superspreaders can sharply increase epidemic size and speed. The hub model is more
dangerous than the strong-infectiousness model under the same `lambda` because it
changes the contact geometry, not only the probability of infection after contact.
