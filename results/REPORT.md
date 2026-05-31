# Member 4 Report: Experiments and Sensitivity Analysis

## Scope

This report covers Member 4 responsibilities: running numerical experiments,
measuring percolation probability, estimating propagation speed, generating epidemic
curves, measuring secondary-infection distributions, and checking sensitivity over
population size, superspreader fraction, density, and random seeds. The focus is the
experiment/result part of the reference paper rather than the mathematical-model
definition plots.

## Reproducibility

- Simulator: `MonteCarloSIRSimulator`
- Random seeds per setting: `30`
- Baseline population: `N = 477`
- Baseline density: `4.77`
- Baseline superspreader fraction for strong and hub models: `lambda = 0.20`
- Percolation event definition: infection reaches the top band of the two-dimensional system
- Generated outputs: `summary_metrics.csv`, `percolation_probability.csv`,
  `critical_density.csv`, `propagation_speed.csv`, `front_distance.csv`,
  `epidemic_curves.csv`, `secondary_distribution.csv`, and `sensitivity_summary.csv`

## How to Read the Metrics

- **Attack rate** is the final fraction of the population that became infected at
  least once. A higher attack rate means a larger outbreak.
- **Percolation probability** is the fraction of simulation runs in which the infection
  reaches the top band of the spatial system. This measures whether the disease can
  spread across the system, not only whether it infects many people locally.
- **Propagation speed** is estimated from the slope of the infection-front distance
  over time. A larger value means the epidemic wave travels faster through space.
- **Duration** is the number of time steps until no infected individuals remain.
- **Mean secondary infections** is the average number of people infected by one infected
  individual. A heavy tail in this distribution indicates superspreading events.
- **Critical density** is the approximate density where percolation probability reaches
  `0.5`. A lower critical density means the model can produce system-wide outbreaks
  even when the population is more sparse.
- The CSV files store raw density `rho = N / L^2`. The percolation and critical-density
  plots use the paper-style normalized density `rho * pi * r0^2` on the axis. Since this
  project uses `r0 = 1`, the plotted density is the raw density multiplied by `pi`.
- The percolation and velocity sweeps use
  `lambda = 0.0, 0.2, 0.4, 0.6, 0.8, 1.0`. The baseline comparison still uses
  `lambda = 0.20` for the two superspreader models.

## Baseline Results

| Model | Mean attack rate | Percolation probability | Mean duration | Mean propagation speed | Mean secondary infections |
| --- | ---: | ---: | ---: | ---: | ---: |
| normal | 0.085 | 0.033 | 9.133 | 0.118 | 0.564 |
| strong | 0.519 | 0.533 | 15.700 | 0.261 | 0.744 |
| hub | 0.594 | 0.600 | 10.733 | 0.449 | 0.748 |

The hub model produces the fastest baseline spread. The strong-infectiousness model
also increases final epidemic size relative to the normal model, but its propagation
speed is lower because superspreaders still infect within the shorter normal range.

The baseline results show a clear ordering of epidemic severity:
`hub > strong > normal`. In the normal model, only about
`0.085` of the population is infected on average, and
percolation probability remains low at `0.033`.
In the strong model, the attack rate increases to about
`0.519`, and percolation probability reaches
`0.533`. In the hub model, the attack rate rises to
about `0.594`, percolation probability reaches
`0.600`, and propagation speed is higher than the
strong-model speed. This supports the main claim that superspreaders do not only
increase the final outbreak size; they also change how quickly the epidemic moves
through space.

## Percolation and Density Sweep

| Model | lambda | Density | Percolation probability | Mean attack rate | Mean propagation speed |
| --- | ---: | ---: | ---: | ---: | ---: |
| hub | 0.20 | 0.50 | 0.00 | 0.00 | 0.02 |
| hub | 0.20 | 0.75 | 0.00 | 0.00 | 0.03 |
| hub | 0.20 | 1.00 | 0.00 | 0.01 | 0.08 |
| hub | 0.20 | 1.25 | 0.00 | 0.01 | 0.08 |
| hub | 0.20 | 1.50 | 0.00 | 0.01 | 0.12 |
| hub | 0.20 | 2.00 | 0.00 | 0.03 | 0.11 |
| hub | 0.20 | 2.50 | 0.10 | 0.10 | 0.19 |
| hub | 0.20 | 3.00 | 0.30 | 0.29 | 0.20 |
| hub | 0.20 | 3.50 | 0.37 | 0.37 | 0.27 |
| hub | 0.20 | 4.00 | 0.50 | 0.48 | 0.32 |
| hub | 0.20 | 4.77 | 0.60 | 0.59 | 0.45 |
| hub | 0.20 | 5.50 | 0.70 | 0.70 | 0.56 |
| hub | 0.20 | 6.50 | 0.87 | 0.86 | 0.68 |
| hub | 0.20 | 8.00 | 0.87 | 0.87 | 0.70 |
| normal | 0.00 | 0.50 | 0.00 | 0.00 | 0.00 |
| normal | 0.00 | 0.75 | 0.00 | 0.00 | 0.00 |
| normal | 0.00 | 1.00 | 0.00 | 0.00 | 0.00 |
| normal | 0.00 | 1.25 | 0.00 | 0.00 | 0.00 |
| normal | 0.00 | 1.50 | 0.00 | 0.00 | 0.01 |
| normal | 0.00 | 2.00 | 0.00 | 0.00 | 0.03 |
| normal | 0.00 | 2.50 | 0.00 | 0.01 | 0.04 |
| normal | 0.00 | 3.00 | 0.00 | 0.01 | 0.07 |
| normal | 0.00 | 3.50 | 0.00 | 0.01 | 0.07 |
| normal | 0.00 | 4.00 | 0.00 | 0.01 | 0.09 |
| normal | 0.00 | 4.77 | 0.03 | 0.08 | 0.12 |
| normal | 0.00 | 5.50 | 0.17 | 0.17 | 0.16 |
| normal | 0.00 | 6.50 | 0.57 | 0.53 | 0.19 |
| normal | 0.00 | 8.00 | 0.83 | 0.81 | 0.28 |
| strong | 0.20 | 0.50 | 0.00 | 0.00 | 0.04 |
| strong | 0.20 | 0.75 | 0.00 | 0.00 | 0.04 |
| strong | 0.20 | 1.00 | 0.00 | 0.01 | 0.06 |
| strong | 0.20 | 1.25 | 0.00 | 0.01 | 0.08 |
| strong | 0.20 | 1.50 | 0.00 | 0.01 | 0.08 |
| strong | 0.20 | 2.00 | 0.00 | 0.01 | 0.09 |
| strong | 0.20 | 2.50 | 0.00 | 0.02 | 0.12 |
| strong | 0.20 | 3.00 | 0.00 | 0.06 | 0.17 |
| strong | 0.20 | 3.50 | 0.07 | 0.16 | 0.17 |
| strong | 0.20 | 4.00 | 0.30 | 0.28 | 0.17 |
| strong | 0.20 | 4.77 | 0.53 | 0.52 | 0.26 |
| strong | 0.20 | 5.50 | 0.70 | 0.69 | 0.36 |
| strong | 0.20 | 6.50 | 0.87 | 0.86 | 0.46 |
| strong | 0.20 | 8.00 | 0.87 | 0.87 | 0.51 |

Percolation probability increases with density because each infectious individual has
more neighbors inside the infection radius. Hub superspreaders reach the high
percolation regime earlier than the strong model because their effective infection range
is larger.

At the baseline superspreader fraction `lambda = 0.20`, the hub model
starts percolating at lower density than the strong model. This means the hub mechanism
lowers the density barrier for large outbreaks. In practical terms, a hub-like
superspreader can connect distant local clusters, so the disease crosses the system
even when ordinary local transmission would still die out.

## Critical Density

| Model | lambda | Critical density |
| --- | ---: | ---: |
| strong | 0.00 | 6.33 |
| strong | 0.20 | 4.66 |
| strong | 0.40 | 3.25 |
| strong | 0.60 | 2.64 |
| strong | 0.80 | 1.97 |
| strong | 1.00 | 1.62 |
| hub | 0.00 | 6.33 |
| hub | 0.20 | 4.00 |
| hub | 0.40 | 2.40 |
| hub | 0.60 | 1.88 |
| hub | 0.80 | 1.25 |
| hub | 1.00 | 1.06 |

The critical-density plot summarizes the percolation curves into one threshold value.
As `lambda` increases, both models need less population density to percolate. The plot
also includes reference critical curves based on the paper's `R0 = Rc` argument, so the
markers show simulation estimates while the solid and dashed lines show the theoretical
trend. Lower critical density means the outbreak is more robust under sparse conditions,
so this result again shows that hub superspreaders are more dangerous than merely
stronger infectious individuals.

## Sensitivity Findings

- Increasing `lambda` raises attack rate and speed for both superspreader models.
- At baseline `N` and density, the strongest hub speed in the sweep occurs at
  `lambda = 1.00` with mean speed
  `1.193`.
- At baseline `N` and density, the strongest strong-model speed in the sweep occurs at
  `lambda = 1.00` with mean speed
  `0.766`.
- Changing `N` while keeping density fixed gives similar attack-rate trends, so the
  model behavior is mainly controlled by density and superspreader fraction rather than
  raw population size alone.
- Lower density is the most important robustness stress test: it reduces outbreak
  probability and makes outcomes more seed-dependent.

The sensitivity results show that `lambda` and density are the dominant controls. When
`lambda` increases, there are more superspreaders, so both attack rate and propagation
speed generally increase. When density increases, individuals have more nearby contacts,
so local clusters are easier to connect. Changing `N` while keeping density fixed has a
smaller effect because the average local neighborhood structure is similar.

## Figure Guide and Meaning

| Figure | What it shows | How to interpret it | Main result |
| --- | --- | --- | --- |
| `plots/percolation_probability.png` | Percolation probability versus normalized density `rho * pi * r0^2` for normal, strong, and hub models at baseline `lambda`. | Curves farther left percolate at lower density. | Hub percolates first, strong second, normal last. |
| `plots/percolation_probability_strong.png` | Strong-model percolation probability for `lambda = 0.0` to `1.0`, using normalized density on the x-axis. | Increasing `lambda` shifts the point cloud left and upward. | More strong superspreaders make system-wide spread possible at lower density. |
| `plots/percolation_probability_hub.png` | Hub-model percolation probability for `lambda = 0.0` to `1.0`, using normalized density on the x-axis. | The left shift is stronger than in the strong model. | Hub superspreaders reduce the percolation threshold more sharply. |
| `plots/critical_density.png` | Approximate normalized density `rho_c * pi * r0^2` where percolation probability reaches `0.5`, plus paper-style `R0 = Rc` reference curves. | Lower values mean easier system-wide spread; markers are simulation, lines are reference curves. | Critical density decreases with `lambda`, and the hub threshold is lower than the strong threshold. |
| `plots/front_distance_strong.png` | Infection-front distance `rf` over time for the strong model. | Steeper growth means faster spatial propagation; a plateau means spread has stopped. | Larger `lambda` makes the front travel farther and faster. |
| `plots/front_distance_hub.png` | Infection-front distance `rf` over time for the hub model. | Compare curve height and slope across `lambda`. | Hub spreading reaches far distances earlier than the strong model. |
| `plots/velocity_vs_lambda.png` | Propagation speed versus superspreader fraction. | Higher speed means the epidemic wave crosses space faster. | Hub speed is consistently higher than strong speed, especially at large `lambda`. |
| `plots/epidemic_curves.png` | Mean new infections per time step. | The peak height shows outbreak intensity; peak timing shows how fast the outbreak develops. | Hub has the tallest and earliest peak; strong is slower; normal remains small. |
| `plots/secondary_distribution.png` | Distribution of the number of secondary infections caused by one infected individual. | A longer tail means rare individuals infect many others. | Superspreader models produce heavier tails than the normal model. |
| `plots/secondary_distribution_normal.png` | Paper Fig. 12 style view for the no-superspreader case. | Most individuals have zero or few secondary infections. | Without superspreaders, the tail is short. |
| `plots/secondary_distribution_superspreaders.png` | Paper Fig. 13 style view comparing strong and hub superspreader cases. | The right tail shows individuals who infect many others. | Both superspreader models produce long-tailed secondary infection distributions. |
| `plots/sensitivity_lambda_attack_rate.png` | Mean final attack rate versus `lambda`. | Higher attack rate means a larger final epidemic. | Attack rate rises sharply as superspreaders are added, then approaches saturation. |
| `plots/infection_route_normal.png` | Spatial infection routes in a representative normal-model run. | Points are individuals, edges are infection events, and color indicates infection time. | Normal spread is mostly local and slower. |
| `plots/infection_route_strong.png` | Spatial infection routes in a representative strong-model run. | Red rings mark superspreaders. | Strong superspreaders create larger local bursts but still mainly transmit nearby. |
| `plots/infection_route_hub.png` | Spatial infection routes in a representative hub-model run. | Long route connections reveal wider effective contacts. | Hub superspreaders bridge distant areas and accelerate spatial spread. |

## Overall Interpretation

Across all experiments, the results support the reference paper's qualitative conclusion:
a small fraction of superspreaders can substantially increase outbreak size, outbreak
probability, and spatial speed. The important distinction is that the two superspreader
mechanisms do not behave the same way. The strong infectiousness model increases
transmission probability after contact, while the hub model changes the effective contact
geometry by allowing wider connections. Because of that, the hub model percolates at
lower density and has higher propagation speed.

## Paper Figures Not Reproduced

The reference paper's Fig. 1 and Fig. 2 are model-definition plots of the infection
probability function, so they are not included in this Member 4 experiment report.
The paper's Fig. 14 and Fig. 15 use empirical SARS Singapore 2003 data, which is also
not generated here because the raw empirical dataset is not included in this repository.
The current outputs focus on the experiment-side result families: percolation curves,
critical density, front-distance/velocity, epidemic curves, infection routes, and
secondary-infection distributions.

## Notes for Final Group Report

The experiments support the qualitative conclusion that a small fraction of
superspreaders can sharply increase epidemic size and speed. The hub model is more
dangerous than the strong-infectiousness model under the same `lambda` because it
changes the contact geometry, not only the probability of infection after contact.
