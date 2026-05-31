# Member 4 Report: Experiments and Sensitivity Analysis

## Scope

This report covers Member 4 responsibilities: percolation probability experiments,
propagation speed, epidemic curves, secondary-infection distribution, and sensitivity
checks over population size, superspreader fraction, density, and random seeds. The
figures are organized to mirror the reference paper: infection probability curves,
percolation curves, critical density, propagation front distance, velocity, infection
routes, epidemic curves, and secondary-infection distributions.

## Reproducibility

- Simulator: `MonteCarloSIRSimulator`
- Random seeds per setting: `10`
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

## Baseline Results

| Model | Mean attack rate | Percolation probability | Mean duration | Mean propagation speed | Mean secondary infections |
| --- | ---: | ---: | ---: | ---: | ---: |
| normal | 0.075 | 0.000 | 9.000 | 0.086 | 0.557 |
| strong | 0.586 | 0.600 | 16.600 | 0.254 | 0.688 |
| hub | 0.690 | 0.700 | 11.600 | 0.493 | 0.699 |

The hub model produces the fastest baseline spread. The strong-infectiousness model
also increases final epidemic size relative to the normal model, but its propagation
speed is lower because superspreaders still infect within the shorter normal range.

The baseline results show a clear ordering of epidemic severity:
`hub > strong > normal`. In the normal model, only about `7.5%` of the population is
infected on average and none of the runs percolate. In the strong model, the attack rate
increases to about `58.6%`, and `60%` of runs reach the top of the system. In the hub
model, the attack rate rises to about `69.0%`, percolation probability reaches `70%`,
and propagation speed is almost twice the strong-model speed. This supports the main
claim that superspreaders do not only increase the final outbreak size; they also change
how quickly the epidemic moves through space.

## Percolation and Density Sweep

| Model | lambda | Density | Percolation probability | Mean attack rate | Mean propagation speed |
| --- | ---: | ---: | ---: | ---: | ---: |
| hub | 0.20 | 1.50 | 0.00 | 0.02 | 0.25 |
| hub | 0.20 | 2.00 | 0.00 | 0.03 | 0.20 |
| hub | 0.20 | 2.50 | 0.00 | 0.05 | 0.28 |
| hub | 0.20 | 3.00 | 0.60 | 0.55 | 0.26 |
| hub | 0.20 | 3.50 | 0.60 | 0.62 | 0.39 |
| hub | 0.20 | 4.00 | 0.70 | 0.67 | 0.41 |
| hub | 0.20 | 4.77 | 0.70 | 0.69 | 0.49 |
| hub | 0.20 | 5.50 | 0.80 | 0.79 | 0.63 |
| hub | 0.20 | 6.50 | 0.80 | 0.80 | 0.64 |
| hub | 0.20 | 8.00 | 0.80 | 0.80 | 0.65 |
| normal | 0.00 | 1.50 | 0.00 | 0.00 | 0.01 |
| normal | 0.00 | 2.00 | 0.00 | 0.00 | 0.01 |
| normal | 0.00 | 2.50 | 0.00 | 0.00 | 0.03 |
| normal | 0.00 | 3.00 | 0.00 | 0.00 | 0.02 |
| normal | 0.00 | 3.50 | 0.00 | 0.01 | 0.03 |
| normal | 0.00 | 4.00 | 0.00 | 0.01 | 0.04 |
| normal | 0.00 | 4.77 | 0.00 | 0.07 | 0.09 |
| normal | 0.00 | 5.50 | 0.20 | 0.15 | 0.19 |
| normal | 0.00 | 6.50 | 0.80 | 0.73 | 0.24 |
| normal | 0.00 | 8.00 | 0.90 | 0.87 | 0.30 |
| strong | 0.20 | 1.50 | 0.00 | 0.01 | 0.15 |
| strong | 0.20 | 2.00 | 0.00 | 0.02 | 0.13 |
| strong | 0.20 | 2.50 | 0.00 | 0.03 | 0.16 |
| strong | 0.20 | 3.00 | 0.00 | 0.08 | 0.21 |
| strong | 0.20 | 3.50 | 0.00 | 0.23 | 0.20 |
| strong | 0.20 | 4.00 | 0.50 | 0.45 | 0.21 |
| strong | 0.20 | 4.77 | 0.60 | 0.59 | 0.25 |
| strong | 0.20 | 5.50 | 0.80 | 0.79 | 0.40 |
| strong | 0.20 | 6.50 | 0.80 | 0.80 | 0.43 |
| strong | 0.20 | 8.00 | 0.80 | 0.80 | 0.48 |

Percolation probability increases with density because each infectious individual has
more neighbors inside the infection radius. Hub superspreaders reach the high
percolation regime earlier than the strong model because their effective infection range
is larger.

At the baseline superspreader fraction `lambda = 0.20`, the hub model starts
percolating around density `3.00`, while the strong model does not reach a `0.50`
percolation probability until density `4.00`. The normal model needs a much denser
population, reaching high percolation only around density `6.50`. This means the hub
mechanism lowers the density barrier for large outbreaks. In practical terms, a hub-like
superspreader can connect distant local clusters, so the disease crosses the system
even when ordinary local transmission would still die out.

## Critical Density

| Model | lambda | Critical density |
| --- | ---: | ---: |
| strong | 0.00 | 6.00 |
| strong | 0.05 | 6.00 |
| strong | 0.10 | 4.51 |
| strong | 0.20 | 4.00 |
| strong | 0.30 | 3.75 |
| strong | 0.40 | 3.00 |
| hub | 0.00 | 6.00 |
| hub | 0.05 | 5.50 |
| hub | 0.10 | 4.19 |
| hub | 0.20 | 2.92 |
| hub | 0.30 | 2.88 |
| hub | 0.40 | 2.17 |

The critical-density plot summarizes the percolation curves into one threshold value.
As `lambda` increases, both models need less population density to percolate. The hub
model drops much faster: at `lambda = 0.20`, its critical density is about `2.92`, while
the strong model needs about `4.00`. Lower critical density means the outbreak is more
robust under sparse conditions, so this result again shows that hub superspreaders are
more dangerous than merely stronger infectious individuals.

## Sensitivity Findings

- Increasing `lambda` raises attack rate and speed for both superspreader models.
- At baseline `N` and density, the strongest hub speed in the sweep occurs at
  `lambda = 0.40` with mean speed
  `0.754`.
- At baseline `N` and density, the strongest strong-model speed in the sweep occurs at
  `lambda = 0.40` with mean speed
  `0.471`.
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
| `plots/infection_probability_strong.png` | Infection probability as a function of distance for normal and strong infectiousness sources. | The strong source has higher transmission probability inside the normal infection radius. | Strong superspreaders increase risk per contact, but they do not extend the spatial range. |
| `plots/infection_probability_hub.png` | Infection probability as a function of distance for normal and hub sources. | The hub curve decays over a larger effective range. | Hub superspreaders can connect individuals that are too far for normal transmission. |
| `plots/percolation_probability.png` | Percolation probability versus density for normal, strong, and hub models at baseline `lambda`. | Curves farther left percolate at lower density. | Hub percolates first, strong second, normal last. |
| `plots/percolation_probability_strong.png` | Strong-model percolation probability for several `lambda` values. | Increasing `lambda` shifts the curve left and upward. | More strong superspreaders make system-wide spread possible at lower density. |
| `plots/percolation_probability_hub.png` | Hub-model percolation probability for several `lambda` values. | The left shift is stronger than in the strong model. | Hub superspreaders reduce the percolation threshold more sharply. |
| `plots/critical_density.png` | Approximate density where percolation probability reaches `0.5`. | Lower values mean easier system-wide spread. | Critical density decreases with `lambda`, and the hub threshold is lower than the strong threshold. |
| `plots/front_distance_strong.png` | Mean infection-front distance over time for the strong model. | Steeper growth means faster spatial propagation; a plateau means spread has stopped. | Larger `lambda` makes the front travel farther and faster. |
| `plots/front_distance_hub.png` | Mean infection-front distance over time for the hub model. | Compare curve height and slope across `lambda`. | Hub spreading reaches far distances earlier than the strong model. |
| `plots/velocity_vs_lambda.png` | Propagation speed versus superspreader fraction. | Higher speed means the epidemic wave crosses space faster. | Hub speed is consistently higher than strong speed, especially at large `lambda`. |
| `plots/epidemic_curves.png` | Mean new infections per time step. | The peak height shows outbreak intensity; peak timing shows how fast the outbreak develops. | Hub has the tallest and earliest peak; strong is slower; normal remains small. |
| `plots/secondary_distribution.png` | Distribution of the number of secondary infections caused by one infected individual. | A longer tail means rare individuals infect many others. | Superspreader models produce heavier tails than the normal model. |
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

## Notes for Final Group Report

The experiments support the qualitative conclusion that a small fraction of
superspreaders can sharply increase epidemic size and speed. The hub model is more
dangerous than the strong-infectiousness model under the same `lambda` because it
changes the contact geometry, not only the probability of infection after contact.
