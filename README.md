# Monte Carlo SIR Simulator with Superspreaders

A simple spatial Monte Carlo SIR simulator with normal, strong superspreader, and hub superspreader models.

## Features

- Generates fixed population positions in an `L × L` square.
- Uses periodic boundary distance.
- Simulates SIR infection and recovery over discrete timesteps.
- Supports normal, strong superspreader, and hub superspreader models.
- Supports reproducible runs with random seeds.
- Returns infection history, recovery history, secondary infection counts, and epidemic curve data.

## Files

- `config.py`: default parameters and SIR state constants.
- `simulator.py`: main simulation code.
- `demo_run.py`: example script showing how to run the simulator.
- `requirements.txt`: required Python package.

## Installation

```bash
pip install -r requirements.txt
```

## Run Demo

```bash
python demo_run.py
```

The demo runs three cases:

- `normal`
- `strong`
- `hub`

It prints the total infected count, duration, new infections per step, and top secondary infection counts.

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

Use one of the following model names:

- `model="normal"`
- `model="strong"`
- `model="hub"`

For `strong` and `hub`, use `lambda_ss` to set the fraction of superspreaders.

Example:

```python
sim = MonteCarloSIRSimulator(
    N=500,
    model="strong",
    lambda_ss=0.4,
    seed=1,
)
```

## SimulationResult

`sim.run()` returns a `SimulationResult` object with:

- `positions`: positions of all individuals.
- `states`: final S/I/R state of each individual.
- `is_superspreader`: superspreader flags.
- `infected_time`: infection time of each individual.
- `recovered_time`: recovery time of each individual.
- `infection_source`: direct source of infection for each individual.
- `secondary_counts`: number of direct infections caused by each individual.
- `new_infections_per_step`: number of new infections at each timestep.
- `total_infected`: total number of individuals ever infected.
- `duration`: number of timesteps executed.
