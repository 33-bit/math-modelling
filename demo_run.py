"""Small demo run for the Monte Carlo SIR simulator."""

from simulator import MonteCarloSIRSimulator


def main() -> None:
    configs = [
        ("normal", 0.0),
        ("strong", 0.2),
        ("hub", 0.2),
    ]

    for model, lambda_ss in configs:
        sim = MonteCarloSIRSimulator(
            N=477,
            model=model,
            lambda_ss=lambda_ss,
            seed=42,
        )
        result = sim.run()

        print("=" * 60)
        print(f"Model: {model}")
        print(f"Lambda: {lambda_ss}")
        print(f"Total infected: {result.total_infected}")
        print(f"Duration: {result.duration}")
        print(f"First 10 new infections per step: {result.new_infections_per_step[:10]}")
        top_secondary_counts = [int(count) for count in sorted(result.secondary_counts, reverse=True)[:10]]
        print(f"Top 10 secondary counts: {top_secondary_counts}")


if __name__ == "__main__":
    main()
