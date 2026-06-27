"""Configuration constants for the Monte Carlo SIR simulator."""

SUSCEPTIBLE = 0
INFECTED = 1
RECOVERED = 2

STATE_LABELS = {
    SUSCEPTIBLE: "S",
    INFECTED: "I",
    RECOVERED: "R",
}

DEFAULT_R0 = 1.0
DEFAULT_W0 = 1.0
DEFAULT_GAMMA = 1.0
DEFAULT_MAX_STEPS = 200
DEFAULT_N = 477

VALID_MODELS = {"normal", "strong", "hub"}


def resolve_L(r0: float, L: float | None = None) -> float:
    """Return the simulation box side length."""
    if L is None:
        return 10.0 * r0
    return float(L)
