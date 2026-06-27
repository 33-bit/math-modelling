"""Backward-compatible wrapper for :mod:`epidemic_model.experiment_pipeline`."""

import importlib
import sys

from epidemic_model.experiment_pipeline import *  # noqa: F401,F403

_SUBMODULES = (
    "cli",
    "files",
    "metrics",
    "plotting",
    "records",
    "scenarios",
    "settings",
    "tables",
)

for _name in _SUBMODULES:
    sys.modules[f"{__name__}.{_name}"] = importlib.import_module(
        f"epidemic_model.experiment_pipeline.{_name}"
    )

