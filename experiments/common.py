"""Shared utilities for the experiments."""
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

RESULTS = os.path.join(ROOT, "results")


def ibexopt_path():
    path = os.environ.get("IBEXOPT") or shutil.which("ibexopt")
    if not path:
        sys.exit("IbexOpt not found: install Ibex with the optim plugin and set IBEXOPT=/path/to/ibexopt")
    return path
