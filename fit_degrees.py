"""Run degree-tail inference for the supplied BRCA1 association network."""

from pathlib import Path
import sys

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from brca1_analysis import load_graph
from powerlaw_analysis import run_analysis

BOOTSTRAPS = 1000
SEED = 2026


def main():
    graph = load_graph()
    degrees = np.array([degree for _, degree in graph.degree()], dtype=int)
    run_analysis(degrees, BASE_DIR / "output", "BRCA1", graph.graph["input_sha256"],
                 bootstraps=BOOTSTRAPS, seed=SEED)


if __name__ == "__main__":
    main()
