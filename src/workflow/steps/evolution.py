"""Evolution stub — writes minimal lineage when not skipped."""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_evolution(paths: dict[str, Path], cfg: WorkflowConfig) -> None:
    t0 = time.perf_counter()
    if cfg.skip_evolution:
        pd.DataFrame(
            [{"generation": 0, "note": "evolution_skipped", "fitness": None, "smiles": None}]
        ).to_parquet(paths["evolution"] / "generations.parquet", index=False)
        log_step(paths, "evolution", time.perf_counter() - t0, output_count=0)
        return
    wt = pd.read_parquet(paths["poses_wt"] / "docking_scores.parquet")
    if wt.empty:
        pd.DataFrame(
            [{"generation": 0, "note": "no_dock_scores", "fitness": None, "smiles": None}]
        ).to_parquet(paths["evolution"] / "generations.parquet", index=False)
        log_step(paths, "evolution", time.perf_counter() - t0, output_count=0)
        return
    k = min(3, len(wt))
    best = wt.nsmallest(k, "score").copy()
    best["generation"] = 1
    best["fitness"] = best["score"]
    best["smiles"] = "placeholder"
    best.to_parquet(paths["evolution"] / "generations.parquet", index=False)
    log_step(paths, "evolution", time.perf_counter() - t0, output_count=len(best))
