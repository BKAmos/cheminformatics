"""Gate + dock pool cap (max_compounds xor top_n_by_similarity)."""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_dock_pool_cap(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    tb = pd.read_parquet(paths["filters"] / "tier_b_rationale.parquet")
    hits = pd.read_parquet(paths["pubchem"] / "hits.parquet")
    m = hits.merge(tb, on="compound_id", how="inner")
    passed = m[m["passed_tier_ab"] == True].copy()  # noqa: E712

    if cfg.dock_pool.mode == "top_n_similarity":
        n = cfg.dock_pool.top_n_by_2d_score or 100
        passed = passed.sort_values("tanimoto_similarity", ascending=False).head(n)
    else:
        cap = cfg.dock_pool.max_compounds_to_dock or 1000
        passed = passed.head(cap)

    passed["cap_mode"] = cfg.dock_pool.mode
    out = paths["filters"] / "dock_pool.parquet"
    passed.to_parquet(out, index=False)
    log_step(
        paths,
        "dock_pool_cap",
        time.perf_counter() - t0,
        input_count=len(m),
        output_count=len(passed),
    )
    return out
