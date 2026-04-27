"""Decision-oriented ranking: WT affinity + mutant selectivity (delta scores)."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_candidate_ranking(paths: dict[str, Path], cfg: WorkflowConfig) -> None:
    t0 = time.perf_counter()
    summary_dir = paths["summary"]
    summary_dir.mkdir(parents=True, exist_ok=True)
    delta_path = paths["poses_mutant"] / "delta_scores.parquet"
    if not delta_path.is_file():
        pd.DataFrame().to_parquet(summary_dir / "candidates_ranked.parquet", index=False)
        (summary_dir / "candidates_ranking_meta.json").write_text(
            json.dumps({"note": "no_delta_scores", "skipped": True}),
            encoding="utf-8",
        )
        log_step(paths, "candidate_ranking", time.perf_counter() - t0, output_count=0)
        return

    delta = pd.read_parquet(delta_path)
    if delta.empty:
        pd.DataFrame().to_parquet(summary_dir / "candidates_ranked.parquet", index=False)
        log_step(paths, "candidate_ranking", time.perf_counter() - t0, output_count=0)
        return

    pool = pd.read_parquet(paths["filters"] / "dock_pool.parquet")
    cols = ["compound_id", "hit_smiles"]
    extra = [c for c in cols if c in pool.columns]
    out = delta.merge(pool[extra], on="compound_id", how="left")

    n = len(out)
    wt_component = (n + 1) - out["score_wt"].rank(ascending=True, method="average")
    del_component = out["delta_score"].rank(ascending=False, method="average")
    out["composite_rank_score"] = wt_component + del_component
    out = out.sort_values("composite_rank_score", ascending=False, kind="mergesort").reset_index(
        drop=True
    )
    out["rationale"] = (
        "Vina-like scores: lower WT is stronger binding. delta=score_mut-score_wt; "
        "more positive delta suggests weaker binding to mutant (selectivity). "
        "composite_rank_score=higher_is_better (rank sum). "
        "WT=" + out["score_wt"].map(lambda x: f"{float(x):.4f}")
        + "; mut=" + out["score_mut"].map(lambda x: f"{float(x):.4f}")
        + "; delta=" + out["delta_score"].map(lambda x: f"{float(x):.4f}")
    )

    out.to_parquet(summary_dir / "candidates_ranked.parquet", index=False)
    (summary_dir / "candidates_ranking_meta.json").write_text(
        json.dumps(
            {
                "n_candidates": len(out),
                "selectivity_definition": "delta_score = score_mut - score_wt",
                "composite": "sum of inverted WT rank and delta rank (higher better)",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    log_step(paths, "candidate_ranking", time.perf_counter() - t0, output_count=len(out))
