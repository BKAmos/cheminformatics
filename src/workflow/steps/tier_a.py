"""Tier A — data-driven stub with provenance columns (ChEMBL/ECOTOX placeholders)."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_tier_a(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    hits = pd.read_parquet(paths["pubchem"] / "hits.parquet")
    rows = []
    for _, r in hits.iterrows():
        cid = r.get("compound_id", "")
        passed = True
        prov = {
            "source": "tier_a_stub",
            "version": "0.1",
            "query": f"stub_lookup/{cid}",
            "threshold": 0.0,
            "payload_id": str(cid),
        }
        rows.append(
            {
                "compound_id": r["compound_id"],
                "passed_tier_a": passed,
                "tier_a_score": 0.0,
                "tier_a_provenance_json": json.dumps(prov),
            }
        )
    df = pd.DataFrame(rows)
    out = paths["filters"] / "tier_a_rationale.parquet"
    chunk = cfg.resources.tier_batch_size
    if len(df) > chunk:
        parts = [df.iloc[i : i + chunk] for i in range(0, len(df), chunk)]
        pd.concat(parts, ignore_index=True).to_parquet(out, index=False)
    else:
        df.to_parquet(out, index=False)
    log_step(paths, "tier_a", time.perf_counter() - t0, input_count=len(hits), output_count=len(df))
    return out
