"""Tier B — structural alerts stub; wide rationale with rule_*_pass columns."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_tier_b(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    ta = pd.read_parquet(paths["filters"] / "tier_a_rationale.parquet")
    hits = pd.read_parquet(paths["pubchem"] / "hits.parquet")
    merged = hits.merge(ta, on="compound_id", how="left")
    rows = []
    for _, r in merged.iterrows():
        smi = str(r.get("hit_smiles", ""))
        rule_alert_pass = len(smi) < 500
        rule_qsar_stub_pass = True
        passed = bool(r.get("passed_tier_a", False)) and rule_alert_pass and rule_qsar_stub_pass
        prov = {"alert_set": "stub_v1", "qsar": "none"}
        rows.append(
            {
                "compound_id": r["compound_id"],
                "passed_tier_a": r.get("passed_tier_a", False),
                "rule_alert_pass": rule_alert_pass,
                "rule_qsar_stub_pass": rule_qsar_stub_pass,
                "passed_tier_ab": passed,
                "tier_b_provenance_json": json.dumps(prov),
            }
        )
    df = pd.DataFrame(rows)
    out = paths["filters"] / "tier_b_rationale.parquet"
    df.to_parquet(out, index=False)
    log_step(paths, "tier_b", time.perf_counter() - t0, input_count=len(merged), output_count=len(df))
    return out
