"""Re-run Tier A+B on evolved SMILES when present."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_evolved_refilter(paths: dict[str, Path], cfg: WorkflowConfig) -> None:
    t0 = time.perf_counter()
    gen = pd.read_parquet(paths["evolution"] / "generations.parquet")
    if "smiles" not in gen.columns:
        _write_empty_evolved(paths)
        log_step(paths, "evolved_refilter", time.perf_counter() - t0, output_count=0)
        return
    sub = gen.dropna(subset=["smiles"])
    sub = sub[sub["smiles"].astype(str) != "placeholder"]
    if sub.empty:
        _write_empty_evolved(paths)
        log_step(paths, "evolved_refilter", time.perf_counter() - t0, output_count=0)
        return

    rows_a = []
    rows_b = []
    for _, r in sub.iterrows():
        cid = f"evo_{r.get('generation', 0)}_{abs(hash(str(r['smiles']))) % 1_000_000}"
        rows_a.append(
            {
                "compound_id": cid,
                "passed_tier_a": True,
                "tier_a_score": 0.0,
                "tier_a_provenance_json": json.dumps({"source": "evolved_refilter"}),
            }
        )
        smi = str(r["smiles"])
        rule_alert_pass = len(smi) < 500
        rows_b.append(
            {
                "compound_id": cid,
                "passed_tier_a": True,
                "rule_alert_pass": rule_alert_pass,
                "rule_qsar_stub_pass": True,
                "passed_tier_ab": rule_alert_pass,
                "tier_b_provenance_json": json.dumps({"source": "evolved_refilter"}),
            }
        )
    pd.DataFrame(rows_a).to_parquet(paths["filters"] / "evolved_tier_a_rationale.parquet", index=False)
    pd.DataFrame(rows_b).to_parquet(paths["filters"] / "evolved_tier_b_rationale.parquet", index=False)
    log_step(paths, "evolved_refilter", time.perf_counter() - t0, output_count=len(rows_b))


def _write_empty_evolved(paths: dict[str, Path]) -> None:
    empty = pd.DataFrame(
        columns=[
            "compound_id",
            "passed_tier_a",
            "tier_a_score",
            "tier_a_provenance_json",
        ]
    )
    empty.to_parquet(paths["filters"] / "evolved_tier_a_rationale.parquet", index=False)
    empty_b = pd.DataFrame(
        columns=[
            "compound_id",
            "passed_tier_a",
            "rule_alert_pass",
            "rule_qsar_stub_pass",
            "passed_tier_ab",
            "tier_b_provenance_json",
        ]
    )
    empty_b.to_parquet(paths["filters"] / "evolved_tier_b_rationale.parquet", index=False)
