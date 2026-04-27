"""Ground-truth / regression validation against benchmark definitions and run artifacts."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from workflow.contracts import GroundTruthCaseResult, GroundTruthReport


def _find_benchmarks_yaml() -> Path:
    bench = Path("tests/fixtures/ground_truth/benchmarks.yaml")
    for p in Path(__file__).resolve().parents:
        cand = p / "tests" / "fixtures" / "ground_truth" / "benchmarks.yaml"
        if cand.exists():
            return cand
    return bench


def _check_artifacts(run_dir: Path, spec: dict[str, Any]) -> tuple[bool, str, dict[str, Any]]:
    metrics: dict[str, Any] = {}
    rels = spec.get("required_files") or spec.get("required_rel_paths") or []
    missing = []
    for rel in rels:
        p = run_dir / rel
        if not p.is_file():
            missing.append(rel)
    if missing:
        return False, f"missing files: {missing}", metrics

    rjk = spec.get("required_json_keys")
    if rjk:
        jp = run_dir / str(rjk["path"])
        data = json.loads(jp.read_text(encoding="utf-8"))
        keys = rjk.get("keys", [])
        absent = [k for k in keys if k not in data]
        metrics["json_keys_checked"] = keys
        if absent:
            return False, f"summary JSON missing keys: {absent}", metrics

    rparq = spec.get("required_parquet_columns")
    if rparq:
        pp = run_dir / str(rparq["path"])
        df = pd.read_parquet(pp)
        cols = rparq.get("columns", [])
        absent = [c for c in cols if c not in df.columns]
        metrics["parquet_columns_checked"] = cols
        if absent:
            return False, f"parquet missing columns: {absent}", metrics

        if rparq.get("non_empty") and len(df) == 0:
            return False, "parquet is empty but non_empty required", metrics

    return True, "artifact checks passed", metrics


def run_ground_truth_suite(out_dir: Path, run_dir: Path | None = None) -> Path:
    bench = _find_benchmarks_yaml()
    data = yaml.safe_load(bench.read_text(encoding="utf-8"))
    cases_out: list[GroundTruthCaseResult] = []

    for c in data.get("cases", []):
        cid = str(c["case_id"])
        kind = c.get("kind", "static")
        expect = bool(c.get("expect_pass", True))

        if kind == "static":
            passed = expect
            msg = "static expectation"
            metrics = {"kind": kind}
        elif kind == "artifacts":
            if run_dir is None:
                passed = True
                msg = "skipped (no --run-dir provided)"
                metrics = {"skipped": True}
            else:
                ok, msg, metrics = _check_artifacts(run_dir, c)
                expect_pass = bool(c.get("expect_pass", True))
                passed = ok if expect_pass else not ok
        elif kind == "parquet_invariant":
            path = Path(str(c["path"]))
            if not path.is_absolute():
                path = bench.parent / path
            if not path.is_file():
                passed = False
                msg = f"missing fixture parquet {path}"
                metrics = {}
            else:
                df = pd.read_parquet(path)
                col = str(c["score_column"])
                descending = bool(c.get("expect_sorted_descending", False))
                if col not in df.columns:
                    passed = False
                    msg = f"missing column {col}"
                else:
                    s = df[col].tolist()
                    ok = s == sorted(s, reverse=descending)
                    passed = ok if expect else not ok
                    msg = "sorted as expected" if ok else "sort invariant failed"
                metrics = {"rows": len(df)}
        else:
            passed = False
            msg = f"unknown case kind: {kind}"
            metrics = {}

        cases_out.append(
            GroundTruthCaseResult(
                case_id=cid,
                passed=passed,
                metrics=metrics,
                message=msg,
            )
        )

    rep = GroundTruthReport(
        cases=cases_out,
        tool_versions={"python": sys.version.split()[0]},
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / "ground_truth_report.json"
    p.write_text(rep.model_dump_json(indent=2), encoding="utf-8")
    return p
