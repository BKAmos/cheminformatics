"""Docking merge and dual-backend wiring."""

from __future__ import annotations

from pathlib import Path

import httpx
import pandas as pd

from workflow.config_load import load_workflow_config
from workflow.contracts import DockingConfig
from workflow.docking_run import (
    merge_docking_results,
    merge_multipocket_merged_frames,
    run_configured_docking,
    write_docking_artifacts,
)
from workflow.pipeline import run_pipeline

_PUBCHEM_JSON = {
    "InformationList": {
        "Information": [
            {"CID": 702, "Similarity": 0.95, "CanonicalSMILES": "CCO"},
        ]
    }
}


def test_merge_multipocket_keeps_best_score_per_compound() -> None:
    a = pd.DataFrame([{"compound_id": "x", "score": -8.0}, {"compound_id": "y", "score": -7.0}])
    b = pd.DataFrame([{"compound_id": "x", "score": -6.0}, {"compound_id": "y", "score": -9.0}])
    m = merge_multipocket_merged_frames([a, b])
    assert set(m["compound_id"]) == {"x", "y"}
    assert float(m.loc[m["compound_id"] == "x", "score"].iloc[0]) == -8.0
    assert float(m.loc[m["compound_id"] == "y", "score"].iloc[0]) == -9.0


def test_merge_both_prefers_physics_score() -> None:
    physics = pd.DataFrame(
        [
            {"compound_id": "a", "score": -9.0, "pose_path": "/p/a.pdb", "backend": "mock"},
        ]
    )
    ml = pd.DataFrame(
        [
            {"compound_id": "a", "score": -5.0, "pose_path": "/m/a.pdb", "backend": "mock_ml"},
        ]
    )
    m = merge_docking_results(physics, ml, "both", "physics")
    assert m["score"].iloc[0] == -9.0
    assert m["pose_path"].iloc[0] == "/p/a.pdb"
    assert "score_physics" in m.columns and "score_ml" in m.columns
    assert m["backend"].iloc[0] == "mock|mock_ml"


def test_merge_both_prefers_ml_score() -> None:
    physics = pd.DataFrame(
        [
            {"compound_id": "a", "score": -9.0, "pose_path": "/p/a.pdb", "backend": "mock"},
        ]
    )
    ml = pd.DataFrame(
        [
            {"compound_id": "a", "score": -5.0, "pose_path": "/m/a.pdb", "backend": "mock_ml"},
        ]
    )
    m = merge_docking_results(physics, ml, "both", "ml")
    assert m["score"].iloc[0] == -5.0
    assert m["pose_path"].iloc[0] == "/m/a.pdb"


def test_run_configured_docking_ml_only(tmp_path: Path, repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(
        update={
            "library_csv": repo_root / "data" / "library.csv",
            "receptor_pdb": repo_root / "data" / "toy.pdb",
            "docking": DockingConfig(mode="ml_only", ml_backend="mock"),
        }
    )
    pool = pd.DataFrame([{"compound_id": "x1", "hit_smiles": "CCO"}])
    pocket = tmp_path / "pocket_spec.json"
    pocket.write_text('{"center": [0,0,0], "size": [10,10,10]}', encoding="utf-8")
    rec = repo_root / "data" / "toy.pdb"
    poses = tmp_path / "poses"
    merged, phys, ml = run_configured_docking(
        cfg,
        dock_pool=pool,
        receptor_pdb=rec,
        pocket_spec_path=pocket,
        poses_root=poses,
    )
    assert phys is None and ml is not None
    assert (ml["backend"] == "mock_ml").all()
    assert merged["score"].equals(ml["score"])
    write_docking_artifacts(poses, cfg, merged, phys, ml)
    assert (poses / "docking_scores_ml.parquet").exists()
    assert not (poses / "docking_scores_physics.parquet").exists()


def test_pipeline_tandem_writes_sidecars(tmp_path: Path, repo_root: Path, httpx_mock) -> None:
    httpx_mock.add_callback(lambda _r: httpx.Response(200, json=_PUBCHEM_JSON))

    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(
        update={
            "library_csv": repo_root / "data" / "library.csv",
            "receptor_pdb": repo_root / "data" / "toy.pdb",
            "docking": DockingConfig(mode="both", physics_backend="mock", ml_backend="mock"),
        }
    )
    run_dir = tmp_path / "tandem"
    run_pipeline(cfg, run_dir, use_langgraph=True)
    wt = run_dir / "poses" / "wt"
    assert (wt / "docking_scores.parquet").exists()
    assert (wt / "docking_scores_physics.parquet").exists()
    assert (wt / "docking_scores_ml.parquet").exists()
    df = pd.read_parquet(wt / "docking_scores.parquet")
    assert "score_physics" in df.columns and "score_ml" in df.columns
