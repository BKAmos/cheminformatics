import json
from pathlib import Path

from workflow.config_load import load_workflow_config
from workflow.pipeline import run_pipeline


def test_pipeline_ci_langgraph(tmp_path: Path, repo_root: Path, mock_pubchem) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(
        update={
            "library_csv": repo_root / "data" / "library.csv",
            "receptor_pdb": repo_root / "data" / "toy.pdb",
        }
    )
    run_dir = tmp_path / "run1"
    run_pipeline(cfg, run_dir, use_langgraph=True)
    assert (run_dir / "summary" / "summary.json").exists()
    assert (run_dir / "poses" / "wt" / "docking_scores.parquet").exists()
    assert (run_dir / "summary" / "candidates_ranked.parquet").exists()
    summary = json.loads((run_dir / "summary" / "summary.json").read_text(encoding="utf-8"))
    assert "n_ranked_candidates" in summary


def test_dry_run(tmp_path: Path, repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(update={"dry_run": True})
    run_pipeline(cfg, tmp_path / "dry", use_langgraph=True)
    s = (tmp_path / "dry" / "summary" / "summary.json").read_text()
    assert "dry_run" in s
