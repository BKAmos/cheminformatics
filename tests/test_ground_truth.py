import json
from pathlib import Path

from workflow.validation.ground_truth import run_ground_truth_suite


def test_ground_truth_report(tmp_path: Path) -> None:
    p = run_ground_truth_suite(tmp_path)
    assert p.exists()
    text = p.read_text(encoding="utf-8")
    assert "toy_1" in text
    assert "passed" in text


def test_ground_truth_with_run_dir(tmp_path: Path, repo_root: Path, mock_pubchem) -> None:
    from workflow.config_load import load_workflow_config
    from workflow.pipeline import run_pipeline
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(
        update={
            "library_csv": repo_root / "data" / "library.csv",
            "receptor_pdb": repo_root / "data" / "toy.pdb",
        }
    )
    run_dir = tmp_path / "e2e"
    run_pipeline(cfg, run_dir, use_langgraph=True)
    rep = run_ground_truth_suite(tmp_path / "val", run_dir=run_dir)
    data = json.loads(rep.read_text(encoding="utf-8"))
    cases = {c["case_id"]: c for c in data["cases"]}
    assert cases["ci_run_artifacts"]["passed"] is True
    assert cases["synthetic_should_fail"]["passed"] is True
