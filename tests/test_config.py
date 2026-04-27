from pathlib import Path

from workflow.config_load import load_workflow_config


def test_load_example_config(repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    assert cfg.library_csv.exists()
    assert cfg.receptor_pdb.exists()
    assert cfg.resources.max_parallel_docks >= 1
    assert cfg.docking.mode == "physics_only"
    assert cfg.docking.physics_backend == "mock"
