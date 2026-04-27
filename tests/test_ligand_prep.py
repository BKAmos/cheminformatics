from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from workflow.config_load import load_workflow_config
from workflow.run_layout import ensure_run_dirs
from workflow.steps.ligand_prep import run_ligand_prep


def test_ligand_prep_enriches_dock_pool_with_pdbqt_paths(
    monkeypatch, tmp_path: Path, repo_root: Path
) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    paths = ensure_run_dirs(tmp_path)
    pool = pd.DataFrame(
        [
            {
                "compound_id": "c1",
                "hit_smiles": "CCO",
                "tanimoto_similarity": 0.9,
                "passed_tier_ab": True,
            },
            {
                "compound_id": "c2",
                "hit_smiles": "CCC",
                "tanimoto_similarity": 0.8,
                "passed_tier_ab": True,
            },
        ]
    )
    pool.to_parquet(paths["filters"] / "dock_pool.parquet", index=False)

    def fake_run_cmd(argv, **kwargs):
        out = Path(argv[argv.index("-O") + 1])
        out.write_text("PDBQT\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("workflow.steps.ligand_prep.obabel_argv0", lambda: ["obabel"])
    monkeypatch.setattr("workflow.steps.ligand_prep.run_cmd", fake_run_cmd)

    run_ligand_prep(paths, cfg)
    out_pool = pd.read_parquet(paths["filters"] / "dock_pool.parquet")
    assert "ligand_pdbqt_path" in out_pool.columns
    assert out_pool["ligand_pdbqt_path"].notna().all()

