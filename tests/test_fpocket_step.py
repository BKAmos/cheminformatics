"""fpocket step: manual box table + subprocess stub."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from workflow.config_load import load_workflow_config
from workflow.contracts import FpocketConfig, WorkflowConfig
from workflow.run_layout import ensure_run_dirs
from workflow.steps.fpocket_parse import run_fpocket


def test_fpocket_manual_writes_ranked_parquet(tmp_path: Path, repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    paths = ensure_run_dirs(tmp_path)
    shutil.copy(repo_root / "data" / "toy.pdb", paths["structure"] / "receptor_prepared.pdb")
    run_fpocket(paths, cfg)
    ranked = paths["structure"] / "pockets_ranked.parquet"
    assert ranked.is_file()
    df = pd.read_parquet(ranked)
    assert len(df) == 1
    assert df.iloc[0]["spec_file"] == "pocket_spec.json"


def test_fpocket_binary_parses_pockets(tmp_path: Path, repo_root: Path, monkeypatch) -> None:
    paths = ensure_run_dirs(tmp_path)
    shutil.copy(repo_root / "data" / "toy.pdb", paths["structure"] / "receptor_prepared.pdb")
    cfg = WorkflowConfig(
        library_csv=repo_root / "data" / "library.csv",
        receptor_pdb=repo_root / "data" / "toy.pdb",
        fpocket=FpocketConfig(
            mode="fpocket",
            use_stub_box=False,
            top_k_pockets=2,
        ),
    )

    def fake_run_fpocket_subprocess(
        *, receptor_pdb: Path, fpocket_raw: Path, executable: str | None = None, timeout_s: float = 600.0
    ) -> None:
        out = fpocket_raw / "toy_out" / "pockets"
        out.mkdir(parents=True, exist_ok=True)
        atm = out / "pocket1_atm.pdb"
        atm.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "ATOM      2  CA  ALA A   1      10.000   0.000   0.000  1.00 20.00           C\n",
            encoding="utf-8",
        )
        (out / "pocket1_info.txt").write_text("Score : 12.5\n", encoding="utf-8")

    monkeypatch.setattr(
        "workflow.steps.fpocket_parse.run_fpocket_subprocess",
        fake_run_fpocket_subprocess,
    )
    monkeypatch.setattr(
        "workflow.steps.fpocket_parse.pick_output_dir",
        lambda raw, stem: raw / "toy_out",
    )

    run_fpocket(paths, cfg)
    specs = paths["structure"] / "pocket_specs"
    assert specs.is_dir()
    assert (paths["structure"] / "pocket_spec.json").is_file()
    df = pd.read_parquet(paths["structure"] / "pockets_ranked.parquet")
    assert len(df) >= 1
    assert (specs / f"pocket_{int(df.iloc[0]['pocket_id'])}.json").is_file()


def test_fpocket_utils_parse_score() -> None:
    from workflow.fpocket_utils import parse_pocket_info_score

    assert parse_pocket_info_score("Score : 3.2\n") == 3.2
