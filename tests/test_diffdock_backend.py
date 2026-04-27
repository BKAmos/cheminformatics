"""DiffDock backend placeholder (and optional external hook) contract tests."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from workflow.backends.diffdock import DiffDockBackend


def test_diffdock_backend_writes_scores(tmp_path: Path) -> None:
    pdb_line = "ATOM      1  N   ALA A   1       0.0   0.0   0.0\n"
    (tmp_path / "rec.pdb").write_text(pdb_line, encoding="utf-8")
    pocket = tmp_path / "pocket_spec.json"
    pocket.write_text(
        json.dumps(
            {
                "center": [0.0, 0.0, 0.0],
                "size": [10.0, 10.0, 10.0],
                "source": "manual",
            }
        ),
        encoding="utf-8",
    )
    pool = pd.DataFrame(
        [
            {"compound_id": "a", "ligand_smiles": "CCO"},
            {"compound_id": "b", "ligand_smiles": "CC"},
        ]
    )
    poses = tmp_path / "poses"
    b = DiffDockBackend()
    out = b.dock_batch(
        receptor_pdb=tmp_path / "rec.pdb",
        pocket_spec_path=pocket,
        dock_pool=pool,
        poses_dir=poses,
        max_parallel=1,
    )
    assert set(out["compound_id"]) == {"a", "b"}
    assert (out["backend"] == "diffdock").all()
    for _, row in out.iterrows():
        assert Path(str(row["pose_path"])).is_file()
        assert float(row["score"]) < 0
