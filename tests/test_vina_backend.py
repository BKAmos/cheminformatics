from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from workflow.backends.vina import VinaDockingBackend, parse_vina_score


def test_parse_vina_score_from_standard_log_table() -> None:
    log = """
-----+------------+----------+----------
mode |   affinity | dist from best mode
-----+------------+----------+----------
   1      -7.4      0.000      0.000
   2      -6.9      1.512      2.015
"""
    assert parse_vina_score(log) == -7.4


def test_parse_vina_score_raises_on_missing_table() -> None:
    bad = "no affinity rows here"
    with pytest.raises(ValueError):
        parse_vina_score(bad)


def test_vina_backend_subprocess_integration_mock(monkeypatch, tmp_path: Path) -> None:
    receptor = tmp_path / "rec.pdb"
    receptor.write_text("ATOM\n", encoding="utf-8")
    pocket = tmp_path / "pocket.json"
    pocket.write_text('{"center":[1.0,2.0,3.0],"size":[10.0,10.0,10.0]}', encoding="utf-8")
    poses = tmp_path / "poses"
    pool = pd.DataFrame([{"compound_id": "c1", "hit_smiles": "CCO"}])

    def fake_which(name: str):
        return "vina" if name == "vina" else None

    def fake_run_cmd(argv, **kwargs):
        if argv[0] == "obabel":
            out = Path(argv[argv.index("-O") + 1])
            out.write_text("PDBQT\n", encoding="utf-8")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if argv[0] == "vina":
            assert "--seed" in argv and "--exhaustiveness" in argv
            out_pose = Path(argv[argv.index("--out") + 1])
            assert "--log" not in argv
            out_pose.write_text("POSE\n", encoding="utf-8")
            return SimpleNamespace(
                returncode=0,
                stdout="   1      -7.8      0.000      0.000\n",
                stderr="",
            )
        return SimpleNamespace(returncode=1, stdout="", stderr="unknown")

    monkeypatch.setattr("workflow.backends.vina.shutil.which", fake_which)
    monkeypatch.setattr("workflow.backends.vina.obabel_argv0", lambda: ["obabel"])
    monkeypatch.setattr("workflow.backends.vina.run_cmd", fake_run_cmd)

    out = VinaDockingBackend().dock_batch(
        receptor_pdb=receptor,
        pocket_spec_path=pocket,
        dock_pool=pool,
        poses_dir=poses,
        max_parallel=1,
    )
    assert len(out) == 1
    assert out.iloc[0]["backend"] == "vina"
    assert float(out.iloc[0]["score"]) == -7.8
    assert (poses / "logs" / "c1.vina.log").exists()


def test_vina_uses_precomputed_ligand_pdbqt_path(monkeypatch, tmp_path: Path) -> None:
    receptor = tmp_path / "rec.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    ligand = tmp_path / "lig.pdbqt"
    ligand.write_text("LIGAND\n", encoding="utf-8")
    pocket = tmp_path / "pocket.json"
    pocket.write_text('{"center":[0,0,0],"size":[8,8,8]}', encoding="utf-8")
    poses = tmp_path / "poses"
    pool = pd.DataFrame(
        [{"compound_id": "c1", "hit_smiles": "CCO", "ligand_pdbqt_path": str(ligand)}]
    )
    calls = []

    def fake_which(name: str):
        return "vina" if name == "vina" else None

    def fake_run_cmd(argv, **kwargs):
        calls.append(argv)
        if argv[0] == "vina":
            out_pose = Path(argv[argv.index("--out") + 1])
            out_pose.write_text("POSE\n", encoding="utf-8")
            return SimpleNamespace(
                returncode=0,
                stdout="   1      -6.2      0.000      0.000\n",
                stderr="",
            )
        if argv[0] == "obabel":
            return SimpleNamespace(returncode=1, stdout="", stderr="should not be called")
        return SimpleNamespace(returncode=1, stdout="", stderr="unknown")

    monkeypatch.setattr("workflow.backends.vina.shutil.which", fake_which)
    monkeypatch.setattr("workflow.backends.vina.run_cmd", fake_run_cmd)

    out = VinaDockingBackend().dock_batch(
        receptor_pdb=receptor,
        pocket_spec_path=pocket,
        dock_pool=pool,
        poses_dir=poses,
        max_parallel=1,
    )
    assert len(out) == 1
    assert all(cmd[0] != "obabel" for cmd in calls)

