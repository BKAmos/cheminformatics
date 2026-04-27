"""AutoDock Vina backend (subprocess) with log parsing."""

from __future__ import annotations

import json
import re
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from workflow.subprocess_runner import run_cmd

_VINA_SCORE_RE = re.compile(r"^\s*\d+\s+(-?\d+(?:\.\d+)?)\s+")


def parse_vina_score(text: str) -> float:
    """Parse best affinity from Vina log/stdout table."""
    for line in text.splitlines():
        m = _VINA_SCORE_RE.match(line)
        if m:
            return float(m.group(1))
    raise ValueError("Could not parse Vina affinity from log output")


class VinaDockingBackend:
    name = "vina"
    family = "physics"

    def dock_batch(
        self,
        *,
        receptor_pdb: Path,
        pocket_spec_path: Path,
        dock_pool: pd.DataFrame,
        poses_dir: Path,
        max_parallel: int = 1,
        vina_seed: int = 42,
        vina_exhaustiveness: int = 8,
    ) -> pd.DataFrame:
        vina_bin = shutil.which("vina")
        if not vina_bin:
            raise RuntimeError(
                "Vina binary not found on PATH. Install AutoDock Vina or switch "
                "docking.physics_backend to 'mock' in config."
            )
        poses_dir.mkdir(parents=True, exist_ok=True)
        logs_dir = poses_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir = poses_dir / "_tmp_vina"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        spec = json.loads(pocket_spec_path.read_text(encoding="utf-8"))
        center = spec.get("center", [0.0, 0.0, 0.0])
        size = spec.get("size", [20.0, 20.0, 20.0])

        receptor_pdbqt = self._prepare_receptor_pdbqt(receptor_pdb, tmp_dir)
        workers = max(1, int(max_parallel))
        rows: list[dict[str, object]] = []
        pool = dock_pool.sort_values("compound_id", kind="mergesort")

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = []
            for row in pool.itertuples():
                futs.append(
                    ex.submit(
                        self._dock_one,
                        vina_bin=vina_bin,
                        receptor_pdbqt=receptor_pdbqt,
                        poses_dir=poses_dir,
                        logs_dir=logs_dir,
                        tmp_dir=tmp_dir,
                        row=row,
                        center=center,
                        size=size,
                        vina_seed=vina_seed,
                        exhaustiveness=vina_exhaustiveness,
                    )
                )
            for f in as_completed(futs):
                rows.append(f.result())
        out = pd.DataFrame(rows)
        return out.sort_values("compound_id", kind="mergesort").reset_index(drop=True)

    def _prepare_receptor_pdbqt(self, receptor_pdb: Path, tmp_dir: Path) -> Path:
        if receptor_pdb.suffix.lower() == ".pdbqt":
            return receptor_pdb
        obabel = shutil.which("obabel")
        if not obabel:
            raise RuntimeError(
                "Receptor is not PDBQT and Open Babel (`obabel`) is not available "
                "to convert receptor to PDBQT."
            )
        out = tmp_dir / f"{receptor_pdb.stem}.pdbqt"
        cp = run_cmd([obabel, "-ipdb", str(receptor_pdb), "-opdbqt", "-O", str(out)], timeout_s=120.0)
        if cp.returncode != 0 or not out.exists():
            raise RuntimeError(f"Failed receptor PDBQT conversion: {cp.stderr}")
        return out

    def _prepare_ligand_pdbqt(self, row, tmp_dir: Path) -> Path:
        prebuilt = getattr(row, "ligand_pdbqt_path", None)
        if prebuilt:
            p = Path(str(prebuilt))
            if p.exists():
                return p
        smiles = getattr(row, "hit_smiles", None)
        if not smiles:
            raise RuntimeError(
                f"Dock row {getattr(row, 'compound_id', 'unknown')} lacks hit_smiles "
                "or ligand_pdbqt_path"
            )
        obabel = shutil.which("obabel")
        if not obabel:
            raise RuntimeError("Open Babel (`obabel`) is required to build ligand PDBQT from SMILES")
        cid = str(getattr(row, "compound_id"))
        smi = tmp_dir / f"{cid}.smi"
        pdbqt = tmp_dir / f"{cid}.pdbqt"
        smi.write_text(f"{smiles}\t{cid}\n", encoding="utf-8")
        cp = run_cmd(
            [obabel, "-ismi", str(smi), "-opdbqt", "-O", str(pdbqt), "--gen3d"],
            timeout_s=120.0,
        )
        if cp.returncode != 0 or not pdbqt.exists():
            raise RuntimeError(f"Failed ligand PDBQT conversion for {cid}: {cp.stderr}")
        return pdbqt

    def _dock_one(
        self,
        *,
        vina_bin: str,
        receptor_pdbqt: Path,
        poses_dir: Path,
        logs_dir: Path,
        tmp_dir: Path,
        row,
        center,
        size,
        vina_seed: int,
        exhaustiveness: int,
    ) -> dict[str, object]:
        cid = str(getattr(row, "compound_id"))
        ligand_pdbqt = self._prepare_ligand_pdbqt(row, tmp_dir)
        out_pose = poses_dir / f"{cid}.pdbqt"
        out_log = logs_dir / f"{cid}.vina.log"
        argv = [
            vina_bin,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--center_x",
            str(float(center[0])),
            "--center_y",
            str(float(center[1])),
            "--center_z",
            str(float(center[2])),
            "--size_x",
            str(float(size[0])),
            "--size_y",
            str(float(size[1])),
            "--size_z",
            str(float(size[2])),
            "--cpu",
            "1",
            "--seed",
            str(int(vina_seed)),
            "--exhaustiveness",
            str(int(exhaustiveness)),
            "--out",
            str(out_pose),
            "--log",
            str(out_log),
        ]
        cp = run_cmd(argv, timeout_s=600.0)
        if cp.returncode != 0:
            raise RuntimeError(
                f"vina failed for {cid} (exit {cp.returncode}). stderr: {cp.stderr[:3000]}"
            )
        log_text = out_log.read_text(encoding="utf-8", errors="replace") if out_log.exists() else (cp.stdout or "")
        score = parse_vina_score(log_text)
        return {
            "compound_id": cid,
            "score": score,
            "pose_path": str(out_pose),
            "backend": self.name,
        }
