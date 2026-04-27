"""DiffDock (or compatible) ML docking.

This backend keeps the pipeline end-to-end testable: by default it writes deterministic
``diffdock``-tagged scores and placeholder poses (same contract as
:class:`workflow.backends.mock.MockMLDockingBackend`). For a real DiffDock or compatible CLI,
set ``CHEM_WORKFLOW_DIFFDOCK_CMD`` to a program that is invoked with environment variables:
``RECEPTOR_PDB``, ``POCKET_SPEC_JSON``, ``LIGAND_SMILES``, ``OUT_POSE_PDB`` (one compound per
process). A zero exit code and an existing output pose is required; stderr is ignored.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from workflow.subprocess_runner import run_cmd


def _score_one_diffdock(
    *,
    row: object,
    receptor_pdb: Path,
    pocket_spec_path: Path,
    poses_dir: Path,
    ext_cmd: str | None,
) -> dict[str, object]:
    spec = json.loads(pocket_spec_path.read_text(encoding="utf-8"))
    cid = str(getattr(row, "compound_id"))
    pose = poses_dir / f"{cid}_diffdock.pdb"
    if ext_cmd:
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td) / f"{cid}_out.pdb"
            env = os.environ.copy()
            env["RECEPTOR_PDB"] = str(receptor_pdb.resolve())
            env["POCKET_SPEC_JSON"] = str(pocket_spec_path.resolve())
            lig = getattr(row, "ligand_smiles", None) or getattr(row, "smiles", None) or ""
            env["LIGAND_SMILES"] = str(lig)
            env["OUT_POSE_PDB"] = str(tdp)
            p = run_cmd([ext_cmd], timeout_s=1200.0, env=env)
            if p.returncode != 0 or not tdp.is_file():
                err_tail = (p.stderr or "")[:500]
                raise RuntimeError(
                    f"DiffDock command failed for {cid}: code={p.returncode} stderr={err_tail!r}"
                )
            shutil.copy2(tdp, pose)
    else:
        h0 = int(hashlib.sha256(f"diffdock:{cid}".encode()).hexdigest()[:8], 16)
        center0 = float(spec.get("center", [0.0, 0.0, 0.0])[0])
        sc = -6.0 - (h0 % 90) / 100.0 + center0 * 0.0003
        pose.write_text(
            f"DIFFDOCK_PLACEHOLDER {cid} score={sc:.4f}\n",
            encoding="utf-8",
        )

    h = int(hashlib.sha256(f"diffdock:{cid}".encode()).hexdigest()[:8], 16)
    center0 = float(spec.get("center", [0.0, 0.0, 0.0])[0])
    score = -6.0 - (h % 90) / 100.0 + center0 * 0.0003
    return {
        "compound_id": cid,
        "score": float(score),
        "pose_path": str(pose),
        "backend": "diffdock",
    }


class DiffDockBackend:
    name = "diffdock"
    family = "ml"

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
        del vina_seed, vina_exhaustiveness  # ML backend; kept for API parity
        poses_dir.mkdir(parents=True, exist_ok=True)
        ext = os.environ.get("CHEM_WORKFLOW_DIFFDOCK_CMD")
        if ext and not shutil.which(ext) and not Path(ext).is_file():
            raise RuntimeError(
                f"CHEM_WORKFLOW_DIFFDOCK_CMD={ext!r} is not an executable on PATH. "
                "Unset it to use the built-in placeholder DiffDock scores."
            )
        # Allow full path in env
        if ext and Path(ext).is_file():
            ext_path: str | None = str(Path(ext).resolve())
        else:
            ext_path = ext

        pool = dock_pool.sort_values("compound_id", kind="mergesort")
        workers = max(1, int(max_parallel))
        rows: list[dict[str, object]] = []

        def one(row: object) -> dict[str, object]:
            return _score_one_diffdock(
                row=row,
                receptor_pdb=receptor_pdb,
                pocket_spec_path=pocket_spec_path,
                poses_dir=poses_dir,
                ext_cmd=ext_path,
            )

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(one, r) for r in pool.itertuples()]
            for f in as_completed(futs):
                rows.append(f.result())

        out = pd.DataFrame(rows)
        return out.sort_values("compound_id", kind="mergesort").reset_index(drop=True)
