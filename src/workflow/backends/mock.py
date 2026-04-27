"""CPU-only mock docking for CI and laptops without Vina/GNINA."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


class MockDockingBackend:
    name = "mock"
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
        poses_dir.mkdir(parents=True, exist_ok=True)
        spec = json.loads(pocket_spec_path.read_text(encoding="utf-8"))
        rows = []
        pool = dock_pool.sort_values("compound_id", kind="mergesort")
        for row in pool.itertuples():
            h = int(hashlib.sha256(str(row.compound_id).encode()).hexdigest()[:8], 16)
            score = -8.0 - (h % 100) / 100.0 + float(spec.get("center", [0])[0]) * 0.001
            pose = poses_dir / f"{row.compound_id}.pdb"
            pose.write_text(f"MOCK_POSE {row.compound_id}\n", encoding="utf-8")
            rows.append(
                {
                    "compound_id": row.compound_id,
                    "score": score,
                    "pose_path": str(pose),
                    "backend": self.name,
                }
            )
        return pd.DataFrame(rows)


class MockMLDockingBackend:
    """CPU-only stand-in for ML docking (distinct scores from physics mock for tandem tests)."""

    name = "mock_ml"
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
        poses_dir.mkdir(parents=True, exist_ok=True)
        spec = json.loads(pocket_spec_path.read_text(encoding="utf-8"))
        rows = []
        pool = dock_pool.sort_values("compound_id", kind="mergesort")
        for row in pool.itertuples():
            h = int(hashlib.sha256(f"ml:{row.compound_id}".encode()).hexdigest()[:8], 16)
            # Different scale than physics mock so tandem merge is visibly dual-sourced.
            score = -6.5 - (h % 80) / 100.0 + float(spec.get("center", [0])[0]) * 0.0005
            pose = poses_dir / f"{row.compound_id}_ml.pdb"
            pose.write_text(f"MOCK_ML_POSE {row.compound_id}\n", encoding="utf-8")
            rows.append(
                {
                    "compound_id": row.compound_id,
                    "score": score,
                    "pose_path": str(pose),
                    "backend": self.name,
                }
            )
        return pd.DataFrame(rows)
