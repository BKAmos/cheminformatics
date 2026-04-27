"""Docking backend protocol."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import pandas as pd


class DockingBackend(Protocol):
    name: str

    def dock_batch(
        self,
        *,
        receptor_pdb: Path,
        pocket_spec_path: Path,
        dock_pool: pd.DataFrame,
        poses_dir: Path,
        max_parallel: int,
        vina_seed: int = 42,
        vina_exhaustiveness: int = 8,
    ) -> pd.DataFrame:
        """Return DataFrame compound_id, score, pose_path, backend."""
        ...
