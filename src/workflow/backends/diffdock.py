"""DiffDock (or compatible) ML docking — optional; subprocess/API not wired in scaffold."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


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
        raise NotImplementedError(
            "Wire DiffDock (or your ML pose engine): receptor, ligands, pocket box, output poses. "
            "Use ml_backend: mock until the CLI/container is available."
        )
