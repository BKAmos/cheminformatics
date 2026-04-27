"""GNINA backend — optional GPU docking; not implemented in scaffold."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


class GninaDockingBackend:
    name = "gnina"
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
        raise NotImplementedError(
            "Add gnina subprocess invocation and score parsing. "
            "Use MockDockingBackend until CUDA/GNINA are available."
        )
