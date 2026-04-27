"""Wild-type receptor docking pass."""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from workflow.artifacts import validate_dock_pool_columns
from workflow.contracts import WorkflowConfig
from workflow.docking_run import run_configured_docking_with_pockets, write_docking_artifacts
from workflow.logging_utils import log_step


def run_dock_wt(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    pool = pd.read_parquet(paths["filters"] / "dock_pool.parquet")
    validate_dock_pool_columns(pool)
    rec = paths["structure"] / "receptor_prepared.pdb"
    poses = paths["poses_wt"]
    merged, phys, ml, multi = run_configured_docking_with_pockets(
        cfg,
        dock_pool=pool,
        receptor_pdb=rec,
        structure_dir=paths["structure"],
        poses_root=poses,
        score_offset=0.0,
    )
    out = write_docking_artifacts(
        poses,
        cfg,
        merged,
        None if multi else phys,
        None if multi else ml,
    )
    log_step(paths, "docking_wt", time.perf_counter() - t0, input_count=len(pool), output_count=len(merged))
    return out
