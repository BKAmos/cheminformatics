"""Receptor prep: copy or placeholder (PDBFixer/OpenMM optional on full install)."""

from __future__ import annotations

import shutil
import time
from pathlib import Path

from workflow.artifacts import validate_receptor_pdb
from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step


def run_receptor_prep(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    src = cfg.receptor_pdb
    validate_receptor_pdb(src)
    dst = paths["structure"] / "receptor_prepared.pdb"
    shutil.copy2(src, dst)
    log_step(paths, "receptor_prep", time.perf_counter() - t0, output_count=1)
    return dst
