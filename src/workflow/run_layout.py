"""Run directory layout (plan: Pipeline artifacts)."""

from __future__ import annotations

from pathlib import Path


def ensure_run_dirs(run_root: Path) -> dict[str, Path]:
    sub = {
        "root": run_root,
        "logs": run_root / "logs",
        "pubchem": run_root / "pubchem",
        "filters": run_root / "filters",
        "ligands": run_root / "ligands",
        "structure": run_root / "structure",
        "poses_wt": run_root / "poses" / "wt",
        "poses_mutant": run_root / "poses" / "mutant",
        "mutations": run_root / "mutations",
        "evolution": run_root / "evolution",
        "checkpoints": run_root / "checkpoints",
        "validation": run_root / "validation",
        "summary": run_root / "summary",
        "plots": run_root / "plots",
    }
    for p in sub.values():
        p.mkdir(parents=True, exist_ok=True)
    (sub["structure"] / "fpocket_raw").mkdir(parents=True, exist_ok=True)
    return sub


def manifest_path(run_root: Path) -> Path:
    return run_root / "run_manifest.json"
