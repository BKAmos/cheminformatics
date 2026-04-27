"""Validate artifacts between pipeline steps (schemas and required files)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def validate_receptor_pdb(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Receptor PDB not found: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()[:8000]
    if not any(
        ln.startswith("ATOM") or ln.startswith("HETATM") or ln.startswith("MODEL")
        for ln in lines
    ):
        raise ValueError(f"Receptor file has no ATOM/HETATM records: {path}")


def validate_dock_pool_columns(df: pd.DataFrame) -> None:
    required = {"compound_id", "hit_smiles"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"dock_pool missing required columns {sorted(missing)}")


def validate_pocket_spec_dict(data: dict[str, Any]) -> None:
    if "center" not in data or "size" not in data:
        raise ValueError("pocket spec must include 'center' and 'size'")
    c, s = data["center"], data["size"]
    if not (isinstance(c, (list, tuple)) and len(c) == 3):
        raise ValueError("pocket spec 'center' must be length-3")
    if not (isinstance(s, (list, tuple)) and len(s) == 3):
        raise ValueError("pocket spec 'size' must be length-3")


def validate_pocket_spec_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Pocket spec not found: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    validate_pocket_spec_dict(data)


def read_json_if_exists(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))
