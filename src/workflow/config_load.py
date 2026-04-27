"""Load workflow config from YAML + env."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from workflow.contracts import WorkflowConfig


def load_workflow_config(path: Path) -> WorkflowConfig:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    return _coerce_config(raw, base_dir=path.parent)


def _coerce_config(raw: dict[str, Any], base_dir: Path) -> WorkflowConfig:
    data = dict(raw)
    if "library_csv" in data:
        data["library_csv"] = (base_dir / data["library_csv"]).resolve()
    if "receptor_pdb" in data:
        data["receptor_pdb"] = (base_dir / data["receptor_pdb"]).resolve()
    if "fpocket" in data and isinstance(data.get("fpocket"), dict):
        sub = data["fpocket"]
        if "stub_center" in sub and isinstance(sub["stub_center"], list):
            sub["stub_center"] = tuple(sub["stub_center"])
        if "stub_size" in sub and isinstance(sub["stub_size"], list):
            sub["stub_size"] = tuple(sub["stub_size"])
    return WorkflowConfig.model_validate(data)
