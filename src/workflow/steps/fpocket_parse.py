"""Pocket specification: manual search box or fpocket-predicted cavities."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import pandas as pd

from workflow.contracts import PocketSpec, WorkflowConfig
from workflow.fpocket_utils import (
    parse_fpocket_output,
    pick_output_dir,
    run_fpocket_subprocess,
)
from workflow.logging_utils import log_step


def _write_pocket_spec(path: Path, spec: PocketSpec) -> None:
    path.write_text(spec.model_dump_json(indent=2), encoding="utf-8")


def _write_ranked_table(structure_dir: Path, rows: list[dict[str, object]]) -> None:
    df = pd.DataFrame(rows)
    df.to_parquet(structure_dir / "pockets_ranked.parquet", index=False)


def run_fpocket(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    fc = cfg.fpocket
    structure_dir = paths["structure"]
    raw_dir = structure_dir / "fpocket_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    receptor = paths["structure"] / "receptor_prepared.pdb"
    mode = fc.resolved_mode()

    if mode == "manual_box":
        spec = PocketSpec(
            center=fc.stub_center,
            size=fc.stub_size,
            source="manual_box",
            raw_path=None,
        )
        out = structure_dir / "pocket_spec.json"
        _write_pocket_spec(out, spec)
        _write_ranked_table(
            structure_dir,
            [
                {
                    "rank": 1,
                    "pocket_id": 1,
                    "fpocket_score": float("nan"),
                    "center_x": spec.center[0],
                    "center_y": spec.center[1],
                    "center_z": spec.center[2],
                    "size_x": spec.size[0],
                    "size_y": spec.size[1],
                    "size_z": spec.size[2],
                    "spec_file": "pocket_spec.json",
                }
            ],
        )
        (raw_dir / "manual_box.txt").write_text(
            "manual_box mode — no fpocket run.\n", encoding="utf-8"
        )
        log_step(paths, "fpocket_parse", time.perf_counter() - t0, output_count=1)
        return out

    # fpocket mode: run binary and parse pockets
    for child in raw_dir.iterdir():
        if child.is_dir():
            shutil.rmtree(child, ignore_errors=True)
        elif child.is_file():
            child.unlink(missing_ok=True)

    run_fpocket_subprocess(
        receptor_pdb=receptor,
        fpocket_raw=raw_dir,
        executable=fc.executable,
    )
    stem = receptor.stem
    out_root = pick_output_dir(raw_dir, stem)
    parsed = parse_fpocket_output(out_root, box_padding_angstrom=fc.box_padding_angstrom)
    if not parsed:
        raise RuntimeError(
            f"fpocket produced no parseable pockets under {out_root / 'pockets'}. "
            "Check receptor structure and fpocket logs."
        )

    spec_dir = structure_dir / "pocket_specs"
    spec_dir.mkdir(parents=True, exist_ok=True)
    ranked_rows: list[dict[str, object]] = []
    for rank, p in enumerate(parsed, start=1):
        fname = f"pocket_{p.pocket_id}.json"
        rel = f"pocket_specs/{fname}"
        spec = PocketSpec(
            center=p.center,
            size=p.size,
            source="fpocket",
            raw_path=p.atm_pdb,
        )
        _write_pocket_spec(spec_dir / fname, spec)
        ranked_rows.append(
            {
                "rank": rank,
                "pocket_id": p.pocket_id,
                "fpocket_score": p.fpocket_score,
                "center_x": p.center[0],
                "center_y": p.center[1],
                "center_z": p.center[2],
                "size_x": p.size[0],
                "size_y": p.size[1],
                "size_z": p.size[2],
                "spec_file": rel,
            }
        )

    primary = spec_dir / f"pocket_{parsed[0].pocket_id}.json"
    shutil.copy2(primary, structure_dir / "pocket_spec.json")
    _write_ranked_table(structure_dir, ranked_rows)

    log_step(paths, "fpocket_parse", time.perf_counter() - t0, output_count=len(ranked_rows))
    return structure_dir / "pocket_spec.json"
