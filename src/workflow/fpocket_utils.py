"""Parse fpocket output directories into pocket boxes (CPU, receptor-only)."""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from workflow.subprocess_runner import run_cmd


@dataclass
class ParsedPocket:
    pocket_id: int
    fpocket_score: float
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    atm_pdb: Path | None = None


_SCORE_LINE_RE = re.compile(r"^\s*Score\s*:\s*([+-]?\d+(?:\.\d+)?)", re.IGNORECASE)


def parse_pocket_info_score(info_text: str) -> float | None:
    best: float | None = None
    for line in info_text.splitlines():
        m = _SCORE_LINE_RE.match(line.strip())
        if m:
            val = float(m.group(1))
            best = val if best is None else max(best, val)
    return best


def _coords_from_pdb_atom_lines(pdb_path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        if len(line) < 54:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        coords.append((x, y, z))
    return coords


def bbox_from_coords(
    coords: list[tuple[float, float, float]],
    padding: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    if not coords:
        raise ValueError("no atomic coordinates in pocket PDB")
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    cx = (min_x + max_x) / 2.0
    cy = (min_y + max_y) / 2.0
    cz = (min_z + max_z) / 2.0
    pad = float(padding)
    sx = max(max_x - min_x + 2 * pad, 12.0)
    sy = max(max_y - min_y + 2 * pad, 12.0)
    sz = max(max_z - min_z + 2 * pad, 12.0)
    return (cx, cy, cz), (sx, sy, sz)


def discover_fpocket_output_dirs(fpocket_raw: Path, receptor_stem: str) -> list[Path]:
    """fpocket writes <stem>_out under cwd; also try nested matches."""
    candidates: list[Path] = []
    direct = fpocket_raw / f"{receptor_stem}_out"
    if direct.is_dir():
        candidates.append(direct)
    for p in sorted(fpocket_raw.glob("*_out")):
        if p.is_dir() and p not in candidates:
            candidates.append(p)
    return candidates


def parse_fpocket_output(
    out_root: Path,
    *,
    box_padding_angstrom: float,
) -> list[ParsedPocket]:
    pockets_dir = out_root / "pockets"
    if not pockets_dir.is_dir():
        return []

    parsed: list[ParsedPocket] = []
    for atm in sorted(pockets_dir.glob("pocket*_atm.pdb")):
        m = re.search(r"pocket(\d+)_atm\.pdb$", atm.name, re.IGNORECASE)
        if not m:
            continue
        pid = int(m.group(1))
        info = pockets_dir / f"pocket{pid}_info.txt"
        score = 0.0
        if info.is_file():
            sc = parse_pocket_info_score(info.read_text(encoding="utf-8", errors="replace"))
            if sc is not None:
                score = sc
        coords = _coords_from_pdb_atom_lines(atm)
        try:
            center, size = bbox_from_coords(coords, box_padding_angstrom)
        except ValueError:
            continue
        parsed.append(
            ParsedPocket(
                pocket_id=pid,
                fpocket_score=score,
                center=center,
                size=size,
                atm_pdb=atm,
            )
        )
    parsed.sort(key=lambda p: (-p.fpocket_score, p.pocket_id))
    return parsed


def run_fpocket_subprocess(
    *,
    receptor_pdb: Path,
    fpocket_raw: Path,
    executable: str | None,
    timeout_s: float = 600.0,
) -> None:
    fpocket_bin = executable or shutil.which("fpocket")
    if not fpocket_bin:
        raise RuntimeError(
            "fpocket binary not found. Install fpocket and ensure it is on PATH, "
            "or set fpocket.executable in config."
        )
    fpocket_raw.mkdir(parents=True, exist_ok=True)
    cp = run_cmd(
        [fpocket_bin, "-f", str(receptor_pdb.resolve())],
        cwd=fpocket_raw,
        timeout_s=timeout_s,
    )
    if cp.returncode != 0:
        raise RuntimeError(
            f"fpocket failed (exit {cp.returncode}). stderr: {(cp.stderr or '')[:4000]}"
        )


def pick_output_dir(fpocket_raw: Path, receptor_stem: str) -> Path:
    dirs = discover_fpocket_output_dirs(fpocket_raw, receptor_stem)
    if not dirs:
        raise RuntimeError(
            f"No fpocket output directory found under {fpocket_raw} "
            f"(expected something like {receptor_stem}_out/)."
        )
    return dirs[-1]
