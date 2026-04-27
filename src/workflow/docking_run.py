"""Run physics and/or ML docking from config and merge results for downstream steps."""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd

from workflow.artifacts import validate_pocket_spec_file
from workflow.backends.factory import get_ml_backend, get_physics_backend
from workflow.contracts import DockingConfig, WorkflowConfig


def merge_docking_results(
    physics: pd.DataFrame | None,
    ml: pd.DataFrame | None,
    mode: str,
    ranking_when_both: str,
) -> pd.DataFrame:
    if mode == "physics_only":
        if physics is None:
            raise ValueError("physics_only mode requires a physics docking result")
        out = physics.copy()
        return out

    if mode == "ml_only":
        if ml is None:
            raise ValueError("ml_only mode requires an ML docking result")
        return ml.copy()

    if mode == "both":
        if physics is None or ml is None:
            raise ValueError("both mode requires physics and ML docking results")
        p = physics.rename(
            columns={
                "score": "score_physics",
                "pose_path": "pose_path_physics",
                "backend": "backend_physics",
            }
        )
        m = ml.rename(
            columns={
                "score": "score_ml",
                "pose_path": "pose_path_ml",
                "backend": "backend_ml",
            }
        )
        merged = p.merge(m, on="compound_id", how="inner")
        if ranking_when_both == "physics":
            merged["score"] = merged["score_physics"]
            merged["pose_path"] = merged["pose_path_physics"]
        else:
            merged["score"] = merged["score_ml"]
            merged["pose_path"] = merged["pose_path_ml"]
        merged["backend"] = merged["backend_physics"].astype(str) + "|" + merged["backend_ml"].astype(str)
        return merged

    raise ValueError(f"Unknown docking mode: {mode!r}")


def effective_top_k_pockets(cfg: WorkflowConfig) -> int:
    k = cfg.fpocket.top_k_pockets
    if cfg.docking.mode == "both" and k > 1:
        warnings.warn(
            "docking.mode=both with fpocket.top_k_pockets>1 is not supported; using top_k=1.",
            stacklevel=2,
        )
        return 1
    return k


def merge_multipocket_merged_frames(frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not frames:
        raise ValueError("no docking frames to merge")
    if len(frames) == 1:
        return frames[0].copy()
    all_df = pd.concat(frames, ignore_index=True)
    idx = all_df.groupby("compound_id")["score"].idxmin()
    return all_df.loc[idx].reset_index(drop=True)


def _poses_dir_for_family(poses_root: Path, docking: DockingConfig, family: str) -> Path:
    if docking.mode == "both":
        d = poses_root / family
        d.mkdir(parents=True, exist_ok=True)
        return d
    poses_root.mkdir(parents=True, exist_ok=True)
    return poses_root


def run_configured_docking(
    cfg: WorkflowConfig,
    *,
    dock_pool: pd.DataFrame,
    receptor_pdb: Path,
    pocket_spec_path: Path,
    poses_root: Path,
    score_offset: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame | None]:
    """Return merged scores table and optional raw physics / ML frames (for sidecar Parquet)."""
    dc = cfg.docking
    max_p = cfg.resources.max_parallel_docks
    physics_df: pd.DataFrame | None = None
    ml_df: pd.DataFrame | None = None
    pool = dock_pool.sort_values("compound_id", kind="mergesort").reset_index(drop=True)
    dock_kw = {
        "vina_seed": dc.vina_seed,
        "vina_exhaustiveness": dc.vina_exhaustiveness,
    }

    if dc.mode in ("physics_only", "both"):
        b = get_physics_backend(dc.physics_backend)
        physics_df = b.dock_batch(
            receptor_pdb=receptor_pdb,
            pocket_spec_path=pocket_spec_path,
            dock_pool=pool,
            poses_dir=_poses_dir_for_family(poses_root, dc, "physics"),
            max_parallel=max_p,
            **dock_kw,
        )
    if dc.mode in ("ml_only", "both"):
        b = get_ml_backend(dc.ml_backend)
        ml_df = b.dock_batch(
            receptor_pdb=receptor_pdb,
            pocket_spec_path=pocket_spec_path,
            dock_pool=pool,
            poses_dir=_poses_dir_for_family(poses_root, dc, "ml"),
            max_parallel=max_p,
            **dock_kw,
        )

    merged = merge_docking_results(physics_df, ml_df, dc.mode, dc.ranking_when_both)
    if score_offset:
        merged["score"] = merged["score"] + score_offset
        if "score_physics" in merged.columns:
            merged["score_physics"] = merged["score_physics"] + score_offset
        if "score_ml" in merged.columns:
            merged["score_ml"] = merged["score_ml"] + score_offset

    return merged, physics_df, ml_df


def run_configured_docking_with_pockets(
    cfg: WorkflowConfig,
    *,
    dock_pool: pd.DataFrame,
    receptor_pdb: Path,
    structure_dir: Path,
    poses_root: Path,
    score_offset: float = 0.0,
) -> tuple[pd.DataFrame, pd.DataFrame | None, pd.DataFrame | None, bool]:
    """Dock using ``pockets_ranked.parquet`` when top-K>1, else single ``pocket_spec.json``.

    Returns merged, physics_df, ml_df, multipocket_ran.
    """
    ranked_path = structure_dir / "pockets_ranked.parquet"
    if not ranked_path.is_file():
        pocket = structure_dir / "pocket_spec.json"
        validate_pocket_spec_file(pocket)
        m, p, l = run_configured_docking(
            cfg,
            dock_pool=dock_pool,
            receptor_pdb=receptor_pdb,
            pocket_spec_path=pocket,
            poses_root=poses_root,
            score_offset=score_offset,
        )
        return m, p, l, False

    ranked = pd.read_parquet(ranked_path).sort_values("rank", kind="mergesort").reset_index(drop=True)
    k = effective_top_k_pockets(cfg)
    if len(ranked) <= 1 or k <= 1:
        row0 = ranked.iloc[0]
        pocket = structure_dir / str(row0["spec_file"])
        validate_pocket_spec_file(pocket)
        m, p, l = run_configured_docking(
            cfg,
            dock_pool=dock_pool,
            receptor_pdb=receptor_pdb,
            pocket_spec_path=pocket,
            poses_root=poses_root,
            score_offset=score_offset,
        )
        return m, p, l, False

    merged_frames: list[pd.DataFrame] = []
    phys_frames: list[pd.DataFrame] = []
    ml_frames: list[pd.DataFrame] = []
    for _, row in ranked.head(k).iterrows():
        pocket = structure_dir / str(row["spec_file"])
        validate_pocket_spec_file(pocket)
        pid = int(row["pocket_id"])
        pr = poses_root / f"_pocket_run_{pid}"
        m, p, l = run_configured_docking(
            cfg,
            dock_pool=dock_pool,
            receptor_pdb=receptor_pdb,
            pocket_spec_path=pocket,
            poses_root=pr,
            score_offset=score_offset,
        )
        m2 = m.copy()
        m2["fpocket_pocket_id"] = pid
        m2["fpocket_pocket_rank"] = int(row["rank"])
        merged_frames.append(m2)
        if p is not None:
            phys_frames.append(p)
        if l is not None:
            ml_frames.append(l)

    merged = merge_multipocket_merged_frames(merged_frames)
    phys_out = pd.concat(phys_frames, ignore_index=True) if phys_frames else None
    ml_out = pd.concat(ml_frames, ignore_index=True) if ml_frames else None
    return merged, phys_out, ml_out, True


def write_docking_artifacts(
    poses_root: Path,
    cfg: WorkflowConfig,
    merged: pd.DataFrame,
    physics_df: pd.DataFrame | None,
    ml_df: pd.DataFrame | None,
) -> Path:
    """Write merged docking_scores plus optional per-family tables."""
    poses_root.mkdir(parents=True, exist_ok=True)
    main = poses_root / "docking_scores.parquet"
    merged.to_parquet(main, index=False)
    dc = cfg.docking
    if physics_df is not None and dc.mode in ("physics_only", "both"):
        physics_df.to_parquet(poses_root / "docking_scores_physics.parquet", index=False)
    if ml_df is not None and dc.mode in ("ml_only", "both"):
        ml_df.to_parquet(poses_root / "docking_scores_ml.parquet", index=False)
    return main
