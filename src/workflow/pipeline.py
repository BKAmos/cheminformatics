"""LangGraph orchestration + sequential fallback."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, TypedDict

import pandas as pd
from langgraph.graph import END, StateGraph

from workflow.contracts import WorkflowConfig
from workflow.run_layout import ensure_run_dirs, manifest_path
from workflow.steps.candidate_ranking import run_candidate_ranking
from workflow.steps.dock_pool_cap import run_dock_pool_cap
from workflow.steps.docking_mutant import run_dock_mutant
from workflow.steps.docking_wt import run_dock_wt
from workflow.steps.evolution import run_evolution
from workflow.steps.evolved_refilter import run_evolved_refilter
from workflow.steps.fpocket_parse import run_fpocket
from workflow.steps.ligand_prep import run_ligand_prep
from workflow.steps.openmm_mutate import run_openmm_mutate
from workflow.steps.pubchem import run_pubchem_step
from workflow.steps.receptor_prep import run_receptor_prep
from workflow.steps.tier_a import run_tier_a
from workflow.steps.tier_b import run_tier_b


class GraphState(TypedDict, total=False):
    ok: bool


def write_manifest(run_dir: Path, cfg: WorkflowConfig, extra: dict[str, Any]) -> None:
    payload = {
        "workflow_version": "0.1.0",
        "config_snapshot": json.loads(cfg.model_dump_json()),
        **extra,
    }
    manifest_path(run_dir).write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")


def write_summary(paths: dict[str, Path], cfg: WorkflowConfig) -> None:
    hits = paths["pubchem"] / "hits.parquet"
    tierb = paths["filters"] / "tier_b_rationale.parquet"
    pool = paths["filters"] / "dock_pool.parquet"
    dock = paths["poses_wt"] / "docking_scores.parquet"
    cand = paths["summary"] / "candidates_ranked.parquet"
    n_pub = len(pd.read_parquet(hits)) if hits.exists() else 0
    n_tier = len(pd.read_parquet(tierb)) if tierb.exists() else 0
    n_pool = len(pd.read_parquet(pool)) if pool.exists() else 0
    n_dock = len(pd.read_parquet(dock)) if dock.exists() else 0
    n_candidates = 0
    if cand.exists():
        cdf = pd.read_parquet(cand)
        n_candidates = len(cdf) if not cdf.empty else 0
    passed = 0
    if tierb.exists() and n_tier:
        passed = int(pd.read_parquet(tierb)["passed_tier_ab"].sum())
    summary = {
        "n_pubchem_hits": n_pub,
        "n_after_tier_ab": passed,
        "n_dock_pool": n_pool,
        "n_docked_wt": n_dock,
        "n_ranked_candidates": n_candidates,
    }
    paths["summary"].mkdir(parents=True, exist_ok=True)
    (paths["summary"] / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


def run_dry_run(cfg: WorkflowConfig, run_dir: Path) -> None:
    paths = ensure_run_dirs(run_dir)
    write_manifest(run_dir, cfg, {"run_dir": str(run_dir.resolve()), "dry_run": True})
    (paths["summary"] / "summary.json").write_text('{"dry_run": true}', encoding="utf-8")


def compile_workflow_graph(cfg: WorkflowConfig, run_dir: Path):
    paths = ensure_run_dirs(run_dir)

    def wrap(step_fn):
        def _inner(_: GraphState) -> GraphState:
            step_fn(paths, cfg)
            return {"ok": True}

        return _inner

    def node_manifest(_: GraphState) -> GraphState:
        write_manifest(run_dir, cfg, {"run_dir": str(run_dir.resolve())})
        return {"ok": True}

    g = StateGraph(GraphState)
    g.add_node("manifest", node_manifest)
    g.add_node("pubchem", wrap(run_pubchem_step))
    g.add_node("tier_a", wrap(run_tier_a))
    g.add_node("tier_b", wrap(run_tier_b))
    g.add_node("dock_pool", wrap(run_dock_pool_cap))
    g.add_node("receptor", wrap(run_receptor_prep))
    g.add_node("fpocket", wrap(run_fpocket))
    g.add_node("ligprep", wrap(run_ligand_prep))
    g.add_node("dock_wt", wrap(run_dock_wt))
    g.add_node("mutate", wrap(run_openmm_mutate))
    g.add_node("dock_mut", wrap(run_dock_mutant))
    g.add_node("evolve", wrap(run_evolution))
    g.add_node("evorefilt", wrap(run_evolved_refilter))
    g.add_node("candidates", wrap(run_candidate_ranking))
    g.add_node("summary", wrap(write_summary))

    g.set_entry_point("manifest")
    for a, b in [
        ("manifest", "pubchem"),
        ("pubchem", "tier_a"),
        ("tier_a", "tier_b"),
        ("tier_b", "dock_pool"),
        ("dock_pool", "receptor"),
        ("receptor", "fpocket"),
        ("fpocket", "ligprep"),
        ("ligprep", "dock_wt"),
        ("dock_wt", "mutate"),
        ("mutate", "dock_mut"),
        ("dock_mut", "evolve"),
        ("evolve", "evorefilt"),
        ("evorefilt", "candidates"),
        ("candidates", "summary"),
    ]:
        g.add_edge(a, b)
    g.add_edge("summary", END)
    return g.compile()


def run_pipeline(cfg: WorkflowConfig, run_dir: Path, *, use_langgraph: bool = True) -> None:
    if cfg.dry_run:
        run_dry_run(cfg, run_dir)
        return
    if use_langgraph:
        app = compile_workflow_graph(cfg, run_dir)
        app.invoke({})
    else:
        paths = ensure_run_dirs(run_dir)
        write_manifest(run_dir, cfg, {"run_dir": str(run_dir.resolve())})
        run_pubchem_step(paths, cfg)
        run_tier_a(paths, cfg)
        run_tier_b(paths, cfg)
        run_dock_pool_cap(paths, cfg)
        run_receptor_prep(paths, cfg)
        run_fpocket(paths, cfg)
        run_ligand_prep(paths, cfg)
        run_dock_wt(paths, cfg)
        run_openmm_mutate(paths, cfg)
        run_dock_mutant(paths, cfg)
        run_evolution(paths, cfg)
        run_evolved_refilter(paths, cfg)
        run_candidate_ranking(paths, cfg)
        write_summary(paths, cfg)
