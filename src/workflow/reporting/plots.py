"""Generate PNG summaries for a run directory (plan Human review plots)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def generate_all(run_dir: Path) -> None:
    run_dir = Path(run_dir)
    plot_dir = run_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    if (run_dir / "summary" / "summary.json").exists():
        plot_funnel(run_dir, plot_dir)
    hits = run_dir / "pubchem" / "hits.parquet"
    if hits.exists():
        plot_similarity_hist(hits, plot_dir)
    wt = run_dir / "poses" / "wt" / "docking_scores.parquet"
    if wt.exists():
        plot_dock_hist(wt, plot_dir)
    val = run_dir / "validation" / "ground_truth_report.json"
    if val.exists():
        plot_ground_truth(val, plot_dir)
    plt.close("all")


def plot_funnel(run_dir: Path, plot_dir: Path) -> None:
    s = json.loads((run_dir / "summary" / "summary.json").read_text(encoding="utf-8"))
    if s.get("dry_run"):
        return
    labels = ["PubChem hits", "After Tier A+B", "Dock pool", "Docked WT"]
    keys = ["n_pubchem_hits", "n_after_tier_ab", "n_dock_pool", "n_docked_wt"]
    values = [s.get(k, 0) for k in keys]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, values, color="steelblue")
    ax.set_ylabel("Count")
    ax.set_title("Pipeline funnel")
    plt.xticks(rotation=15, ha="right")
    fig.tight_layout()
    fig.savefig(plot_dir / "funnel.png", dpi=150)
    plt.close(fig)


def plot_similarity_hist(hits_path: Path, plot_dir: Path) -> None:
    df = pd.read_parquet(hits_path, columns=["tanimoto_similarity"])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["tanimoto_similarity"], bins=20, color="seagreen", edgecolor="white")
    ax.set_xlabel("Tanimoto (2D)")
    ax.set_ylabel("Count")
    ax.set_title("PubChem similarity")
    fig.tight_layout()
    fig.savefig(plot_dir / "similarity_hist.png", dpi=150)
    plt.close(fig)


def plot_dock_hist(scores_path: Path, plot_dir: Path) -> None:
    df = pd.read_parquet(scores_path, columns=["score"])
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df["score"], bins=20, color="slateblue", edgecolor="white")
    ax.set_xlabel("Score")
    ax.set_title("WT docking (primary score)")
    fig.tight_layout()
    fig.savefig(plot_dir / "dock_wt_hist.png", dpi=150)
    plt.close(fig)


def plot_ground_truth(report_path: Path, plot_dir: Path) -> None:
    rep = json.loads(report_path.read_text(encoding="utf-8"))
    df = pd.DataFrame(rep["cases"])
    vc = df["passed"].map({True: "pass", False: "fail"}).value_counts()
    vc = vc.reindex(["pass", "fail"]).fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(4, 4))
    vc.plot.bar(ax=ax, color=["green", "red"])
    ax.set_xticklabels(vc.index, rotation=0)
    ax.set_ylabel("Cases")
    ax.set_title("Ground truth")
    fig.tight_layout()
    fig.savefig(plot_dir / "ground_truth_pass_fail.png", dpi=150)
    plt.close(fig)
