from pathlib import Path

import pandas as pd

from workflow.config_load import load_workflow_config
from workflow.contracts import DockPoolCapConfig
from workflow.run_layout import ensure_run_dirs
from workflow.steps.dock_pool_cap import run_dock_pool_cap
from workflow.steps.tier_a import run_tier_a
from workflow.steps.tier_b import run_tier_b


def _minimal_hits(tmp: Path) -> None:
    pub = tmp / "pubchem"
    pub.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "compound_id": "c1",
                "query_smiles": "CCO",
                "cid": "1",
                "hit_smiles": "CCO",
                "tanimoto_similarity": 0.99,
                "lineage_query_id": "q1",
            },
            {
                "compound_id": "c2",
                "query_smiles": "CCO",
                "cid": "2",
                "hit_smiles": "CCCP",
                "tanimoto_similarity": 0.5,
                "lineage_query_id": "q1",
            },
        ]
    ).to_parquet(pub / "hits.parquet")


def test_top_n_cap(tmp_path: Path, repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(
        update={
            "dock_pool": DockPoolCapConfig(
                mode="top_n_similarity",
                top_n_by_2d_score=1,
                max_compounds_to_dock=None,
            )
        }
    )
    paths = ensure_run_dirs(tmp_path)
    _minimal_hits(tmp_path)
    run_tier_a(paths, cfg)
    run_tier_b(paths, cfg)
    run_dock_pool_cap(paths, cfg)
    pool = pd.read_parquet(paths["filters"] / "dock_pool.parquet")
    assert len(pool) == 1
    assert pool.iloc[0].tanimoto_similarity >= 0.9
