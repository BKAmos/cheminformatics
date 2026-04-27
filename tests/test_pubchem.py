from pathlib import Path

import httpx
import pytest

from workflow.config_load import load_workflow_config
from workflow.run_layout import ensure_run_dirs
from workflow.steps.pubchem import fetch_similarity_hits, run_pubchem_step


def test_fetch_similarity_parses_information_list(repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    body = {
        "InformationList": {
            "Information": [
                {"CID": 702, "Similarity": 1.0, "CanonicalSMILES": "CCO"},
            ]
        }
    }
    transport = httpx.MockTransport(lambda r: httpx.Response(200, json=body))
    with httpx.Client(transport=transport) as client:
        hits = fetch_similarity_hits(client, "CCO", cfg)
    assert len(hits) >= 1
    assert hits[0].hit_smiles == "CCO"


def test_run_pubchem_step_mocked(tmp_path: Path, repo_root: Path) -> None:
    cfg = load_workflow_config(repo_root / "configs" / "ci.yaml")
    cfg = cfg.model_copy(update={"library_csv": repo_root / "data" / "library.csv"})
    paths = ensure_run_dirs(tmp_path)
    body = {
        "InformationList": {
            "Information": [
                {"CID": 702, "Similarity": 0.99, "CanonicalSMILES": "CCO"},
            ]
        }
    }
    transport = httpx.MockTransport(lambda r: httpx.Response(200, json=body))
    with httpx.Client(transport=transport) as client:
        out = run_pubchem_step(paths, cfg, client=client)
    assert out.exists()
    import pandas as pd

    df = pd.read_parquet(out)
    assert len(df) >= 1
