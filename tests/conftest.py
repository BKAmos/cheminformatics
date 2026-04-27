"""Pytest configuration."""

import re
from pathlib import Path

import pandas as pd
import pytest

pytest_plugins = ["pytest_httpx"]

_PUBCHEM_BODY = {
    "InformationList": {
        "Information": [
            {"CID": 702, "Similarity": 0.95, "CanonicalSMILES": "CCO"},
        ]
    }
}

_PUBCHEM_RE = re.compile(r"^https://pubchem\.ncbi\.nlm\.nih\.gov/.+")


@pytest.fixture
def mock_pubchem(httpx_mock) -> object:
    httpx_mock.add_response(
        url=_PUBCHEM_RE,
        json=_PUBCHEM_BODY,
        is_reusable=True,
    )
    return httpx_mock


@pytest.fixture(scope="session", autouse=True)
def _ensure_ground_truth_parquet_fixtures() -> None:
    d = Path(__file__).resolve().parent / "fixtures" / "ground_truth"
    d.mkdir(parents=True, exist_ok=True)
    sample = d / "ranking_sample.parquet"
    if not sample.is_file():
        pd.DataFrame({"composite_rank_score": [10.0, 7.5, 3.0]}).to_parquet(sample, index=False)
    pockets = d / "pockets_ranked_sample.parquet"
    if not pockets.is_file():
        pd.DataFrame(
            {
                "rank": [1, 2],
                "pocket_id": [1, 2],
                "fpocket_score": [12.0, 9.0],
                "center_x": [0.0, 1.0],
                "center_y": [0.0, 0.0],
                "center_z": [0.0, 0.0],
                "size_x": [10.0, 10.0],
                "size_y": [10.0, 10.0],
                "size_z": [10.0, 10.0],
                "spec_file": ["pocket_spec.json", "pocket_specs/pocket_2.json"],
            }
        ).to_parquet(pockets, index=False)


@pytest.fixture
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]
