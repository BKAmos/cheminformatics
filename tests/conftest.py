"""Pytest configuration."""

from pathlib import Path

import httpx
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


@pytest.fixture
def mock_pubchem(httpx_mock):
    def _handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=_PUBCHEM_BODY)

    httpx_mock.add_callback(_handler)
    return httpx_mock


@pytest.fixture(scope="session", autouse=True)
def _ensure_ground_truth_parquet_fixtures() -> None:
    d = Path(__file__).resolve().parent / "fixtures" / "ground_truth"
    d.mkdir(parents=True, exist_ok=True)
    sample = d / "ranking_sample.parquet"
    if not sample.is_file():
        pd.DataFrame({"composite_rank_score": [10.0, 7.5, 3.0]}).to_parquet(sample, index=False)


@pytest.fixture
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]
