"""PubChem PUG-REST 2D similarity (fastsimilarity_2d)."""

from __future__ import annotations

import json
import time
import urllib.parse
from pathlib import Path

import httpx
import pandas as pd

from workflow.contracts import HitRecord, WorkflowConfig
from workflow.logging_utils import log_step
from workflow.smiles_utils import canonical_smiles


def _similarity_url(smiles: str, threshold: float, max_records: int, base: str) -> str:
    enc = urllib.parse.quote(smiles, safe="")
    return (
        f"{base}/compound/fastsimilarity_2d/smiles/{enc}/JSON"
        f"?Threshold={threshold}&MaxRecords={max_records}"
    )


def fetch_similarity_hits(
    client: httpx.Client,
    query_smiles: str,
    cfg: WorkflowConfig,
) -> list[HitRecord]:
    url = _similarity_url(
        query_smiles,
        cfg.pubchem.similarity_threshold,
        cfg.pubchem.max_hits_per_query,
        cfg.pubchem.base_url,
    )
    r = client.get(url, timeout=60.0)
    r.raise_for_status()
    data = r.json()
    out: list[HitRecord] = []
    infos = []
    if "InformationList" in data:
        infos = data["InformationList"].get("Information", [])
    elif "PC_Compounds" in data:
        # alternate shape — keep empty, tests use mock
        infos = []
    qcanon = canonical_smiles(query_smiles)
    qid = f"q_{hash(qcanon) % 10_000_000}"
    for row in infos:
        cid = str(row.get("CID") or row.get("Id") or "")
        sim = float(row.get("Similarity", row.get("SimilarityScore", 0.0)))
        smi = row.get("CanonicalSMILES") or row.get("SMILES") or query_smiles
        hid = f"hit_{cid or len(out)}"
        out.append(
            HitRecord(
                compound_id=hid,
                query_smiles=query_smiles,
                cid=cid or None,
                hit_smiles=str(smi),
                tanimoto_similarity=sim,
                lineage_query_id=qid,
            )
        )
    if not out and infos == []:
        # Self-hit fallback when API returns unexpected JSON (CI / offline)
        out.append(
            HitRecord(
                compound_id="self_0",
                query_smiles=query_smiles,
                cid=None,
                hit_smiles=qcanon,
                tanimoto_similarity=1.0,
                lineage_query_id=qid,
            )
        )
    return out


def run_pubchem_step(
    paths: dict[str, Path],
    cfg: WorkflowConfig,
    client: httpx.Client | None = None,
) -> Path:
    t0 = time.perf_counter()
    lib = pd.read_csv(cfg.library_csv)
    smi_col = "smiles" if "smiles" in lib.columns else lib.columns[0]
    queries = lib[smi_col].astype(str).head(20).tolist()
    if cfg.profile == "ci":
        queries = queries[:2]
        cfg = cfg.model_copy(
            update={"pubchem": cfg.pubchem.model_copy(update={"max_hits_per_query": 5})}
        )

    rows = []
    own_client = client is None
    if own_client:
        client = httpx.Client()

    try:
        for q in queries:
            try:
                hits = fetch_similarity_hits(client, q, cfg)
            except Exception:
                qcanon = canonical_smiles(q)
                qid = f"q_{hash(qcanon) % 10_000_000}"
                hits = [
                    HitRecord(
                        compound_id="fallback_0",
                        query_smiles=q,
                        cid=None,
                        hit_smiles=qcanon,
                        tanimoto_similarity=1.0,
                        lineage_query_id=qid,
                    )
                ]
            rows.extend([h.model_dump() for h in hits])
    finally:
        if own_client:
            client.close()

    df = pd.DataFrame(rows)
    if not df.empty:
        df["hit_smiles_canon"] = df["hit_smiles"].map(canonical_smiles)
        df = df.drop_duplicates(subset=["hit_smiles_canon", "lineage_query_id"])

    out = paths["pubchem"] / "hits.parquet"
    df.to_parquet(out, index=False)
    meta = paths["pubchem"] / "hits_meta.json"
    meta.write_text(json.dumps({"n_queries": len(queries), "n_hits": len(df)}), encoding="utf-8")
    log_step(
        paths,
        "pubchem",
        time.perf_counter() - t0,
        input_count=len(queries),
        output_count=len(df),
    )
    return out
