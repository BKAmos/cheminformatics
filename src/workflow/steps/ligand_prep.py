"""Ligand prep for docking: write multi-SMILES file (SDF optional with RDKit)."""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from workflow.contracts import WorkflowConfig
from workflow.logging_utils import log_step
from workflow.obabel_resolve import obabel_argv0
from workflow.subprocess_runner import run_cmd


def run_ligand_prep(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    pool = pd.read_parquet(paths["filters"] / "dock_pool.parquet")
    lig_dir = paths["ligands"] / "prepared_sdf"
    lig_dir.mkdir(parents=True, exist_ok=True)
    out_txt = lig_dir / "ligands.smi"
    lines = [f"{row.hit_smiles}\t{row.compound_id}" for row in pool.itertuples()]
    out_txt.write_text("\n".join(lines), encoding="utf-8")
    index_rows: list[dict[str, str | None]] = [
        {"compound_id": str(r.compound_id), "ligand_pdbqt_path": None} for r in pool.itertuples()
    ]

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        w = Chem.SDWriter(str(lig_dir / "ligands.sdf"))
        for row in pool.itertuples():
            m = Chem.MolFromSmiles(row.hit_smiles)
            if m:
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE)
                w.write(m)
        w.close()
    except ImportError:
        pass

    # Optional fast path for Vina: precompute ligand PDBQT and persist paths.
    ob = obabel_argv0()
    if ob:
        pdbqt_dir = paths["ligands"] / "prepared_pdbqt"
        pdbqt_dir.mkdir(parents=True, exist_ok=True)
        smi_dir = paths["ligands"] / "prepared_smi"
        smi_dir.mkdir(parents=True, exist_ok=True)
        idx_map: dict[str, str] = {}
        for row in pool.itertuples():
            cid = str(row.compound_id)
            smi = smi_dir / f"{cid}.smi"
            pdbqt = pdbqt_dir / f"{cid}.pdbqt"
            smi.write_text(f"{row.hit_smiles}\t{cid}\n", encoding="utf-8")
            cp = run_cmd(
                [*ob, "-ismi", str(smi), "-opdbqt", "-O", str(pdbqt), "--gen3d"],
                timeout_s=120.0,
            )
            if cp.returncode == 0 and pdbqt.exists():
                idx_map[cid] = str(pdbqt)
        if pool.shape[0] and not idx_map:
            raise RuntimeError(
                "Open Babel failed to produce PDBQT for every ligand. "
                "Check `obabel`, SMILES validity, and run logs."
            )
        for rec in index_rows:
            cid = str(rec["compound_id"])
            rec["ligand_pdbqt_path"] = idx_map.get(cid)

    idx_df = pd.DataFrame(index_rows)
    idx_df.to_parquet(paths["ligands"] / "ligands_index.parquet", index=False)
    pool2 = pool.merge(idx_df, on="compound_id", how="left")
    pool2.to_parquet(paths["filters"] / "dock_pool.parquet", index=False)

    log_step(
        paths,
        "ligand_prep",
        time.perf_counter() - t0,
        input_count=len(pool),
        output_count=len(pool),
    )
    return out_txt
