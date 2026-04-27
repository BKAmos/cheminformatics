"""Receptor mutation: copy (skip) or apply point mutations via PDBFixer/OpenMM."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

from workflow.contracts import MutationResidue, WorkflowConfig
from workflow.logging_utils import log_step


def _resname_at(pdb_path: Path, residue: MutationResidue) -> str:
    want_chain = residue.chain.strip()
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 54:
            continue
        chain = line[21].strip()
        if chain != want_chain:
            continue
        try:
            seq = int(line[22:26])
        except ValueError:
            continue
        if seq != residue.resseq:
            continue
        return line[17:20].strip()
    raise ValueError(
        f"No ATOM/HETATM residue found for chain={want_chain!r} "
        f"resseq={residue.resseq} in {pdb_path}"
    )


def run_openmm_mutate(paths: dict[str, Path], cfg: WorkflowConfig) -> Path:
    t0 = time.perf_counter()
    wt = paths["structure"] / "receptor_prepared.pdb"
    mut = paths["structure"] / "mutant_receptor.pdb"
    meta = paths["mutations"] / "residue_map.json"

    if cfg.skip_mutation or not cfg.mutation.residues:
        shutil.copy2(wt, mut)
        meta.write_text(
            json.dumps({"note": "no_mutation", "skipped": bool(cfg.skip_mutation)}),
            encoding="utf-8",
        )
        log_step(paths, "openmm_mutate", time.perf_counter() - t0, output_count=1)
        return mut

    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as e:
        raise RuntimeError(
            "Real mutations require optional dependencies. Install with: "
            'pip install -e ".[mutate]"'
        ) from e

    mut_strings: list[str] = []
    provenance: list[dict[str, str | int]] = []
    for r in cfg.mutation.residues:
        old = _resname_at(wt, r)
        ch = r.chain.strip()
        mut_strings.append(f"{old}-{r.resseq}-{ch}-{r.to_aa}")
        provenance.append(
            {
                "chain": ch,
                "resseq": r.resseq,
                "from_aa": old,
                "to_aa": r.to_aa,
            }
        )

    fixer = PDBFixer(filename=str(wt))
    fixer.applyMutations(mut_strings)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with mut.open("w", encoding="utf-8") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    meta.write_text(
        json.dumps(
            {
                "note": "pdbfixer_applyMutations",
                "mutations": provenance,
                "skipped": False,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    log_step(paths, "openmm_mutate", time.perf_counter() - t0, output_count=1)
    return mut
