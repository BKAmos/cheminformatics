"""Minimal `obabel` stand-in for environments where `obabel.exe` is missing (e.g. conda win-64).

Converts PDB/SMILES to PDBQT via RDKit (3D coordinates) and Open Babel `pybel` using the
`xyz` format (always available in minimal Open Babel builds).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem


def _mol_to_pdbqt_obabel(mol: Chem.Mol, out: Path) -> None:
    xyz = Chem.MolToXYZBlock(mol)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".xyz", delete=False, encoding="utf-8"
    ) as tf:
        tf.write(xyz)
        tf.flush()
        tmp = Path(tf.name)
    try:
        mols = list(pybel.readfile("xyz", str(tmp)))
        if not mols:
            raise RuntimeError("pybel could not read XYZ from RDKit output")
        mols[0].write("pdbqt", str(out), overwrite=True)
    finally:
        tmp.unlink(missing_ok=True)


def _main() -> int:
    a = sys.argv[1:]
    # -ipdb <file> -opdbqt -O <out>
    if len(a) >= 5 and a[0] == "-ipdb" and a[2] == "-opdbqt" and a[3] == "-O":
        inp = Path(a[1])
        out = Path(a[4])
        mol = Chem.MolFromPDBFile(str(inp), removeHs=False, sanitize=False)
        if mol is None:
            mol = Chem.MolFromPDBFile(str(inp), removeHs=False)
        if mol is None:
            return 1
        _mol_to_pdbqt_obabel(mol, out)
        return 0
    # -ismi <file> -opdbqt -O <out> --gen3d
    if len(a) >= 6 and a[0] == "-ismi" and a[2] == "-opdbqt" and a[3] == "-O" and a[-1] == "--gen3d":
        inp = Path(a[1])
        out = Path(a[4])
        text = inp.read_text(encoding="utf-8", errors="replace").splitlines()
        if not text:
            return 1
        smi = text[0].split()[0] if text[0].strip() else ""
        if not smi:
            return 1
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return 1
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
            return 1
        _mol_to_pdbqt_obabel(mol, out)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(_main())
