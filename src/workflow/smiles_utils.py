"""Canonical SMILES when RDKit is available; fallback normalization."""

from __future__ import annotations


def canonical_smiles(smiles: str) -> str:
    try:
        from rdkit import Chem
        from rdkit.Chem import MolToSmiles

        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return smiles.strip()
        return MolToSmiles(m)
    except ImportError:
        return smiles.strip()
