"""How we invoke Open Babel: bundled `tools/obabel_cli.py` or the real `obabel` binary."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
_SHIM = _REPO / "tools" / "obabel_cli.py"


def _shim_imports_ok() -> bool:
    """`tools/obabel_cli.py` needs RDKit + Open Babel pybel (e.g. ``pip install -e '.[chem]'``)."""
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    try:
        from openbabel import pybel  # noqa: F401
    except ImportError:
        return False
    return True


def obabel_argv0() -> list[str] | None:
    """Argv prefix: either ``[sys.executable, shim]`` or ``[obabel]`` from PATH.

    The bundled shim is not used in minimal envs (e.g. GitHub CI with only
    ``.[dev]``) so optional ligand PDBQT precompute is skipped, matching
    environments without a real ``obabel`` binary.
    """
    if _SHIM.is_file() and _shim_imports_ok():
        return [sys.executable, str(_SHIM)]
    o = shutil.which("obabel")
    return [o] if o else None
