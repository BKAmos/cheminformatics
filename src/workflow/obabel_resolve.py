"""How we invoke Open Babel: bundled `tools/obabel_cli.py` or the real `obabel` binary."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
_SHIM = _REPO / "tools" / "obabel_cli.py"


def obabel_argv0() -> list[str] | None:
    """Argv prefix: either ``[sys.executable, shim]`` or ``[obabel]`` from PATH."""
    if _SHIM.is_file():
        return [sys.executable, str(_SHIM)]
    o = shutil.which("obabel")
    return [o] if o else None
