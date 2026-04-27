"""Safe subprocess wrapper (plan § Robustness): no shell=True, cwd under run dir."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Sequence


def run_cmd(
    argv: Sequence[str],
    *,
    cwd: Path | None = None,
    timeout_s: float | None = None,
    capture: bool = True,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv),
        cwd=cwd,
        timeout=timeout_s,
        capture_output=capture,
        text=True,
        shell=False,
        check=False,
    )
