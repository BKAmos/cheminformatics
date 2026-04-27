"""Safe subprocess wrapper (plan § Robustness): no shell=True, cwd under run dir."""

from __future__ import annotations

import os
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path


def run_cmd(
    argv: Sequence[str],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    timeout_s: float | None = None,
    capture: bool = True,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv),
        cwd=cwd,
        env=os.environ if env is None else dict(env),
        timeout=timeout_s,
        capture_output=capture,
        text=True,
        shell=False,
        check=False,
    )
