"""Structured step logging (JSON lines)."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any


def log_step(
    paths: dict[str, Path],
    step: str,
    duration_s: float,
    input_count: int | None = None,
    output_count: int | None = None,
    exit_code: int = 0,
    extra: dict[str, Any] | None = None,
) -> None:
    line = {
        "step": step,
        "duration_s": round(duration_s, 4),
        "input_count": input_count,
        "output_count": output_count,
        "exit_code": exit_code,
        **(extra or {}),
    }
    logf = paths["logs"] / "steps.jsonl"
    with logf.open("a", encoding="utf-8") as f:
        f.write(json.dumps(line) + "\n")


def timed_step(paths: dict[str, Path], step: str):
    def decorator(fn):
        def wrapper(*args, **kwargs):
            t0 = time.perf_counter()
            try:
                out = fn(*args, **kwargs)
                log_step(
                    paths,
                    step,
                    time.perf_counter() - t0,
                    exit_code=0,
                )
                return out
            except Exception as e:
                log_step(
                    paths,
                    step,
                    time.perf_counter() - t0,
                    exit_code=1,
                    extra={"error": str(e)},
                )
                raise

        return wrapper

    return decorator
