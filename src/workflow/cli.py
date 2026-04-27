"""CLI: workflow run | plot | validate."""

from __future__ import annotations

from pathlib import Path

import typer

from workflow.config_load import load_workflow_config
from workflow.pipeline import run_pipeline

app = typer.Typer(no_args_is_help=True)


@app.command()
def run(
    config: Path = typer.Option(..., "--config", exists=True, dir_okay=False),
    run_dir: Path | None = typer.Option(None, "--run-dir"),
    dry_run: bool = typer.Option(False, "--dry-run"),
    no_langgraph: bool = typer.Option(False, "--no-langgraph"),
) -> None:
    cfg = load_workflow_config(config)
    if dry_run:
        cfg = cfg.model_copy(update={"dry_run": True})
    rd = run_dir or Path(".runs") / "latest"
    rd.mkdir(parents=True, exist_ok=True)
    run_pipeline(cfg, rd, use_langgraph=not no_langgraph)
    typer.echo(f"Run complete. Artifacts: {rd.resolve()}")


@app.command()
def plot(
    run_dir: Path = typer.Option(..., "--run", exists=True, file_okay=False),
) -> None:
    try:
        from workflow.reporting import plots as plotmod
    except ImportError as e:
        typer.echo(f"Plotting requires matplotlib and pandas: {e}", err=True)
        raise typer.Exit(1) from e
    plotmod.generate_all(run_dir)
    typer.echo(f"Plots written to {run_dir / 'plots'}")


@app.command("validate")
def validate_cmd(
    suite: str = typer.Option("ground_truth", "--suite"),
    run_dir: Path | None = typer.Option(
        None,
        "--run-dir",
        exists=True,
        file_okay=False,
        help="Optional workflow run directory for artifact-based benchmark cases.",
    ),
) -> None:
    if suite != "ground_truth":
        typer.echo("Only ground_truth suite implemented", err=True)
        raise typer.Exit(1)
    from workflow.validation import ground_truth as gt

    out = Path(".runs") / "validation"
    out.mkdir(parents=True, exist_ok=True)
    report_path = gt.run_ground_truth_suite(out, run_dir=run_dir)
    typer.echo(f"Report: {report_path}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
