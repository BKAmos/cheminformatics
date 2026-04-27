# Chem workflow

Python **3.11+** pipeline: **PubChem 2D similarity → Tier A/B filters → dock pool cap → receptor + ligand prep → fpocket (stub) → WT docking → mutant docking (stub) → evolution (stub) → evolved refilter**, with **LangGraph** orchestration, **Pydantic** contracts, run-directory artifacts, optional **matplotlib** plots, and a **ground-truth** report stub.

> Full herbicide science, ChEMBL/ECOTOX integration, real GNINA/Vina/fpocket/OpenMM, and production evolution are **not** finished here—this repo is a **runnable scaffold** aligned with `multi-agent_chem_workflow_85a47944.plan.md`.

## Quick start

```bash
python -m pip install -e ".[dev]"
workflow run --config configs/ci.yaml --run-dir .runs/demo
workflow plot --run .runs/demo
workflow validate --suite ground_truth
python -m pytest -q
```

- **`configs/ci.yaml`**: small limits, suitable for CI and laptops.
- **`configs/example.yaml`**: larger defaults; uses **mock docking** by default.

### Binding site (manual box vs fpocket)

- **`fpocket.mode: manual_box`** (default when `use_stub_box: true`): uses `stub_center` / `stub_size` and writes `structure/pockets_ranked.parquet` with a single pocket for downstream steps.
- **`fpocket.mode: fpocket`** (or `use_stub_box: false` without an explicit mode): runs the **fpocket** binary on `receptor_prepared.pdb`, parses cavities, writes `structure/pocket_specs/pocket_*.json`, copies the top site to `pocket_spec.json`, and ranks pockets in `pockets_ranked.parquet`.
- **`fpocket.top_k_pockets`**: dock against the top *K* ranked pockets (scores merged by best affinity per compound). Not supported together with `docking.mode: both` (falls back to `K=1` with a warning).

### Docking (physics vs ML)

`WorkflowConfig.docking` selects **physics-based** backends (`mock`, `vina`, `gnina`) and/or **ML-based** backends (`mock`, `diffdock` stub). Set `mode` to `physics_only`, `ml_only`, or `both`. In `both`, the pipeline runs each family and writes `docking_scores_physics.parquet` / `docking_scores_ml.parquet` plus a merged `docking_scores.parquet` whose primary `score` follows `ranking_when_both` (`physics` or `ml`) for evolution and summaries. Pose files live under `poses/wt/physics/` and `poses/wt/ml/` when `mode` is `both`; otherwise they stay directly under `poses/wt/`.

### Vina backend notes

- Set `docking.physics_backend: vina` in config.
- Requires `vina` on PATH.
- If receptor/ligands are not already `.pdbqt`, backend uses `obabel` to convert from PDB/SMILES.
- Vina logs are written per ligand to `poses/.../logs/*.vina.log`; best affinity is parsed from the Vina table.
- Fast path: if dock pool rows include `ligand_pdbqt_path`, Vina uses those directly and skips ligand SMILES conversion.
- `run_ligand_prep` now writes `ligands/ligands_index.parquet` and enriches `filters/dock_pool.parquet` with `ligand_pdbqt_path` when Open Babel is available.
- AutoDock Vina runs with fixed **`--seed`** and **`--exhaustiveness`** from `WorkflowConfig.docking` for more reproducible scores.
- After WT/mutant docking, **`summary/candidates_ranked.parquet`** combines potency + selectivity (`delta_score`) with a **`composite_rank_score`** and text **`rationale`**.
- Optional real point mutations: configure `mutation.residues` and set `skip_mutation: false`, then install extras with `pip install -e ".[mutate]"` (OpenMM + PDBFixer).

## Layout

- `src/workflow/` — config, contracts, LangGraph `pipeline.py`, steps, `backends/mock.py`, `reporting/`, `validation/`.
- `configs/` — YAML inputs (paths relative to each YAML file).
- `data/` — toy library + PDB for smoke tests.
- `runs/` or `.runs/` — gitignored run outputs (Parquet, JSON, logs).

## Resource defaults

`WorkflowConfig.resources` caps **HTTP** and **parallel docks** (mock backend ignores parallelism but honors the API). Prefer raising **wall time** over **RAM** by keeping `max_parallel_docks: 1` on shared machines.

## License

MIT — see `LICENSE`. Third-party binaries are **not** vendored; see `THIRD_PARTY_NOTICES.md`.
