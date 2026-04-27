"""Versioned Pydantic contracts between pipeline steps (plan § Robustness)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator


class ResourceLimits(BaseModel):
    max_parallel_http: int = Field(default=2, ge=1, le=64)
    max_parallel_docks: int = Field(default=1, ge=1, le=32)
    filter_chunk_size: int = Field(default=500, ge=1)
    tier_batch_size: int = Field(default=500, ge=1)
    pubchem_page_or_batch_limit: int = Field(default=10_000, ge=1)
    ligand_prep_batch_size: int = Field(default=200, ge=1)
    evolution_eval_batch_size: int = Field(default=50, ge=1)


class PubChemConfig(BaseModel):
    similarity_threshold: float = Field(default=0.85, ge=0.0, le=1.0)
    max_hits_per_query: int = Field(default=100, ge=1, le=10_000)
    base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class DockPoolCapConfig(BaseModel):
    mode: Literal["max_compounds", "top_n_similarity"] = "max_compounds"
    max_compounds_to_dock: int | None = Field(default=1000, ge=1)
    top_n_by_2d_score: int | None = Field(default=None, ge=1)


class FpocketConfig(BaseModel):
    """Pocket definition: manual box or fpocket prediction (receptor-only).

    Legacy YAML uses ``use_stub_box``; when ``mode`` is omitted, ``use_stub_box=False``
    selects fpocket mode.
    """

    mode: Literal["manual_box", "fpocket"] | None = None
    use_stub_box: bool = True
    stub_center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    stub_size: tuple[float, float, float] = (20.0, 20.0, 20.0)
    executable: str | None = None
    top_k_pockets: int = Field(default=1, ge=1, le=20)
    box_padding_angstrom: float = Field(default=4.0, ge=0.0)

    def resolved_mode(self) -> Literal["manual_box", "fpocket"]:
        if self.mode is not None:
            return self.mode
        return "manual_box" if self.use_stub_box else "fpocket"


class DockingConfig(BaseModel):
    """Physics-based (e.g. Vina/GNINA) vs ML-based (e.g. DiffDock) docking.

    * ``physics_only`` / ``ml_only`` — single family; artifacts match the classic layout when one
      family.
    * ``both`` — run both in tandem; writes per-family tables plus a merged
      ``docking_scores.parquet`` with a primary ``score`` chosen by ``ranking_when_both`` for
      evolution and summaries.
    """

    mode: Literal["physics_only", "ml_only", "both"] = "physics_only"
    physics_backend: Literal["mock", "vina", "gnina"] = "mock"
    ml_backend: Literal["mock", "diffdock"] = "mock"
    ranking_when_both: Literal["physics", "ml"] = "physics"
    vina_seed: int = Field(default=42)
    vina_exhaustiveness: int = Field(default=8, ge=1, le=32)


class MutationResidue(BaseModel):
    chain: str = Field(min_length=1, max_length=4)
    resseq: int = Field(ge=-999, le=9999)
    to_aa: str = Field(min_length=3, max_length=3, pattern=r"^[A-Z]{3}$")

    @field_validator("to_aa", mode="before")
    @classmethod
    def uppercase_aa(cls, v: object) -> object:
        if isinstance(v, str):
            return v.strip().upper()
        return v

    @field_validator("chain", mode="before")
    @classmethod
    def strip_chain(cls, v: object) -> object:
        if isinstance(v, str):
            return v.strip()
        return v


class MutationConfig(BaseModel):
    """Point mutations applied to the prepared receptor (optional OpenMM PDBFixer)."""

    residues: list[MutationResidue] = Field(default_factory=list)


class PairConfig(BaseModel):
    protein_id: str
    ligand_smiles_or_ids: list[str] = Field(default_factory=list)


class WorkflowConfig(BaseModel):
    library_csv: Path
    receptor_pdb: Path
    pairs: list[PairConfig] = Field(default_factory=list)
    resources: ResourceLimits = Field(default_factory=ResourceLimits)
    pubchem: PubChemConfig = Field(default_factory=PubChemConfig)
    dock_pool: DockPoolCapConfig = Field(default_factory=DockPoolCapConfig)
    fpocket: FpocketConfig = Field(default_factory=FpocketConfig)
    docking: DockingConfig = Field(default_factory=DockingConfig)
    mutation: MutationConfig = Field(default_factory=MutationConfig)
    dry_run: bool = False
    profile: Literal["default", "ci"] = "default"
    skip_evolution: bool = True
    skip_mutation: bool = True
    generate_plots: bool = False


class HitRecord(BaseModel):
    compound_id: str
    query_smiles: str
    cid: str | None = None
    hit_smiles: str
    tanimoto_similarity: float
    lineage_query_id: str


class PocketSpec(BaseModel):
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    source: str = "fpocket"
    raw_path: Path | None = None


class DockScoreRecord(BaseModel):
    compound_id: str
    score: float
    pose_path: str | None = None
    backend: str


class GroundTruthCaseResult(BaseModel):
    case_id: str
    passed: bool
    metrics: dict[str, Any] = Field(default_factory=dict)
    message: str = ""


class GroundTruthReport(BaseModel):
    cases: list[GroundTruthCaseResult]
    tool_versions: dict[str, str] = Field(default_factory=dict)
