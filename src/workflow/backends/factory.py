"""Instantiate docking backends from config names."""

from __future__ import annotations

from workflow.backends.base import DockingBackend


def get_physics_backend(name: str) -> DockingBackend:
    if name == "mock":
        from workflow.backends.mock import MockDockingBackend

        return MockDockingBackend()
    if name == "vina":
        from workflow.backends.vina import VinaDockingBackend

        return VinaDockingBackend()
    if name == "gnina":
        from workflow.backends.gnina import GninaDockingBackend

        return GninaDockingBackend()
    raise ValueError(f"Unknown physics docking backend: {name!r}")


def get_ml_backend(name: str) -> DockingBackend:
    if name == "mock":
        from workflow.backends.mock import MockMLDockingBackend

        return MockMLDockingBackend()
    if name == "diffdock":
        from workflow.backends.diffdock import DiffDockBackend

        return DiffDockBackend()
    raise ValueError(f"Unknown ML docking backend: {name!r}")
