"""BioFold AI: Protein Structure Prediction & 3D Interactive Visualization."""

from __future__ import annotations

from biofold.confidence import ConfidenceSummary, extract_plddt_from_pdb
from biofold.esm_client import ESMFoldError, predict_structure_esm
from biofold.presets import PRESETS, PresetProtein, get_preset_pdb
from biofold.rcsb_client import RCSBError, extract_sequence_from_pdb_text, fetch_rcsb_pdb
from biofold.reports import generate_json_report, generate_text_report
from biofold.sequence_analyzer import (
    ProteinProperties,
    analyze_protein_sequence,
    create_residue_distribution_dataframe,
    validate_and_clean_sequence,
)
from biofold.visualizer import build_3d_viewer, render_3d_viewer

__version__ = "1.0.0"
__all__ = [
    "ProteinProperties",
    "validate_and_clean_sequence",
    "analyze_protein_sequence",
    "create_residue_distribution_dataframe",
    "ConfidenceSummary",
    "extract_plddt_from_pdb",
    "predict_structure_esm",
    "ESMFoldError",
    "fetch_rcsb_pdb",
    "extract_sequence_from_pdb_text",
    "RCSBError",
    "PRESETS",
    "PresetProtein",
    "get_preset_pdb",
    "build_3d_viewer",
    "render_3d_viewer",
    "generate_json_report",
    "generate_text_report",
]

