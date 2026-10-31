"""Structural summary report generator in JSON and plain-text formats."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any

from biofold.confidence import ConfidenceSummary
from biofold.sequence_analyzer import ProteinProperties


def generate_json_report(
    source_label: str,
    properties: ProteinProperties,
    confidence: ConfidenceSummary,
    additional_metadata: dict[str, Any] | None = None,
) -> str:
    """Generate a structured, machine-readable JSON summary report."""
    report = {
        "generator": "BioFold AI - Protein Structure Prediction & 3D Interactive Visualization",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": source_label,
        "metadata": additional_metadata or {},
        "biophysical_properties": properties.to_dict(),
        "structure_confidence": confidence.to_dict(),
    }
    return json.dumps(report, indent=2)


def generate_text_report(
    source_label: str,
    properties: ProteinProperties,
    confidence: ConfidenceSummary,
    additional_metadata: dict[str, Any] | None = None,
) -> str:
    """Generate an executive, human-readable plain text / markdown summary report."""
    now_str = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = [
        "=" * 72,
        "  BIOFOLD AI: PROTEIN STRUCTURE PREDICTION & ANALYSIS REPORT",
        "=" * 72,
        f"Generated At      : {now_str}",
        f"Structure Source  : {source_label}",
    ]

    if additional_metadata:
        for k, v in additional_metadata.items():
            lines.append(f"{k.replace('_', ' ').title():<18}: {v}")

    lines.extend([
        "-" * 72,
        "1. PRIMARY SEQUENCE & COMPOSITION",
        "-" * 72,
        f"Length            : {properties.length} amino acids",
        f"Molecular Weight  : {properties.molecular_weight_kda:.2f} kDa",
        f"Isoelectric Point : pH {properties.isoelectric_point:.2f}",
        f"Net Charge (pH 7) : {properties.charge_at_ph_7:+.2f}",
        f"Instability Index : {properties.instability_index:.2f} ({'Stable' if properties.is_stable else 'Unstable'})",
        f"Aliphatic Index   : {properties.aliphatic_index:.2f}",
        f"GRAVY Hydropathy  : {properties.gravy:+.3f}",
        "",
        "Secondary Structure Fractions (Estimated):",
        f"  - Alpha Helix (α) : {properties.secondary_structure.get('Helix (α)', 0.0):.1f}%",
        f"  - Turn            : {properties.secondary_structure.get('Turn', 0.0):.1f}%",
        f"  - Beta Sheet (β)  : {properties.secondary_structure.get('Sheet (β)', 0.0):.1f}%",
        "",
        "Sequence (FASTA):",
    ])

    # Wrap sequence at 60 chars
    seq = properties.sequence
    for i in range(0, len(seq), 60):
        lines.append(f"  {seq[i:i+60]}")

    lines.extend([
        "",
        "-" * 72,
        "2. 3D STRUCTURE CONFIDENCE ASSESSMENT (pLDDT / B-Factor)",
        "-" * 72,
        f"Total Atoms       : {confidence.total_atoms}",
        f"Total Residues    : {confidence.total_residues}",
        f"Mean pLDDT        : {confidence.mean_plddt:.2f}",
        f"Median pLDDT      : {confidence.median_plddt:.2f}",
        f"Min / Max pLDDT   : {confidence.min_plddt:.2f} / {confidence.max_plddt:.2f}",
        f"Std Dev           : {confidence.std_plddt:.2f}",
        "",
        "Confidence Class Breakdown:",
        f"  - Very High (pLDDT > 90) : {confidence.class_counts.get('Very High', 0):>3} residues ({confidence.class_percentages.get('Very High', 0.0):.1f}%)",
        f"  - Confident (70-90)      : {confidence.class_counts.get('Confident', 0):>3} residues ({confidence.class_percentages.get('Confident', 0.0):.1f}%)",
        f"  - Low (50-70)            : {confidence.class_counts.get('Low', 0):>3} residues ({confidence.class_percentages.get('Low', 0.0):.1f}%)",
        f"  - Very Low (≤ 50)        : {confidence.class_counts.get('Very Low', 0):>3} residues ({confidence.class_percentages.get('Very Low', 0.0):.1f}%)",
        "",
        "=" * 72,
        "End of BioFold AI Report",
        "=" * 72,
    ])

    return "\n".join(lines)

