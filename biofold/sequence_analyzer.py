"""Sequence validation, FASTA parsing, and biophysical property analytics."""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Any

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from biofold.constants import (
    ALIPHATIC_WEIGHTS,
    AMINO_ACID_NAMES,
    GROUP_COLORS,
    RESIDUE_GROUPS,
    STANDARD_AMINO_ACIDS,
)


@dataclass(frozen=True)
class ProteinProperties:
    sequence: str
    header: str | None
    length: int
    molecular_weight_kda: float
    isoelectric_point: float
    instability_index: float
    is_stable: bool
    aliphatic_index: float
    gravy: float
    secondary_structure: dict[str, float]  # Percentage: Helix, Turn, Sheet
    amino_acid_counts: dict[str, int]
    amino_acid_percentages: dict[str, float]
    charge_at_ph_7: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def validate_and_clean_sequence(raw: str) -> tuple[str, str | None]:
    """Validate and clean FASTA or raw amino acid sequence string.

    Returns:
        tuple of (cleaned_sequence, fasta_header_or_None)

    Raises:
        ValueError: If sequence is empty or contains non-standard amino acids.
    """
    if not raw or not raw.strip():
        raise ValueError("Sequence input is empty. Please paste a FASTA or amino acid sequence.")

    header: str | None = None
    lines = raw.strip().splitlines()
    seq_lines: list[str] = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if header is None:
                header = stripped[1:].strip()
            continue
        seq_lines.append(stripped)

    combined = "".join(seq_lines).upper()
    cleaned = re.sub(r"[\s\d\-_.*]", "", combined)

    if not cleaned:
        raise ValueError("No valid amino acid characters found in the input.")

    invalid_chars = sorted(set(cleaned) - STANDARD_AMINO_ACIDS)
    if invalid_chars:
        invalid_list = ", ".join(f"'{c}'" for c in invalid_chars)
        raise ValueError(
            f"Invalid amino acid character(s) detected: {invalid_list}. "
            "BioFold AI supports the 20 standard single-letter amino acid codes (ACDEFGHIKLMNPQRSTVWY)."
        )

    return cleaned, header


def calculate_aliphatic_index(sequence: str) -> float:
    """Calculate the Aliphatic Index (Ikai, 1980) of a globular protein."""
    if not sequence:
        return 0.0
    total = len(sequence)
    counts = {aa: sequence.count(aa) for aa in ("A", "V", "I", "L")}
    score = (
        counts["A"] * ALIPHATIC_WEIGHTS["A"]
        + counts["V"] * ALIPHATIC_WEIGHTS["V"]
        + (counts["I"] + counts["L"]) * ALIPHATIC_WEIGHTS["I"]
    )
    return round((100.0 * score) / total, 2)


def analyze_protein_sequence(sequence: str, header: str | None = None) -> ProteinProperties:
    """Perform biophysical and sequence analysis using Bio.SeqUtils.ProtParam."""
    anal = ProteinAnalysis(sequence)
    counts = anal.count_amino_acids()
    length = len(sequence)

    percentages = {aa: round((cnt / length) * 100.0, 2) for aa, cnt in counts.items()}

    # Secondary structure fraction estimates (Helix, Turn, Sheet)
    helix, turn, sheet = anal.secondary_structure_fraction()
    sec_struct = {
        "Helix (α)": round(helix * 100.0, 2),
        "Turn": round(turn * 100.0, 2),
        "Sheet (β)": round(sheet * 100.0, 2),
    }

    instability = anal.instability_index()
    is_stable = instability < 40.0

    return ProteinProperties(
        sequence=sequence,
        header=header,
        length=length,
        molecular_weight_kda=round(anal.molecular_weight() / 1000.0, 2),
        isoelectric_point=round(anal.isoelectric_point(), 2),
        instability_index=round(instability, 2),
        is_stable=is_stable,
        aliphatic_index=calculate_aliphatic_index(sequence),
        gravy=round(anal.gravy(), 3),
        secondary_structure=sec_struct,
        amino_acid_counts={k: counts[k] for k in sorted(counts)},
        amino_acid_percentages={k: percentages[k] for k in sorted(percentages)},
        charge_at_ph_7=round(anal.charge_at_pH(7.0), 2),
    )


def create_residue_distribution_dataframe(counts: dict[str, int], length: int) -> pd.DataFrame:
    """Build a rich DataFrame suitable for residue frequency charts."""
    records = []
    for aa in sorted(STANDARD_AMINO_ACIDS):
        cnt = counts.get(aa, 0)
        pct = round((cnt / length) * 100.0, 2) if length > 0 else 0.0
        group = RESIDUE_GROUPS.get(aa, "Other")
        color = GROUP_COLORS.get(group, "#6B7280")
        name = AMINO_ACID_NAMES.get(aa, aa)
        records.append(
            {
                "Residue": aa,
                "Name": name,
                "Count": cnt,
                "Frequency (%)": pct,
                "Chemical Group": group,
                "Color": color,
            }
        )
    return pd.DataFrame(records)

