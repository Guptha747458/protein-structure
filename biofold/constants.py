"""Constants and reference data for BioFold AI."""

from __future__ import annotations

# Standard 20 amino acids
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# 3-letter to 1-letter amino acid code mapping
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

# 1-letter to 3-letter mapping
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

# Amino acid full names
AMINO_ACID_NAMES = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartate",
    "C": "Cysteine",
    "Q": "Glutamine",
    "E": "Glutamate",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
}

# Biochemical classification of residues
RESIDUE_GROUPS = {
    "A": "Hydrophobic",
    "V": "Hydrophobic",
    "I": "Hydrophobic",
    "L": "Hydrophobic",
    "M": "Hydrophobic",
    "F": "Aromatic",
    "W": "Aromatic",
    "Y": "Aromatic",
    "S": "Polar Uncharged",
    "T": "Polar Uncharged",
    "N": "Polar Uncharged",
    "Q": "Polar Uncharged",
    "C": "Special / Disulfide",
    "G": "Special / Flexible",
    "P": "Special / Rigid",
    "K": "Positively Charged",
    "R": "Positively Charged",
    "H": "Positively Charged",
    "D": "Negatively Charged",
    "E": "Negatively Charged",
}

# Color palette for chemical groups
GROUP_COLORS = {
    "Hydrophobic": "#F59E0B",      # Amber
    "Aromatic": "#D97706",         # Dark Amber
    "Polar Uncharged": "#10B981",  # Emerald Green
    "Positively Charged": "#3B82F6", # Bright Blue
    "Negatively Charged": "#EF4444", # Vibrant Red
    "Special / Disulfide": "#8B5CF6", # Purple
    "Special / Flexible": "#9CA3AF",  # Gray
    "Special / Rigid": "#6B7280",     # Darker Gray
}

# Aliphatic index relative volume weights (Ikai, 1980)
ALIPHATIC_WEIGHTS = {
    "A": 1.0,
    "V": 2.9,
    "I": 3.9,
    "L": 3.9,
}

# Standard ESMFold public REST API URL
DEFAULT_ESMFOLD_URL = "https://api.esmatlas.com/v1/sequence/"

# RCSB PDB download template
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# pLDDT Confidence Tier Colors (standard AlphaFold / ESMFold convention)
PLDDT_TIERS = {
    "Very High": {"min": 90.0, "max": 100.0, "color": "#0053D6", "label": "pLDDT > 90 (Very High)"},
    "Confident": {"min": 70.0, "max": 90.0, "color": "#65CBF3", "label": "70 < pLDDT ≤ 90 (Confident)"},
    "Low": {"min": 50.0, "max": 70.0, "color": "#FFDB13", "label": "50 < pLDDT ≤ 70 (Low)"},
    "Very Low": {"min": 0.0, "max": 50.0, "color": "#FF7D45", "label": "pLDDT ≤ 50 (Very Low)"},
}

