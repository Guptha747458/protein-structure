"""RCSB PDB fetch client and structure-to-sequence extraction."""

from __future__ import annotations

import re
from io import StringIO
from pathlib import Path
import requests
import streamlit as st
from Bio.PDB import PDBParser

from biofold.constants import RCSB_PDB_URL, THREE_TO_ONE
from biofold.sequence_analyzer import validate_and_clean_sequence

LOCAL_DATA_DIR = Path(__file__).parent / "data"

PDB_ID_REGEX = re.compile(r"^[0-9A-Za-z]{4}$")


class RCSBError(Exception):
    """Exception raised for RCSB PDB retrieval or parsing errors."""
    pass


@st.cache_data(show_spinner=False, ttl=86400)
def fetch_rcsb_pdb(pdb_id: str, timeout: int = 25) -> str:
    """Download a verified PDB structure from RCSB Protein Data Bank.

    Includes automatic fallback to local bundled PDB data if offline.
    """
    clean_id = pdb_id.strip().upper()
    if not PDB_ID_REGEX.match(clean_id):
        raise RCSBError(
            f"Invalid PDB ID '{clean_id}'. PDB IDs must be exactly 4 alphanumeric characters (e.g., '1UBQ', '1CRN')."
        )

    # Check local pre-cached data first if network fails
    local_mapping = {
        "1UBQ": "ubiquitin_clean.pdb",
        "1CRN": "crambin_clean.pdb",
        "1L2Y": "trp_cage_clean.pdb",
        "4INS": "insulin_clean.pdb",
    }

    url = RCSB_PDB_URL.format(pdb_id=clean_id.lower())
    try:
        response = requests.get(url, timeout=(10, timeout))
        if response.status_code == 200:
            return response.text
        elif response.status_code == 404:
            raise RCSBError(f"PDB ID '{clean_id}' was not found in the RCSB database.")
        else:
            raise RCSBError(f"RCSB PDB returned HTTP status {response.status_code}.")
    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as err:
        # Fallback to local copy if available
        if clean_id in local_mapping:
            local_file = LOCAL_DATA_DIR / local_mapping[clean_id]
            if local_file.exists():
                return local_file.read_text(encoding="utf-8")
        raise RCSBError(
            f"Failed to connect to RCSB PDB server for entry '{clean_id}'. "
            "Please check network connection or try a bundled preset sequence."
        )
    except requests.RequestException as err:
        raise RCSBError(f"Network error while fetching PDB '{clean_id}': {err}")


def extract_sequence_from_pdb_text(pdb_text: str, chain_id: str | None = None) -> tuple[str, list[str]]:
    """Parse PDB ATOM coordinates and extract the canonical amino acid sequence.

    Returns:
        tuple of (primary_sequence, available_chain_ids)
    """
    if not pdb_text or "ATOM" not in pdb_text:
        raise RCSBError("PDB text does not contain ATOM records.")

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("reference", StringIO(pdb_text))
    except Exception as err:
        raise RCSBError(f"Biopython failed to parse PDB structure: {err}")

    models = list(structure.get_models())
    if not models:
        raise RCSBError("Structure contains no models.")

    model = models[0]
    chains = list(model.get_chains())
    if not chains:
        raise RCSBError("Model contains no peptide chains.")

    available_chains = [c.id.strip() or "A" for c in chains]

    # Select requested chain or first chain
    target_chain = None
    if chain_id and any(c.id == chain_id for c in chains):
        target_chain = next(c for c in chains if c.id == chain_id)
    else:
        target_chain = chains[0]

    residue_letters: list[str] = []
    for residue in target_chain.get_residues():
        # Standard residues have blank hetero-flag (' ')
        if residue.id[0] == " ":
            resname = residue.resname.strip().upper()
            if resname in THREE_TO_ONE:
                residue_letters.append(THREE_TO_ONE[resname])

    if not residue_letters:
        raise RCSBError(
            f"Chain '{target_chain.id}' contains no standard amino acid residues."
        )

    full_seq = "".join(residue_letters)
    cleaned_seq, _ = validate_and_clean_sequence(full_seq)
    return cleaned_seq, available_chains

