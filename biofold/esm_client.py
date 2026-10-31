"""ESMFold API client for single-sequence de novo protein structure prediction."""

from __future__ import annotations

import logging
import requests
import streamlit as st

from biofold.constants import DEFAULT_ESMFOLD_URL

logger = logging.getLogger(__name__)


class ESMFoldError(Exception):
    """Custom exception raised when structure prediction fails."""
    pass


@st.cache_data(show_spinner=False, ttl=3600)
def predict_structure_esm(
    sequence: str,
    api_url: str = DEFAULT_ESMFOLD_URL,
    api_token: str | None = None,
    timeout: int = 60,
) -> str:
    """Submit a protein amino acid sequence to the ESMFold REST API.

    Args:
        sequence: Validated uppercase single-letter amino acid sequence.
        api_url: ESMFold API endpoint URL.
        api_token: Optional Bearer token for authenticated ESM endpoints.
        timeout: Network request timeout in seconds.

    Returns:
        PDB file contents as string.

    Raises:
        ESMFoldError: If the API request times out or returns an error.
    """
    clean_seq = sequence.strip().upper()
    if len(clean_seq) > 400:
        logger.warning(
            f"Sequence length is {len(clean_seq)} aa. Sequences longer than 400 aa "
            "may experience extended latency or timeouts on public ESMFold endpoints."
        )

    headers = {"Content-Type": "text/plain"}
    if api_token:
        headers["Authorization"] = f"Bearer {api_token.strip()}"

    try:
        response = requests.post(
            api_url.strip(),
            data=clean_seq,
            headers=headers,
            timeout=(15, timeout),
        )
    except requests.exceptions.Timeout:
        raise ESMFoldError(
            f"Structure prediction timed out after {timeout} seconds. "
            "The ESMFold server might be experiencing high load. "
            "Try again shortly or reduce sequence length."
        )
    except requests.exceptions.ConnectionError as err:
        raise ESMFoldError(
            f"Unable to connect to ESMFold endpoint ({api_url}). "
            "Please check your internet connection or network firewall.\n"
            f"Details: {err}"
        )
    except requests.RequestException as err:
        raise ESMFoldError(f"Network error during ESMFold API call: {err}")

    if response.status_code == 200:
        pdb_text = response.text.strip()
        if not pdb_text or "ATOM" not in pdb_text:
            raise ESMFoldError(
                "ESMFold returned a 200 OK status, but the response contained no PDB ATOM coordinates."
            )
        return pdb_text

    if response.status_code == 403:
        raise ESMFoldError(
            "ESMFold API returned HTTP 403 (Forbidden / Missing Token). "
            "The public Meta ESMAtlas endpoint requires an API token. "
            "You can enter an API Token under 'Advanced Settings' in the sidebar, "
            "or explore verified ground-truth structures using 'Preset Demo' or 'PDB ID Lookup'."
        )
    elif response.status_code == 429:
        raise ESMFoldError(
            "ESMFold API rate limit reached (HTTP 429). Please wait a moment before trying again."
        )
    elif response.status_code == 400:
        raise ESMFoldError(
            f"ESMFold rejected the request (HTTP 400): {response.text[:200]}"
        )
    elif response.status_code >= 500:
        raise ESMFoldError(
            f"ESMFold server error (HTTP {response.status_code}). The ESMAtlas service may be temporarily down."
        )
    else:
        raise ESMFoldError(
            f"ESMFold API returned HTTP {response.status_code}: {response.text[:200]}"
        )

