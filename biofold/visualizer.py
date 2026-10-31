"""3D Molecular Structure Visualization using py3Dmol and Streamlit components."""

from __future__ import annotations

import html
from typing import Any

import py3Dmol
import streamlit as st

try:
    import stmol  # type: ignore
    HAS_STMOL = True
except Exception:
    stmol = None
    HAS_STMOL = False


def build_3d_viewer(
    pdb_text: str,
    representation: str = "Cartoon",
    color_scheme: str = "By Confidence",
    spinning: bool = False,
    background_color: str = "#182230",
    surface_opacity: float = 0.65,
    width: int = 800,
    height: int = 560,
) -> Any:
    """Build and configure a py3Dmol viewer object."""
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(pdb_text, "pdb")
    viewer.setBackgroundColor(background_color)

    # Determine base representation keyword
    rep_key = "cartoon"
    base_style: dict[str, Any] = {}

    if representation == "Ribbon":
        rep_key = "cartoon"
        base_style = {"ribbon": True, "thickness": 0.4}
    elif representation == "Stick":
        rep_key = "stick"
        base_style = {"radius": 0.22}
    elif representation == "Sphere":
        rep_key = "sphere"
        base_style = {"scale": 0.32}
    elif representation == "Surface":
        rep_key = "cartoon"
        base_style = {"opacity": 0.35}

    # Apply Color Schemes
    if color_scheme == "By Confidence":
        # AlphaFold / ESMFold discrete tiers
        viewer.setStyle({}, {rep_key: {**base_style, "color": "#FF7D45"}})
        viewer.setStyle({"b": {"between": [50.0, 70.0]}}, {rep_key: {**base_style, "color": "#FFDB13"}})
        viewer.setStyle({"b": {"between": [70.0, 90.0]}}, {rep_key: {**base_style, "color": "#65CBF3"}})
        viewer.setStyle({"b": {"between": [90.0, 100.01]}}, {rep_key: {**base_style, "color": "#0053D6"}})
    elif color_scheme == "Secondary Structure":
        # ssPyMol: Helices red, Sheets yellow, Loops green
        viewer.setStyle({}, {rep_key: {**base_style, "color": "ssPyMol"}})
    elif color_scheme == "Residue Type":
        # Biochemical classification
        viewer.setStyle({}, {rep_key: {**base_style, "color": "#9CA3AF"}})  # Fallback
        viewer.setStyle(
            {"resn": ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO"]},
            {rep_key: {**base_style, "color": "#F59E0B"}},
        )
        viewer.setStyle(
            {"resn": ["SER", "THR", "ASN", "GLN", "CYS"]},
            {rep_key: {**base_style, "color": "#10B981"}},
        )
        viewer.setStyle(
            {"resn": ["LYS", "ARG", "HIS"]},
            {rep_key: {**base_style, "color": "#3B82F6"}},
        )
        viewer.setStyle(
            {"resn": ["ASP", "GLU"]},
            {rep_key: {**base_style, "color": "#EF4444"}},
        )
    elif color_scheme == "Rainbow":
        # N- to C-terminus spectrum
        viewer.setStyle({}, {rep_key: {**base_style, "color": "spectrum"}})
    else:
        viewer.setStyle({}, {rep_key: base_style})

    # Optional surface overlay
    if representation == "Surface":
        surf_style = {"opacity": surface_opacity}
        if color_scheme == "By Confidence":
            surf_style["colorscheme"] = {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}
        viewer.addSurface(py3Dmol.VDW, surf_style)

    viewer.zoomTo()
    if spinning:
        viewer.spin(True)
    else:
        viewer.spin(False)

    return viewer


def render_3d_viewer(
    viewer: Any,
    height: int = 580,
    use_stmol: bool = False,
) -> None:
    """Render the py3Dmol viewer into the Streamlit app.

    Uses stmol if requested and available, otherwise uses native
    streamlit.components.v1.html to ensure 100% stability.
    """
    if use_stmol and HAS_STMOL:
        try:
            stmol.showmol(viewer, height=height, width=None)
            return
        except Exception:
            pass  # Fallback to direct html

    raw_html = viewer._make_html()
    # Embed in iframe with responsive container
    st.components.v1.html(raw_html, height=height + 20, scrolling=False)

