"""3D Molecular Structure Visualization using py3Dmol and Streamlit components."""

from __future__ import annotations

import re
from typing import Any

import py3Dmol
import streamlit as st

from biofold.confidence import ConfidenceSummary, extract_plddt_from_pdb

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
    confidence_summary: ConfidenceSummary | None = None,
) -> Any:
    """Build and configure a py3Dmol viewer object with responsive layout and centered model.

    Initializes py3Dmol with data directly so camera center and bounding box are computed
    from the protein coordinates, ensuring natural, centered 3D mouse rotation and dragging.
    """
    # 1. Initialize viewer with data so 3Dmol centers on the model centroid immediately
    viewer = py3Dmol.view(
        data=pdb_text,
        format="pdb",
        width="100%",
        height="100%",
    )
    viewer.setBackgroundColor(background_color)

    # 2. Configure base representation style
    rep_key = "cartoon"
    base_style: dict[str, Any] = {}

    if representation == "Ribbon":
        rep_key = "cartoon"
        base_style = {"ribbon": True, "thickness": 0.45}
    elif representation == "Stick":
        rep_key = "stick"
        base_style = {"radius": 0.22}
    elif representation == "Sphere":
        rep_key = "sphere"
        base_style = {"scale": 0.32}
    elif representation == "Surface":
        rep_key = "cartoon"
        base_style = {"opacity": 0.4}

    # 3. Apply Color Schemes
    if color_scheme == "By Confidence":
        # Obtain residue-level confidence
        conf = confidence_summary if confidence_summary else extract_plddt_from_pdb(pdb_text)

        if conf and conf.residues:
            # Set default fallback color
            viewer.setStyle({}, {rep_key: {**base_style, "color": "#65CBF3"}})
            # Apply exact AlphaFold / ESMFold discrete tier colors to every residue
            for res in conf.residues:
                viewer.setStyle(
                    {"resi": res.position},
                    {rep_key: {**base_style, "color": res.color}},
                )
        else:
            # Fallback to B-factor spectrum
            viewer.setStyle(
                {},
                {
                    rep_key: {
                        **base_style,
                        "colorscheme": {
                            "prop": "b",
                            "gradient": "roygb",
                            "min": 50,
                            "max": 90,
                        },
                    }
                },
            )

    elif color_scheme == "Secondary Structure":
        # PyMOL secondary structure coloring: Helices red, Sheets yellow, Loops green
        viewer.setStyle({}, {rep_key: {**base_style, "color": "ssPyMol"}})

    elif color_scheme == "Residue Type":
        # Biochemical classification
        viewer.setStyle({}, {rep_key: {**base_style, "color": "#9CA3AF"}})
        # Hydrophobic
        viewer.setStyle(
            {"resn": ["ALA", "VAL", "ILE", "LEU", "MET", "PRO"]},
            {rep_key: {**base_style, "color": "#F59E0B"}},
        )
        # Aromatic
        viewer.setStyle(
            {"resn": ["PHE", "TRP", "TYR"]},
            {rep_key: {**base_style, "color": "#D97706"}},
        )
        # Polar uncharged
        viewer.setStyle(
            {"resn": ["SER", "THR", "ASN", "GLN", "CYS"]},
            {rep_key: {**base_style, "color": "#10B981"}},
        )
        # Positively charged
        viewer.setStyle(
            {"resn": ["LYS", "ARG", "HIS"]},
            {rep_key: {**base_style, "color": "#3B82F6"}},
        )
        # Negatively charged
        viewer.setStyle(
            {"resn": ["ASP", "GLU"]},
            {rep_key: {**base_style, "color": "#EF4444"}},
        )

    elif color_scheme == "Rainbow":
        # N- to C-terminus spectrum
        viewer.setStyle({}, {rep_key: {**base_style, "color": "spectrum"}})
    else:
        viewer.setStyle({}, {rep_key: base_style})

    # 4. Optional surface overlay
    if representation == "Surface":
        surf_style: dict[str, Any] = {"opacity": surface_opacity}
        if color_scheme == "By Confidence":
            surf_style["colorscheme"] = {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}
        viewer.addSurface(py3Dmol.VDW, surf_style)

    # 5. Center camera precisely on the protein and zoom to fit viewport
    viewer.zoomTo()
    if spinning:
        viewer.spin(True)
    else:
        viewer.spin(False)

    return viewer


def render_3d_viewer(
    viewer: Any,
    height: int = 560,
    background_color: str = "#182230",
) -> None:
    """Render the py3Dmol viewer into Streamlit with interactive HUD controls.

    Injects zero-margin iframe styling, auto-resize listeners, and on-canvas
    HUD buttons (Reset Camera, Zoom In, Zoom Out, Spin toggle) to prevent dragging issues.
    """
    raw_html = viewer._make_html()

    # Capture the viewer instance variable from py3Dmol generated script
    match = re.search(r"(viewer_\d+)\s*=\s*\$3Dmol\.createViewer", raw_html)
    viewer_var = match.group(1) if match else "viewer"

    # Wrap inside full-bleed, zero-margin HTML container with interactive HUD overlay
    enhanced_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            html, body {{
                margin: 0 !important;
                padding: 0 !important;
                width: 100% !important;
                height: 100% !important;
                overflow: hidden !important;
                background-color: {background_color} !important;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            }}
            #viewer-wrapper {{
                position: relative;
                width: 100%;
                height: {height}px;
                border-radius: 12px;
                overflow: hidden;
                box-shadow: 0 4px 20px rgba(0, 0, 0, 0.25);
            }}
            /* HUD Toolbar Overlay */
            .hud-bar {{
                position: absolute;
                top: 12px;
                right: 12px;
                z-index: 100;
                display: flex;
                gap: 6px;
                background: rgba(15, 23, 42, 0.75);
                backdrop-filter: blur(8px);
                border: 1px solid rgba(255, 255, 255, 0.15);
                padding: 4px 6px;
                border-radius: 8px;
            }}
            .hud-btn {{
                background: rgba(255, 255, 255, 0.1);
                color: #F8FAFC;
                border: 1px solid rgba(255, 255, 255, 0.2);
                border-radius: 6px;
                padding: 5px 10px;
                font-size: 11px;
                font-weight: 600;
                cursor: pointer;
                transition: all 0.15s ease-in-out;
                display: flex;
                align-items: center;
                gap: 4px;
            }}
            .hud-btn:hover {{
                background: rgba(59, 130, 246, 0.8);
                border-color: #3B82F6;
                transform: translateY(-1px);
            }}
            /* Mouse Gesture Hint */
            .mouse-hint {{
                position: absolute;
                bottom: 12px;
                left: 12px;
                z-index: 100;
                background: rgba(15, 23, 42, 0.75);
                backdrop-filter: blur(8px);
                color: #CBD5E1;
                border: 1px solid rgba(255, 255, 255, 0.12);
                border-radius: 6px;
                padding: 4px 10px;
                font-size: 10.5px;
                pointer-events: none;
                user-select: none;
            }}
        </style>
    </head>
    <body>
        <div id="viewer-wrapper">
            <!-- HUD Control Bar -->
            <div class="hud-bar">
                <button class="hud-btn" onclick="resetCamera()" title="Re-center and fit structure">🎯 Reset</button>
                <button class="hud-btn" onclick="zoomIn()" title="Zoom In">🔍+</button>
                <button class="hud-btn" onclick="zoomOut()" title="Zoom Out">🔍-</button>
                <button class="hud-btn" onclick="toggleSpin()" title="Toggle Rotation">🔄 Spin</button>
            </div>

            <!-- Gesture Instructions -->
            <div class="mouse-hint">
                🖱️ Left-drag: <b>Rotate</b> &bull; Right-drag: <b>Pan</b> &bull; Scroll: <b>Zoom</b>
            </div>

            {raw_html}
        </div>

        <script>
            var isSpinning = false;

            function getActiveViewer() {{
                try {{
                    if (typeof {viewer_var} !== 'undefined' && {viewer_var} !== null) {{
                        return {viewer_var};
                    }}
                }} catch(e) {{}}
                return null;
            }}

            function resetCamera() {{
                var v = getActiveViewer();
                if (v) {{
                    v.zoomTo();
                    v.render();
                }}
            }}

            function zoomIn() {{
                var v = getActiveViewer();
                if (v) {{
                    v.zoom(1.25);
                    v.render();
                }}
            }}

            function zoomOut() {{
                var v = getActiveViewer();
                if (v) {{
                    v.zoom(0.8);
                    v.render();
                }}
            }}

            function toggleSpin() {{
                var v = getActiveViewer();
                if (v) {{
                    isSpinning = !isSpinning;
                    v.spin(isSpinning);
                }}
            }}

            // Resize WebGL canvas dynamically when iframe dimension adjusts
            window.addEventListener('resize', function() {{
                var v = getActiveViewer();
                if (v) {{
                    v.resize();
                    v.render();
                }}
            }});
        </script>
    </body>
    </html>
    """

    st.components.v1.html(enhanced_html, height=height + 15, scrolling=False)
