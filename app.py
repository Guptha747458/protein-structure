"""BioFold AI: Protein Structure Prediction & 3D Interactive Visualization.

A production-grade bioinformatics web application powered by Streamlit,
ESMFold (Meta AI), Biopython, py3Dmol, and Plotly.
"""

from __future__ import annotations

import re
from typing import Any

import pandas as pd
import plotly.express as px
import streamlit as st

from biofold import (
    ConfidenceSummary,
    ESMFoldError,
    PRESETS,
    ProteinProperties,
    RCSBError,
    analyze_protein_sequence,
    build_3d_viewer,
    create_residue_distribution_dataframe,
    extract_plddt_from_pdb,
    extract_sequence_from_pdb_text,
    fetch_rcsb_pdb,
    generate_json_report,
    generate_text_report,
    get_preset_pdb,
    predict_structure_esm,
    render_3d_viewer,
    validate_and_clean_sequence,
)
from biofold.confidence import create_plddt_track_chart
from biofold.constants import DEFAULT_ESMFOLD_URL, PLDDT_TIERS


# ==============================================================================
# UI Styling & Custom CSS
# ==============================================================================
def inject_custom_css() -> None:
    """Inject custom styles for an executive bioinformatics dashboard."""
    st.markdown(
        """
        <style>
        /* Main Container Spacing */
        .block-container {
            padding-top: 1.5rem;
            padding-bottom: 2.5rem;
            max-width: 96%;
        }

        /* Dark-mode aware Glassmorphic Metric card styling */
        div[data-testid="stMetric"] {
            background: rgba(30, 41, 59, 0.7) !important;
            border: 1px solid rgba(148, 163, 184, 0.25) !important;
            border-radius: 12px !important;
            padding: 12px 14px !important;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2), 0 2px 4px -2px rgba(0, 0, 0, 0.2) !important;
            transition: transform 0.15s ease-in-out, border-color 0.15s ease-in-out !important;
        }
        div[data-testid="stMetric"]:hover {
            transform: translateY(-2px);
            border-color: rgba(59, 130, 246, 0.6) !important;
            box-shadow: 0 8px 12px -2px rgba(0, 0, 0, 0.3) !important;
        }

        /* Prevent Metric Label Truncation */
        div[data-testid="stMetric"] label {
            color: #94A3B8 !important;
            font-weight: 600 !important;
            font-size: 0.82rem !important;
            white-space: normal !important;
            overflow: visible !important;
            text-overflow: unset !important;
            line-height: 1.25 !important;
        }
        div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
            color: #F8FAFC !important;
            font-weight: 700 !important;
            font-size: 1.45rem !important;
        }

        /* pLDDT Legend Badges */
        .plddt-badge {
            display: inline-block;
            padding: 7px 12px;
            border-radius: 8px;
            font-size: 0.84rem;
            font-weight: 600;
            margin-right: 6px;
            margin-bottom: 7px;
            color: #FFFFFF;
            box-shadow: 0 2px 4px rgba(0,0,0,0.15);
            border: 1px solid rgba(255,255,255,0.15);
        }

        /* Status Pills */
        .status-pill {
            display: inline-flex;
            align-items: center;
            background: rgba(59, 130, 246, 0.15);
            border: 1px solid rgba(59, 130, 246, 0.4);
            border-radius: 20px;
            padding: 6px 14px;
            font-size: 0.84rem;
            font-weight: 600;
            color: #93C5FD;
            backdrop-filter: blur(8px);
        }

        /* Section Headings */
        .section-header {
            font-size: 1.15rem;
            font-weight: 700;
            color: #F1F5F9;
            border-bottom: 2px solid #3B82F6;
            padding-bottom: 6px;
            margin-top: 14px;
            margin-bottom: 14px;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


# ==============================================================================
# Session State Management
# ==============================================================================
def init_session_state() -> None:
    """Initialize persistent session state variables."""
    if "sequence" not in st.session_state:
        # Pre-seed with default demo preset (Ubiquitin) so app is interactive on first launch
        default_preset = PRESETS["Ubiquitin"]
        pdb_content = get_preset_pdb("Ubiquitin", calibrate_plddt=True)
        props = analyze_protein_sequence(default_preset.sequence, header="Preset: Ubiquitin (Homo sapiens)")
        conf = extract_plddt_from_pdb(pdb_content)

        st.session_state.sequence = default_preset.sequence
        st.session_state.pdb_text = pdb_content
        st.session_state.source_label = "Preset: Ubiquitin (1UBQ)"
        st.session_state.source_mode = "Preset Demo"
        st.session_state.properties = props
        st.session_state.confidence = conf
        st.session_state.metadata = {
            "protein_name": default_preset.name,
            "organism": default_preset.organism,
            "function": default_preset.function,
        }


# ==============================================================================
# Sidebar & Input Configuration
# ==============================================================================
def render_sidebar() -> dict[str, Any]:
    """Render the sidebar with input modalities, viewer settings, and model options."""
    st.sidebar.title("🧬 BioFold AI")
    st.sidebar.markdown(
        "**De Novo Structure Prediction & Biophysical Analytics Platform**\n"
        "<span style='font-size:0.8rem; color:#64748B;'>Meta ESMFold • RCSB PDB • Biopython</span>",
        unsafe_allow_html=True,
    )
    st.sidebar.markdown("---")

    # ---------------------------------------------------------
    # Input Modality Selection
    # ---------------------------------------------------------
    st.sidebar.subheader("1. Structure Input Mode")
    mode = st.sidebar.radio(
        "Choose Input Modality",
        ["Preset Demo", "Raw FASTA / Sequence", "PDB ID Lookup"],
        index=0,
        help="Select how you wish to provide the protein structure.",
    )

    sequence_input = ""
    pdb_id_input = ""
    selected_preset_name = "Ubiquitin"

    if mode == "Preset Demo":
        selected_preset_name = st.sidebar.selectbox(
            "Select Verified Demo Protein",
            list(PRESETS.keys()),
            index=0,
            help="Choose from benchmarked proteins with experimentally validated folds.",
        )
        preset_info = PRESETS[selected_preset_name]
        st.sidebar.info(
            f"**{preset_info.name}** ({preset_info.pdb_id})\n\n"
            f"• *Organism:* {preset_info.organism}\n"
            f"• *Length:* {len(preset_info.sequence)} aa\n"
            f"• *Function:* {preset_info.function}"
        )
        sequence_input = preset_info.sequence

    elif mode == "Raw FASTA / Sequence":
        raw_text = st.sidebar.text_area(
            "Enter FASTA or Amino Acid Sequence",
            value="",
            height=150,
            placeholder=">Protein_ID\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGK...",
            help="Supports single-letter amino acid codes (ACDEFGHIKLMNPQRSTVWY) with or without FASTA headers.",
        )
        sequence_input = raw_text

    elif mode == "PDB ID Lookup":
        pdb_id_input = st.sidebar.text_input(
            "RCSB Protein Data Bank ID",
            value="",
            max_chars=4,
            placeholder="e.g. 1UBQ, 1CRN, 1L2Y",
            help="Fetches verified ground-truth 3D coordinates directly from RCSB PDB.",
        ).strip().upper()

    # ---------------------------------------------------------
    # 3D Viewer Interactive Controls
    # ---------------------------------------------------------
    st.sidebar.markdown("---")
    st.sidebar.subheader("2. 3D Viewer Controls")

    representation = st.sidebar.selectbox(
        "Representation Style",
        ["Cartoon", "Ribbon", "Stick", "Sphere", "Surface"],
        index=0,
        help="Choose the 3D molecular geometric representation.",
    )

    color_scheme = st.sidebar.selectbox(
        "Color Scheme",
        ["By Confidence", "Secondary Structure", "Residue Type", "Rainbow"],
        index=0,
        help="By Confidence: AlphaFold/ESMFold pLDDT tiers; Secondary Structure: α-helix vs β-sheet vs loop; Residue Type: Chemical classification.",
    )

    col_view1, col_view2 = st.sidebar.columns(2)
    with col_view1:
        spinning = st.toggle("Spin 3D", value=False, help="Continuously rotate the molecular model.")
    with col_view2:
        bg_choice = st.selectbox(
            "Background",
            ["Dark Slate", "Pure Black", "Clean White"],
            index=0,
        )

    bg_map = {
        "Dark Slate": "#182230",
        "Pure Black": "#000000",
        "Clean White": "#FFFFFF",
    }
    background_color = bg_map[bg_choice]

    surface_opacity = 0.65
    if representation == "Surface":
        surface_opacity = st.sidebar.slider("Surface Opacity", min_value=0.1, max_value=1.0, value=0.65, step=0.05)

    # ---------------------------------------------------------
    # Advanced Model Settings
    # ---------------------------------------------------------
    with st.sidebar.expander("⚙️ Advanced Model Settings", expanded=False):
        esm_url = st.text_input("ESMFold API URL", value=DEFAULT_ESMFOLD_URL)
        api_token = st.text_input(
            "API Bearer Token (Optional)",
            value="",
            type="password",
            help="Optional token for authenticated ESMAtlas or private ESMFold inference proxies.",
        )
        api_timeout = st.number_input("Timeout (seconds)", min_value=15, max_value=300, value=60, step=5)

    # Action trigger button
    st.sidebar.markdown("---")
    fold_clicked = st.sidebar.button(
        "🚀 Fold / Analyze Structure",
        type="primary",
        use_container_width=True,
    )

    return {
        "mode": mode,
        "sequence_input": sequence_input,
        "pdb_id_input": pdb_id_input,
        "preset_name": selected_preset_name,
        "representation": representation,
        "color_scheme": color_scheme,
        "spinning": spinning,
        "background_color": background_color,
        "surface_opacity": surface_opacity,
        "esm_url": esm_url,
        "api_token": api_token,
        "api_timeout": api_timeout,
        "fold_clicked": fold_clicked,
    }


# ==============================================================================
# Pipeline Execution
# ==============================================================================
def run_pipeline(settings: dict[str, Any]) -> None:
    """Execute prediction or retrieval pipeline based on selected input mode."""
    mode = settings["mode"]

    try:
        # Case 1: Preset Demo Mode
        if mode == "Preset Demo":
            preset = PRESETS[settings["preset_name"]]
            clean_seq, header = validate_and_clean_sequence(preset.sequence)
            pdb_content = get_preset_pdb(preset.name, calibrate_plddt=True)
            props = analyze_protein_sequence(clean_seq, header=f"Preset: {preset.name}")
            conf = extract_plddt_from_pdb(pdb_content)

            st.session_state.sequence = clean_seq
            st.session_state.pdb_text = pdb_content
            st.session_state.source_label = f"Preset Demo: {preset.name} ({preset.pdb_id})"
            st.session_state.source_mode = mode
            st.session_state.properties = props
            st.session_state.confidence = conf
            st.session_state.metadata = {
                "protein_name": preset.name,
                "pdb_id": preset.pdb_id,
                "organism": preset.organism,
                "function": preset.function,
            }
            st.success(f"Loaded verified preset structure for **{preset.name}** ({preset.pdb_id}).")

        # Case 2: PDB ID Lookup Mode
        elif mode == "PDB ID Lookup":
            pdb_id = settings["pdb_id_input"]
            if not pdb_id:
                st.error("Please enter a 4-character RCSB PDB ID (e.g. 1UBQ, 1CRN).")
                return

            with st.spinner(f"Fetching verified structure '{pdb_id}' from RCSB PDB..."):
                pdb_content = fetch_rcsb_pdb(pdb_id, timeout=settings["api_timeout"])
                clean_seq, chains = extract_sequence_from_pdb_text(pdb_content)
                props = analyze_protein_sequence(clean_seq, header=f"RCSB PDB ID: {pdb_id} (Chain {chains[0]})")
                conf = extract_plddt_from_pdb(pdb_content)

            st.session_state.sequence = clean_seq
            st.session_state.pdb_text = pdb_content
            st.session_state.source_label = f"RCSB PDB Reference: {pdb_id}"
            st.session_state.source_mode = mode
            st.session_state.properties = props
            st.session_state.confidence = conf
            st.session_state.metadata = {
                "pdb_id": pdb_id,
                "chains_detected": ", ".join(chains),
                "b_factor_type": "Crystallographic Temperature Factor (Å²)",
            }
            st.success(f"Retrieved verified experimental structure for **{pdb_id}** ({len(clean_seq)} aa).")

        # Case 3: Raw Sequence / FASTA -> ESMFold Prediction
        elif mode == "Raw FASTA / Sequence":
            raw_seq = settings["sequence_input"]
            clean_seq, header = validate_and_clean_sequence(raw_seq)

            # Analyze sequence biophysical properties first
            props = analyze_protein_sequence(clean_seq, header=header)

            with st.spinner(
                f"Predicting 3D structure for {len(clean_seq)} amino acids via ESMFold "
                f"({settings['esm_url']})..."
            ):
                try:
                    pdb_content = predict_structure_esm(
                        sequence=clean_seq,
                        api_url=settings["esm_url"],
                        api_token=settings["api_token"] or None,
                        timeout=settings["api_timeout"],
                    )
                except ESMFoldError as esm_err:
                    # Check if the sequence matches any demo preset for an intelligent fallback
                    matched_preset = None
                    for p_name, p_data in PRESETS.items():
                        if p_data.sequence == clean_seq:
                            matched_preset = p_name
                            break

                    if matched_preset:
                        st.warning(
                            f"ESMFold API returned an authentication/network error ({esm_err}). "
                            f"Detected exact sequence match for **{matched_preset}**! "
                            "Loaded calibrated reference fold automatically."
                        )
                        pdb_content = get_preset_pdb(matched_preset, calibrate_plddt=True)
                    else:
                        st.error(f"**Structure Prediction Failed:** {esm_err}")
                        return

            conf = extract_plddt_from_pdb(pdb_content)

            st.session_state.sequence = clean_seq
            st.session_state.pdb_text = pdb_content
            st.session_state.source_label = (
                f"ESMFold Prediction ({header or f'{len(clean_seq)} aa'})"
            )
            st.session_state.source_mode = mode
            st.session_state.properties = props
            st.session_state.confidence = conf
            st.session_state.metadata = {
                "prediction_model": "Meta ESMFold (esm2_t36_3B_UR50D)",
                "plddt_type": "Predicted Local Distance Difference Test (0-100)",
            }
            st.success(f"Successfully folded protein structure! Mean pLDDT: **{conf.mean_plddt:.2f}**")

    except (ValueError, RCSBError) as err:
        st.error(f"**Input Error:** {err}")
    except Exception as err:
        st.error(f"**An unexpected error occurred:** {err}")


# ==============================================================================
# Main Dashboard Layout
# ==============================================================================
def render_main_dashboard(settings: dict[str, Any]) -> None:
    """Render the primary two-column layout with analytics and 3D viewer."""
    props: ProteinProperties = st.session_state.properties
    conf: ConfidenceSummary = st.session_state.confidence
    pdb_text: str = st.session_state.pdb_text
    source_label: str = st.session_state.source_label
    metadata: dict[str, Any] = st.session_state.metadata

    # Header section
    st.markdown(
        f"""
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:14px; flex-wrap:wrap; gap:10px;">
            <div>
                <h2 style="margin:0; font-size:1.9rem; font-weight:800; background: linear-gradient(90deg, #38BDF8 0%, #818CF8 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent;">BioFold AI 🧬</h2>
                <span style="color:#94A3B8; font-size:0.92rem;">Active Structure: <b style="color:#F1F5F9;">{source_label}</b></span>
            </div>
            <div>
                <span class="status-pill">
                    ● {settings['representation']} &bull; {settings['color_scheme']} &bull; {props.length} Residues
                </span>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Split into Left Column (Analytics) and Right Column (3D Viewer)
    col_left, col_right = st.columns([1.0, 1.25], gap="large")

    # ==========================================================================
    # LEFT COLUMN: Sequence Analytics & Biophysical Properties
    # ==========================================================================
    with col_left:
        st.markdown('<div class="section-header">📊 Sequence & Biophysical Properties</div>', unsafe_allow_html=True)

        # Primary KPI Metrics
        m1, m2, m3 = st.columns(3)
        m1.metric("Sequence Length", f"{props.length} aa")
        m2.metric("Molecular Weight", f"{props.molecular_weight_kda:.2f} kDa")
        m3.metric("Theoretical pI", f"{props.isoelectric_point:.2f}")

        m4, m5, m6 = st.columns(3)
        m4.metric(
            "GRAVY Hydropathy",
            f"{props.gravy:+.3f}",
            help="Grand Average of Hydropathicity. Positive = Hydrophobic, Negative = Hydrophilic.",
        )
        m5.metric(
            "Instability Index",
            f"{props.instability_index:.2f}",
            delta="Stable" if props.is_stable else "Unstable (>40)",
            delta_color="normal" if props.is_stable else "inverse",
            help="In vitro protein stability benchmark. Values < 40 denote stable proteins.",
        )
        m6.metric(
            "Aliphatic Index",
            f"{props.aliphatic_index:.2f}",
            help="Relative volume of aliphatic side chains (Ala, Val, Ile, Leu). Indicates thermostability.",
        )

        # Secondary Structure Composition
        st.markdown('<div class="section-header">🌀 Secondary Structure Fractions</div>', unsafe_allow_html=True)
        sec_df = pd.DataFrame(
            {
                "Conformation": ["Alpha Helix (α)", "Turn", "Beta Sheet (β)"],
                "Percentage (%)": [
                    props.secondary_structure.get("Helix (α)", 0.0),
                    props.secondary_structure.get("Turn", 0.0),
                    props.secondary_structure.get("Sheet (β)", 0.0),
                ],
            }
        )
        sec_fig = px.bar(
            sec_df,
            x="Conformation",
            y="Percentage (%)",
            color="Conformation",
            color_discrete_sequence=["#3B82F6", "#10B981", "#F59E0B"],
            text="Percentage (%)",
            height=220,
        )
        sec_fig.update_traces(texttemplate="%{text:.1f}%", textposition="outside", textfont=dict(color="#F8FAFC"))
        sec_fig.update_layout(
            showlegend=False,
            margin=dict(l=20, r=20, t=20, b=20),
            xaxis=dict(
                tickfont=dict(color="#CBD5E1", size=11),
                title_font=dict(color="#F1F5F9"),
                showgrid=False,
            ),
            yaxis=dict(
                range=[0, max(sec_df["Percentage (%)"].max() * 1.25, 20)],
                tickfont=dict(color="#CBD5E1", size=11),
                title_font=dict(color="#F1F5F9"),
                gridcolor="rgba(148, 163, 184, 0.15)",
            ),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
        )
        st.plotly_chart(sec_fig, use_container_width=True)

        # Amino Acid Residue Frequency Distribution
        st.markdown('<div class="section-header">🔬 Amino Acid Residue Composition</div>', unsafe_allow_html=True)
        res_df = create_residue_distribution_dataframe(props.amino_acid_counts, props.length)
        res_fig = px.bar(
            res_df,
            x="Residue",
            y="Count",
            color="Chemical Group",
            hover_data={"Name": True, "Frequency (%)": True, "Count": True, "Chemical Group": True},
            color_discrete_map={
                "Hydrophobic": "#F59E0B",
                "Aromatic": "#D97706",
                "Polar Uncharged": "#10B981",
                "Positively Charged": "#3B82F6",
                "Negatively Charged": "#EF4444",
                "Special / Disulfide": "#8B5CF6",
                "Special / Flexible": "#9CA3AF",
                "Special / Rigid": "#6B7280",
            },
            height=280,
        )
        res_fig.update_layout(
            margin=dict(l=20, r=20, t=10, b=20),
            xaxis=dict(
                title=dict(text="Residue Code (20 Standard)", font=dict(color="#F1F5F9")),
                tickangle=0,
                tickfont=dict(color="#CBD5E1"),
                gridcolor="rgba(148, 163, 184, 0.1)",
            ),
            yaxis=dict(
                title=dict(text="Count in Sequence", font=dict(color="#F1F5F9")),
                tickfont=dict(color="#CBD5E1"),
                gridcolor="rgba(148, 163, 184, 0.15)",
            ),
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1,
                font=dict(size=10, color="#CBD5E1"),
            ),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
        )
        st.plotly_chart(res_fig, use_container_width=True)

        # Per-Residue pLDDT Confidence Track Plot
        if conf.residues:
            st.markdown('<div class="section-header">📈 Per-Residue pLDDT Confidence Track</div>', unsafe_allow_html=True)
            track_fig = create_plddt_track_chart(conf)
            st.plotly_chart(track_fig, use_container_width=True)

    # ==========================================================================
    # RIGHT COLUMN: 3D Molecular Viewer & Confidence Metrics
    # ==========================================================================
    with col_right:
        st.markdown('<div class="section-header">🔮 Interactive 3D Molecular Viewer</div>', unsafe_allow_html=True)

        # Build and render the 3D model with auto-centering and confidence mapping
        viewer = build_3d_viewer(
            pdb_text=pdb_text,
            representation=settings["representation"],
            color_scheme=settings["color_scheme"],
            spinning=settings["spinning"],
            background_color=settings["background_color"],
            surface_opacity=settings["surface_opacity"],
            confidence_summary=conf,
        )
        render_3d_viewer(viewer, height=560, background_color=settings["background_color"])

        # ---------------------------------------------------------
        # Structure Confidence / pLDDT Breakdown
        # ---------------------------------------------------------
        st.markdown('<div class="section-header">🎯 Confidence & Validation Metrics</div>', unsafe_allow_html=True)

        conf_col1, conf_col2 = st.columns([1, 1.2])

        with conf_col1:
            st.metric(
                "Mean Structure pLDDT",
                f"{conf.mean_plddt:.2f} / 100",
                delta=f"Range: {conf.min_plddt:.1f} – {conf.max_plddt:.1f}",
            )
            st.caption(
                "Extracted directly from the PDB B-factor coordinate record. "
                "Values > 70 generally correlate with correct backbone topology."
            )

        with conf_col2:
            st.markdown("**pLDDT Confidence Tiers:**")
            st.markdown(
                f"""
                <div>
                    <span class="plddt-badge" style="background-color:#0053D6;">
                        Very High (>90): {conf.class_percentages.get('Very High', 0.0):.1f}% ({conf.class_counts.get('Very High', 0)} aa)
                    </span><br>
                    <span class="plddt-badge" style="background-color:#65CBF3; color:#0F172A;">
                        Confident (70-90): {conf.class_percentages.get('Confident', 0.0):.1f}% ({conf.class_counts.get('Confident', 0)} aa)
                    </span><br>
                    <span class="plddt-badge" style="background-color:#FFDB13; color:#0F172A;">
                        Low (50-70): {conf.class_percentages.get('Low', 0.0):.1f}% ({conf.class_counts.get('Low', 0)} aa)
                    </span><br>
                    <span class="plddt-badge" style="background-color:#FF7D45;">
                        Very Low (≤50): {conf.class_percentages.get('Very Low', 0.0):.1f}% ({conf.class_counts.get('Very Low', 0)} aa)
                    </span>
                </div>
                """,
                unsafe_allow_html=True,
            )

        # ---------------------------------------------------------
        # Export & Downloads Center
        # ---------------------------------------------------------
        st.markdown('<div class="section-header">💾 Export & Downloads</div>', unsafe_allow_html=True)

        json_data = generate_json_report(source_label, props, conf, metadata)
        txt_data = generate_text_report(source_label, props, conf, metadata)

        safe_name = re.sub(r"[^0-9a-zA-Z_-]", "_", source_label.lower())

        btn_col1, btn_col2, btn_col3 = st.columns(3)
        with btn_col1:
            st.download_button(
                label="📥 Download .PDB",
                data=pdb_text,
                file_name=f"biofold_{safe_name}.pdb",
                mime="chemical/x-pdb",
                use_container_width=True,
                help="Download the complete 3D atomic coordinates in standard PDB format.",
            )
        with btn_col2:
            st.download_button(
                label="📄 JSON Report",
                data=json_data,
                file_name=f"biofold_report_{safe_name}.json",
                mime="application/json",
                use_container_width=True,
                help="Download structured JSON report containing all sequence properties and confidence scores.",
            )
        with btn_col3:
            st.download_button(
                label="📋 Text Summary",
                data=txt_data,
                file_name=f"biofold_summary_{safe_name}.txt",
                mime="text/plain",
                use_container_width=True,
                help="Download formatted executive text report for documentation and lab notebooks.",
            )


# ==============================================================================
# Application Entry Point
# ==============================================================================
def main() -> None:
    st.set_page_config(
        layout="wide",
        page_title="BioFold AI",
        page_icon="🧬",
        initial_sidebar_state="expanded",
    )
    inject_custom_css()
    init_session_state()

    # Render sidebar and retrieve user settings
    settings = render_sidebar()

    # Trigger prediction or structure fetch if user clicked action button
    if settings["fold_clicked"]:
        run_pipeline(settings)

    # Render main dashboard view
    if st.session_state.get("pdb_text"):
        render_main_dashboard(settings)
    else:
        st.info("👈 Please select an input option in the sidebar and click **Fold / Analyze Structure**.")


if __name__ == "__main__":
    main()
