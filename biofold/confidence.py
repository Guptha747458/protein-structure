"""pLDDT confidence assessment, B-factor parsing, and confidence track plots."""

from __future__ import annotations

import statistics
from collections import Counter
from dataclasses import asdict, dataclass
from typing import Any

import plotly.graph_objects as go

from biofold.constants import PLDDT_TIERS


@dataclass(frozen=True)
class ResidueConfidence:
    position: int
    res_name: str
    plddt: float
    tier: str
    color: str


@dataclass(frozen=True)
class ConfidenceSummary:
    total_atoms: int
    total_residues: int
    mean_plddt: float
    median_plddt: float
    min_plddt: float
    max_plddt: float
    std_plddt: float
    class_counts: dict[str, int]
    class_percentages: dict[str, float]
    residues: list[ResidueConfidence]

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["residues"] = [asdict(r) for r in self.residues]
        return data


def classify_plddt(score: float) -> tuple[str, str]:
    """Map a pLDDT score to tier category and standard color."""
    if score > 90.0:
        return "Very High", PLDDT_TIERS["Very High"]["color"]
    if score > 70.0:
        return "Confident", PLDDT_TIERS["Confident"]["color"]
    if score > 50.0:
        return "Low", PLDDT_TIERS["Low"]["color"]
    return "Very Low", PLDDT_TIERS["Very Low"]["color"]


def extract_plddt_from_pdb(pdb_text: str) -> ConfidenceSummary:
    """Extract B-factor values from PDB records as pLDDT scores.

    Groups atoms by residue sequence ID (preferring CA atom B-factor when available,
    otherwise computing residue average atom B-factor).
    """
    if not pdb_text:
        return ConfidenceSummary(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, {}, [])

    atom_b_factors: list[float] = []
    residue_data: dict[str, dict[str, Any]] = {}

    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                res_seq_str = line[22:26].strip()
                b_factor = float(line[60:66].strip())
            except (ValueError, IndexError):
                continue

            # Ignore crystallographic waters and common buffer/ion additives
            if res_name in {"HOH", "WAT", "H2O", "DOD", "SO4", "PO4", "GOL", "EDO", "NA", "CL", "ZN", "MG", "CA"}:
                continue

            atom_b_factors.append(b_factor)
            key = f"{chain_id}:{res_seq_str}"
            if key not in residue_data:
                residue_data[key] = {
                    "res_name": res_name,
                    "ca_b": None,
                    "all_b": [],
                    "order": len(residue_data) + 1,
                }
            residue_data[key]["all_b"].append(b_factor)
            if atom_name == "CA":
                residue_data[key]["ca_b"] = b_factor

    if not atom_b_factors:
        return ConfidenceSummary(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, {}, [])

    residues: list[ResidueConfidence] = []
    res_scores: list[float] = []

    for idx, (key, data) in enumerate(residue_data.items(), start=1):
        score = data["ca_b"] if data["ca_b"] is not None else statistics.mean(data["all_b"])
        score = round(max(0.0, min(100.0, score)), 2)
        res_scores.append(score)
        tier, color = classify_plddt(score)
        residues.append(
            ResidueConfidence(
                position=idx,
                res_name=data["res_name"],
                plddt=score,
                tier=tier,
                color=color,
            )
        )

    eval_scores = res_scores if res_scores else atom_b_factors
    tier_counts = Counter(classify_plddt(s)[0] for s in eval_scores)
    total_items = len(eval_scores)

    all_tiers = ["Very High", "Confident", "Low", "Very Low"]
    counts_map = {t: tier_counts.get(t, 0) for t in all_tiers}
    pcts_map = {t: round((counts_map[t] / total_items) * 100.0, 1) for t in all_tiers}

    return ConfidenceSummary(
        total_atoms=len(atom_b_factors),
        total_residues=len(residues),
        mean_plddt=round(statistics.mean(eval_scores), 2),
        median_plddt=round(statistics.median(eval_scores), 2),
        min_plddt=round(min(eval_scores), 2),
        max_plddt=round(max(eval_scores), 2),
        std_plddt=round(statistics.stdev(eval_scores), 2) if len(eval_scores) > 1 else 0.0,
        class_counts=counts_map,
        class_percentages=pcts_map,
        residues=residues,
    )


def create_plddt_track_chart(summary: ConfidenceSummary) -> go.Figure:
    """Generate an interactive per-residue pLDDT confidence line & scatter plot."""
    if not summary.residues:
        fig = go.Figure()
        fig.update_layout(title="No pLDDT confidence data available.")
        return fig

    x_vals = [r.position for r in summary.residues]
    y_vals = [r.plddt for r in summary.residues]
    hover_texts = [
        f"<b>Residue {r.position} ({r.res_name})</b><br>pLDDT: <b>{r.plddt:.2f}</b><br>Confidence: {r.tier}"
        for r in summary.residues
    ]
    colors = [r.color for r in summary.residues]

    fig = go.Figure()

    # Background color bands for AlphaFold tiers
    fig.add_hrect(y0=90, y1=100, fillcolor="#0053D6", opacity=0.08, line_width=0, annotation_text="Very High (>90)", annotation_position="top right")
    fig.add_hrect(y0=70, y1=90, fillcolor="#65CBF3", opacity=0.08, line_width=0, annotation_text="Confident (70-90)", annotation_position="top right")
    fig.add_hrect(y0=50, y1=70, fillcolor="#FFDB13", opacity=0.08, line_width=0, annotation_text="Low (50-70)", annotation_position="top right")
    fig.add_hrect(y0=0, y1=50, fillcolor="#FF7D45", opacity=0.08, line_width=0, annotation_text="Very Low (≤50)", annotation_position="top right")

    # Trend line connecting residues
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="lines",
            line=dict(color="#4B5563", width=1.5, dash="dot"),
            showlegend=False,
            hoverinfo="skip",
        )
    )

    # Residue markers colored by tier
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="markers",
            marker=dict(size=7, color=colors, line=dict(width=1, color="#1F2937")),
            text=hover_texts,
            hoverinfo="text",
            name="pLDDT Score",
        )
    )

    # Average line
    fig.add_hline(
        y=summary.mean_plddt,
        line_dash="dash",
        line_color="#E11D48",
        line_width=2,
        annotation_text=f"Mean pLDDT: {summary.mean_plddt:.1f}",
        annotation_position="bottom left",
    )

    fig.update_layout(
        title=dict(text="<b>Per-Residue pLDDT Confidence Profile</b>", font=dict(size=14, color="#F8FAFC")),
        xaxis=dict(
            title="Residue Sequence Index",
            dtick=max(1, len(x_vals) // 10),
            showgrid=True,
            gridcolor="rgba(148, 163, 184, 0.15)",
            tickfont=dict(color="#CBD5E1"),
            title_font=dict(color="#F1F5F9"),
        ),
        yaxis=dict(
            title="pLDDT Confidence (0-100)",
            range=[0, 105],
            showgrid=True,
            gridcolor="rgba(148, 163, 184, 0.15)",
            tickfont=dict(color="#CBD5E1"),
            title_font=dict(color="#F1F5F9"),
        ),
        height=300,
        margin=dict(l=40, r=20, t=40, b=40),
        plot_bgcolor="rgba(15, 23, 42, 0.4)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


def create_confidence_donut_chart(summary: ConfidenceSummary) -> go.Figure:
    """Generate a donut chart showing proportions of residues in each pLDDT tier."""
    labels = ["Very High (>90)", "Confident (70-90)", "Low (50-70)", "Very Low (≤50)"]
    keys = ["Very High", "Confident", "Low", "Very Low"]
    values = [summary.class_counts.get(k, 0) for k in keys]
    colors = [PLDDT_TIERS[k]["color"] for k in keys]

    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.6,
                marker=dict(colors=colors, line=dict(color="#FFFFFF", width=2)),
                textinfo="percent+label",
                textposition="outside",
                hoverinfo="label+value+percent",
            )
        ]
    )
    fig.update_layout(
        showlegend=False,
        height=260,
        margin=dict(l=20, r=20, t=20, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
        annotations=[
            dict(
                text=f"<b>{summary.mean_plddt:.1f}</b><br><span style='font-size:10px'>Mean pLDDT</span>",
                x=0.5,
                y=0.5,
                font_size=18,
                showarrow=False,
            )
        ],
    )
    return fig
