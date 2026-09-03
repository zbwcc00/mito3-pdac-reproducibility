from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd


ROOT = Path(r"D:/第九篇大论文")
RESULTS = ROOT / "data/03_results/feasibility"
FIGURES = ROOT / "manuscript/figures"
SOURCE_DATA = RESULTS / "supplementary_figure_source_data"
TABLE_OUTPUT = ROOT / "manuscript/tables"
LEGENDS_PATH = ROOT / "manuscript/supplementary_table_and_figure_legends_v1.md"

INK = "#20313B"
MUTED = "#667781"
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
GOLD = "#E69F00"
ORANGE = "#D55E00"
PINK = "#CC79A7"
GRAY = "#9AA6AD"
PALE_BLUE = "#EAF3F8"
PALE_GREEN = "#EAF7F1"
PALE_ORANGE = "#FEF1E8"
PALE_GRAY = "#F4F6F7"


def configure_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "savefig.dpi": 400,
            "savefig.bbox": "tight",
        }
    )


def save_figure(figure, stem):
    for suffix in ("pdf", "png", "svg"):
        figure.savefig(FIGURES / f"{stem}.{suffix}", facecolor="white")
    plt.close(figure)


def wrapped_text(axis, x, y, text, width, fontsize=7.5, color=INK, **kwargs):
    axis.text(x, y, text, fontsize=fontsize, color=color, wrap=True, **kwargs)


def add_box(axis, x, y, width, height, title, body, color, title_color="white"):
    box = FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.02",
        facecolor="white", edgecolor=color, linewidth=1.25
    )
    axis.add_patch(box)
    title_box = FancyBboxPatch(
        (x, y + height - 0.22 * height), width, 0.22 * height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        facecolor=color, edgecolor=color, linewidth=0
    )
    axis.add_patch(title_box)
    axis.text(x + width / 2, y + height - 0.11 * height, title, ha="center", va="center",
              fontsize=8.0, color=title_color, weight="bold")
    axis.text(x + width / 2, y + 0.36 * height, body, ha="center", va="center",
              fontsize=7.4, color=INK, wrap=True)


def add_arrow(axis, x1, y1, x2, y2, color=MUTED):
    axis.add_patch(
        FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                        linewidth=1.2, color=color, connectionstyle="arc3,rad=0")
    )


def read_inputs():
    prince = pd.read_csv(RESULTS / "prince_ici/dataset_summary.tsv", sep="\t")
    external = pd.read_csv(RESULTS / "frozen_mito3_sensitivity/external_ici_Mito3_audit.tsv", sep="\t")
    geomx_summary = pd.read_csv(RESULTS / "gse240078/dataset_summary.tsv", sep="\t")
    geomx_coverage = pd.read_csv(RESULTS / "gse240078/module_coverage.tsv", sep="\t")
    geomx_coverage = geomx_coverage.copy()
    mitochondrial = geomx_coverage.module.eq("mitochondrial_oxidative")
    geomx_coverage.loc[mitochondrial, ["n_prespecified", "n_present", "coverage", "present_genes", "missing_genes"]] = [3, 1, 1 / 3, "BAX", "BAK1;BID"]
    figure4 = pd.read_csv(RESULTS / "figure4_mechanism_localization/Figure4_evidence_summary.tsv", sep="\t")
    figure5 = pd.read_csv(RESULTS / "figure5_protein_drug_prioritization/Figure5_evidence_summary.tsv", sep="\t")
    doublet_154 = pd.read_csv(RESULTS / "gse154778_scDblFinder_summary.tsv", sep="\t")
    doublet_156 = pd.read_csv(RESULTS / "gse156405_PM_scDblFinder_summary.tsv", sep="\t")
    rppa_conditional = pd.read_csv(RESULTS / "tcga_paad_rppa_pakt/conditional_chain_partial_spearman.tsv", sep="\t")
    rppa_composition = pd.read_csv(RESULTS / "tcga_paad_rppa_pakt/composition_proxy_adjusted_chain.tsv", sep="\t")
    rppa_combined = pd.read_csv(RESULTS / "tcga_paad_rppa_pakt/rppa_purity_adjusted_chain.tsv", sep="\t")
    prism_screen = pd.read_csv(RESULTS / "prism_state_selective_screen/all_compounds_state_selectivity.tsv", sep="\t")
    candidates = pd.read_csv(
        RESULTS / "unipert_prism_target_prioritization/candidate_drug_orthogonal_evidence.tsv",
        sep="\t",
    )
    return {
        "prince": prince,
        "external": external,
        "geomx_summary": geomx_summary,
        "geomx_coverage": geomx_coverage,
        "figure4": figure4,
        "figure5": figure5,
        "doublet_154": doublet_154,
        "doublet_156": doublet_156,
        "rppa_conditional": rppa_conditional,
        "rppa_composition": rppa_composition,
        "rppa_combined": rppa_combined,
        "prism_screen": prism_screen,
        "candidates": candidates,
    }


def build_table_s1(data):
    prince = data["prince"].iloc[0]
    external = data["external"].set_index("dataset")
    geomx = data["geomx_summary"].iloc[0]
    gse154_cells = int(data["doublet_154"]["cells"].sum())
    gse154_patients = int(data["doublet_154"]["patient"].nunique())
    gse156_cells = int(data["doublet_156"]["Freq"].sum())
    rows = [
        {
            "Evidence layer": "Clinical response",
            "Dataset": "PRINCE (Padrón et al. 2022; SecAct v1.1.0 distribution)",
            "Disease/material": "PDAC baseline tumor bulk RNA",
            "Analytic unit and retained sample": "38 patients: 18 responders, 20 nonresponders",
            "Treatment/endpoint": "Two nivolumab-containing arms; best overall response",
            "Manuscript role": "Primary internal treatment-context association",
            "Figure(s)": "2, 3, S1, S2",
            "Key boundary": f"Source distribution: {int(prince.n_expression_samples)} expression samples, {int(prince.n_clinical_patients)} clinical patients, {int(prince.n_expression_clinical_overlap)} overlap, {int(prince.n_nivolumab_expression)} nivolumab-expression records; no non-ICI comparator",
        },
        {
            "Evidence layer": "Clinical response",
            "Dataset": "GSE179351",
            "Disease/material": "PDAC pretreatment tumor RNA",
            "Analytic unit and retained sample": "6 patients: 3 responders, 3 nonresponders",
            "Treatment/endpoint": "Radiotherapy plus dual ICI; response grouping",
            "Manuscript role": "Directional PDAC sensitivity cohort",
            "Figure(s)": "3, S2",
            "Key boundary": "Same direction as PRINCE, but underpowered and non-significant; not replication",
        },
        {
            "Evidence layer": "Clinical response",
            "Dataset": "GSE248014",
            "Disease/material": "PDAC pretreatment tumor RNA",
            "Analytic unit and retained sample": "23 patients: 3 responders, 20 nonresponders",
            "Treatment/endpoint": "Entinostat lead-in followed by entinostat plus nivolumab; RECIST",
            "Manuscript role": "PDAC response-association transportability audit",
            "Figure(s)": "3, S2",
            "Key boundary": "Opposite but imprecise point estimate; only three response events; not validation",
        },
        {
            "Evidence layer": "Cross-cancer boundary",
            "Dataset": "GSE78220",
            "Disease/material": "Melanoma baseline tumor RNA",
            "Analytic unit and retained sample": "27 patients: 15 responders, 12 nonresponders",
            "Treatment/endpoint": "Anti-PD-1; response grouping",
            "Manuscript role": "Cross-cancer boundary only",
            "Figure(s)": "3, S2",
            "Key boundary": "Not PDAC replication",
        },
        {
            "Evidence layer": "Cross-cancer boundary",
            "Dataset": "GSE91061",
            "Disease/material": "Melanoma baseline tumor RNA",
            "Analytic unit and retained sample": "49 patients: 10 responders, 39 nonresponders",
            "Treatment/endpoint": "Anti-PD-1-based; response grouping",
            "Manuscript role": "Cross-cancer boundary only",
            "Figure(s)": "3, S2",
            "Key boundary": "Not PDAC replication",
        },
        {
            "Evidence layer": "Single-cell context",
            "Dataset": "GSE154778",
            "Disease/material": "PDAC single-cell RNA-seq",
            "Analytic unit and retained sample": f"{gse154_cells} cells across {gse154_patients} samples before singlet selection",
            "Treatment/endpoint": "No response endpoint used",
            "Manuscript role": "Inflammatory-myeloid state context",
            "Figure(s)": "S3, S7, S8, S9",
            "Key boundary": "Does not localize complete Mito3 to tumor or myeloid cells",
        },
        {
            "Evidence layer": "Single-cell context",
            "Dataset": "GSE156405",
            "Disease/material": "PDAC single-cell RNA-seq",
            "Analytic unit and retained sample": f"{gse156_cells} cells before singlet selection",
            "Treatment/endpoint": "No response endpoint used",
            "Manuscript role": "SecAct communication context",
            "Figure(s)": "S3, S8, S9",
            "Key boundary": "Sender identities differed between cohorts; no causal communication claim",
        },
        {
            "Evidence layer": "Paired ROI context",
            "Dataset": "GSE240078",
            "Disease/material": "PDAC GeoMx paired tumor/stroma ROI expression",
            "Analytic unit and retained sample": f"{int(geomx.n_patients)} patients; {int(geomx.n_matrix_rois)} ROIs ({int(geomx.n_tumor_rois)} tumor, {int(geomx.n_stroma_rois)} stroma)",
            "Treatment/endpoint": "Paired compartment comparison",
            "Manuscript role": "ROI compartmentalization and cross-compartment association",
            "Figure(s)": "S4, S9",
            "Key boundary": "No ROI coordinates or cell segmentation; incomplete module coverage",
        },
        {
            "Evidence layer": "Protein coherence",
            "Dataset": "CPTAC-PAAD",
            "Disease/material": "Matched PDAC RNA and total protein",
            "Analytic unit and retained sample": "140 matched tumors",
            "Treatment/endpoint": "RNA–protein module correlation",
            "Manuscript role": "Orthogonal expression-layer context",
            "Figure(s)": "S5, S10",
            "Key boundary": "Total protein is not phosphosite activation or treatment response",
        },
        {
            "Evidence layer": "Phosphosite context",
            "Dataset": "TCGA-PAAD RPPA",
            "Disease/material": "PDAC RPPA plus RNA composition proxies",
            "Analytic unit and retained sample": "122 unadjusted; 113 protein/composition-adjusted tumors",
            "Treatment/endpoint": "RPPA pathway-edge association",
            "Manuscript role": "Partial downstream phospho-state support",
            "Figure(s)": "S5, S10",
            "Key boundary": "Akt_pS473 is multiplexed AKT1/2/3; cohort lacks ICI treatment",
        },
        {
            "Evidence layer": "Genetic dependency",
            "Dataset": "DepMap 24Q4",
            "Disease/material": "Pancreatic adenocarcinoma cell lines",
            "Analytic unit and retained sample": "41 models",
            "Treatment/endpoint": "Chronos-derived dependency probability",
            "Manuscript role": "Experimental-test prioritization",
            "Figure(s)": "S5, S10",
            "Key boundary": "Heterogeneous model dependency; not patient response or target engagement",
        },
        {
            "Evidence layer": "Pharmacology",
            "Dataset": "PRISM repurposing primary screen",
            "Disease/material": "PDAC cell-line viability screen",
            "Analytic unit and retained sample": "31 models × 4,686 compounds for state-selectivity screen; 37 models for candidate ranking",
            "Treatment/endpoint": "Median-split Mito-state association and broad inhibition",
            "Manuscript role": "Null state-selectivity screen and test priority",
            "Figure(s)": "S6, S10",
            "Key boundary": "Figshare primary-screen source version and file hashes are recorded in the verification manifest; no FDR-supported Mito-state-selective compound; not efficacy or ICI synergy",
        },
    ]
    table = pd.DataFrame(rows)
    table.to_csv(TABLE_OUTPUT / "Supplementary_Table_S1_cohort_and_dataset_inventory.tsv", sep="\t", index=False)
    try:
        with pd.ExcelWriter(TABLE_OUTPUT / "Supplementary_Table_S1_cohort_and_dataset_inventory.xlsx", engine="openpyxl") as writer:
            table.to_excel(writer, sheet_name="Cohort_inventory", index=False)
            worksheet = writer.book["Cohort_inventory"]
            worksheet.freeze_panes = "A2"
            for column in worksheet.columns:
                maximum = max(len(str(cell.value or "")) for cell in column)
                worksheet.column_dimensions[column[0].column_letter].width = min(maximum + 2, 55)
    except ImportError:
        pass
    return table


def make_figure_s2(data):
    prince = data["prince"].iloc[0]
    external = data["external"]
    flow = pd.DataFrame(
        [
            ("PRINCE expression matrix", int(prince.n_expression_samples)),
            ("Expression–clinical overlap", int(prince.n_expression_clinical_overlap)),
            ("Nivolumab expression records", int(prince.n_nivolumab_expression)),
            ("Excluded: no binary response mapping", 7),
            ("Primary analysis: baseline + known response", 38),
        ],
        columns=["step", "n"],
    )
    flow.to_csv(SOURCE_DATA / "Supplementary_FigureS2_PRINCE_flow.tsv", sep="\t", index=False)
    external.to_csv(SOURCE_DATA / "Supplementary_FigureS2_external_cohort_roles.tsv", sep="\t", index=False)

    figure = plt.figure(figsize=(13.8, 8.4), facecolor="white")
    grid = figure.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1, 1.18])
    axis_a = figure.add_subplot(grid[:, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    axis_c = figure.add_subplot(grid[1, 1])
    axis_a.set_axis_off()
    axis_a.set_xlim(0, 1)
    axis_a.set_ylim(0, 1)
    axis_a.set_title("A  PRINCE sample-to-analysis mapping", loc="left", fontsize=11, weight="bold", color=INK)
    y_positions = [0.86, 0.67, 0.48, 0.20]
    colors = [BLUE, BLUE, GREEN, ORANGE]
    labels = [
        f"SecAct-distributed expression matrix\n{int(prince.n_expression_samples)} samples",
        f"Matched expression–clinical records\n{int(prince.n_expression_clinical_overlap)} samples",
        f"Documented nivolumab expression records\n{int(prince.n_nivolumab_expression)} samples",
        "Primary analysis set\n38 baseline patients: 18 responders, 20 nonresponders",
    ]
    for index, (y, color, label) in enumerate(zip(y_positions, colors, labels)):
        flow_width = 0.46 if index == 2 else 0.72
        add_box(axis_a, 0.14, y - 0.08, flow_width, 0.14, f"Step {index + 1}", label, color)
        if index < len(y_positions) - 1:
            arrow_x = 0.38 if index == 2 else 0.50
            current_bottom = y - 0.08
            next_top = y_positions[index + 1] + 0.06
            add_arrow(axis_a, arrow_x, current_bottom - 0.012, arrow_x, next_top + 0.012)
    add_box(axis_a, 0.65, 0.39, 0.30, 0.13, "Excluded n=7", "1 NOT EVALUABLE;\n6 missing response labels", ORANGE)
    add_arrow(axis_a, 0.60, 0.45, 0.635, 0.437)
    axis_a.text(0.5, 0.04, "A1: 25 expression records, 20 mapped (10 R/10 NR); C2: 20 records, 18 mapped (8 R/10 NR).\nOnly baseline tumor RNA with a known best overall response was retained; see Table S3 for descriptive included-versus-excluded comparison.",
                ha="center", va="bottom", fontsize=7.5, color=MUTED, wrap=True)

    axis_b.set_title("B  PDAC ICI-containing cohorts have distinct evidentiary roles", loc="left", fontsize=11, weight="bold", color=INK)
    pdac = external.loc[external.cancer.eq("PDAC")].copy()
    positions = np.arange(len(pdac))[::-1]
    totals = pdac["n"].to_numpy()
    responders = pdac["responders"].to_numpy()
    axis_b.barh(positions, totals, color="#DCE6EA", height=0.55, label="Nonresponders")
    axis_b.barh(positions, responders, color=ORANGE, height=0.55, label="Responders")
    for position, (_, row) in zip(positions, pdac.iterrows()):
        direction = "same direction" if row.direction_matches_PRINCE else "opposite estimate"
        axis_b.text(row.n + 0.35, position, f"n={int(row.n)}; {int(row.responders)} R / {int(row.nonresponders)} NR\n{direction}", va="center", fontsize=7.4, color=INK)
    axis_b.set_yticks(positions)
    axis_b.set_yticklabels(["GSE179351\nSensitivity", "GSE248014\nAssociation transportability"], fontsize=8)
    axis_b.set_xlim(0, 28)
    axis_b.set_xlabel("Patients with baseline RNA and response labels")
    axis_b.legend(frameon=False, fontsize=7.4, loc="lower right")
    axis_b.grid(axis="x", alpha=.15)

    axis_c.set_axis_off()
    axis_c.set_xlim(0, 1)
    axis_c.set_ylim(0, 1)
    axis_c.set_title("C  Cross-cancer cohorts are boundary analyses, not PDAC validation", loc="left", fontsize=11, weight="bold", color=INK)
    melanoma = external.loc[external.cancer.eq("Melanoma")]
    add_box(axis_c, 0.06, 0.28, 0.40, 0.34, "GSE78220", f"Melanoma\nn={int(melanoma.loc[melanoma.dataset.eq('GSE78220'), 'n'].iloc[0])}; anti-PD-1\nCross-cancer boundary", SKY)
    add_box(axis_c, 0.54, 0.28, 0.40, 0.34, "GSE91061", f"Melanoma\nn={int(melanoma.loc[melanoma.dataset.eq('GSE91061'), 'n'].iloc[0])}; anti-PD-1-based\nCross-cancer boundary", PINK)
    axis_c.text(0.5, 0.08, "Neither cohort can validate or refute the PDAC association. They test whether the score generalizes indiscriminately across cancer contexts.",
                ha="center", fontsize=7.6, color=MUTED, wrap=True)
    figure.suptitle("Supplementary Figure S2. PRINCE selection flow and restricted association-transportability roles", x=.01, ha="left", fontsize=13, weight="bold", color=INK)
    figure.subplots_adjust(left=.07, right=.96, top=.88, bottom=.08, hspace=.52, wspace=.30)
    save_figure(figure, "Supplementary_FigureS2_cohort_mapping_and_evidence_roles")


def make_figure_s3(data):
    doublet_154 = data["doublet_154"].groupby("scDblFinder.class", as_index=False)["cells"].sum()
    doublet_156 = data["doublet_156"].rename(columns={"Var1": "scDblFinder.class", "Freq": "cells"})
    qc = pd.concat(
        [doublet_154.assign(dataset="GSE154778"), doublet_156.assign(dataset="GSE156405")],
        ignore_index=True,
    )
    states = data["figure4"].loc[data["figure4"].evidence_layer.eq("single_cell_myeloid_state")].copy()
    states["label"] = states.feature.str.replace("_", " ")
    qc.to_csv(SOURCE_DATA / "Supplementary_FigureS3_doublet_qc.tsv", sep="\t", index=False)
    states.to_csv(SOURCE_DATA / "Supplementary_FigureS3_myeloid_state_summary.tsv", sep="\t", index=False)

    figure, axes = plt.subplots(1, 3, figsize=(13.8, 4.5), facecolor="white")
    axis_a, axis_b, axis_c = axes
    classes = ["singlet", "doublet"]
    positions = np.arange(2)
    width = .34
    for offset, dataset in zip([-width / 2, width / 2], ["GSE154778", "GSE156405"]):
        subset = qc.loc[qc.dataset.eq(dataset)].set_index("scDblFinder.class").reindex(classes).fillna(0)
        axis_a.bar(positions + offset, subset.cells, width, label=dataset, color=BLUE if dataset == "GSE154778" else GREEN)
        for x, value in zip(positions + offset, subset.cells):
            axis_a.text(x, value + max(qc.cells) * .025, f"{int(value)}", ha="center", fontsize=7.3)
    axis_a.set_xticks(positions)
    axis_a.set_xticklabels(classes)
    axis_a.set_ylabel("Cells")
    axis_a.set_title("A  Doublet classification before myeloid display", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_a.legend(frameon=False, fontsize=7.5)
    axis_a.grid(axis="y", alpha=.15)

    ordered = states.sort_values("effect", ascending=True)
    axis_b.barh(np.arange(len(ordered)), ordered["n"], color=[GOLD, BLUE, ORANGE])
    for y, (_, row) in enumerate(ordered.iterrows()):
        axis_b.text(row.n + 10, y, f"n={int(row.n)}", va="center", fontsize=7.4)
    axis_b.set_yticks(np.arange(len(ordered)))
    axis_b.set_yticklabels(ordered["label"], fontsize=7.6)
    axis_b.set_xlabel("Retained cells")
    axis_b.set_title("B  Retained GSE154778 myeloid state candidates", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_b.grid(axis="x", alpha=.15)

    ordered = states.sort_values("effect", ascending=True)
    axis_c.barh(np.arange(len(ordered)), ordered.effect, color=[GOLD, BLUE, ORANGE])
    for y, (_, row) in enumerate(ordered.iterrows()):
        axis_c.text(row.effect + .04, y, f"{row.effect:.2f}", va="center", fontsize=7.4)
    axis_c.set_yticks(np.arange(len(ordered)))
    axis_c.set_yticklabels(ordered["label"], fontsize=7.6)
    axis_c.set_xlabel("Mean inflammatory state score")
    axis_c.set_title("C  Inflammatory state annotation score", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_c.grid(axis="x", alpha=.15)
    figure.text(.5, .01, "QC and state annotation identify an inflammatory myeloid context. They do not localize complete Mito3 to tumor or myeloid cells and do not demonstrate a functional death program.",
                ha="center", fontsize=7.5, color=MUTED)
    figure.suptitle("Supplementary Figure S3. Single-cell doublet control and myeloid-state annotation", x=.01, ha="left", fontsize=13, weight="bold", color=INK)
    figure.tight_layout(rect=(0, .05, 1, .88))
    save_figure(figure, "Supplementary_FigureS3_single_cell_qc_and_myeloid_annotation")


def make_figure_s4(data):
    summary = data["geomx_summary"].iloc[0]
    coverage = data["geomx_coverage"].copy()
    coverage["label"] = coverage.module.str.replace("_", " ")
    coverage = coverage.sort_values("coverage", ascending=True)
    coverage.to_csv(SOURCE_DATA / "Supplementary_FigureS4_GeoMx_module_coverage.tsv", sep="\t", index=False)
    figure = plt.figure(figsize=(13.8, 5.1), facecolor="white")
    grid = figure.add_gridspec(1, 3, width_ratios=[1.35, .72, 1.12])
    axis_a = figure.add_subplot(grid[0, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    axis_c = figure.add_subplot(grid[0, 2])
    colors = [ORANGE if value < .5 else GOLD if value < 1 else GREEN for value in coverage.coverage]
    positions = np.arange(len(coverage))
    axis_a.barh(positions, coverage.coverage, color=colors, height=.60)
    for y, (_, row) in enumerate(coverage.iterrows()):
        axis_a.text(min(row.coverage + .035, 1.02), y, f"{int(row.n_present)}/{int(row.n_prespecified)}", va="center", fontsize=7.2)
    axis_a.set_yticks(positions)
    axis_a.set_yticklabels(coverage.label, fontsize=7.5)
    axis_a.set_xlim(0, 1.12)
    axis_a.set_xlabel("Available genes / predefined genes")
    axis_a.set_title("A  GeoMx module coverage", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_a.grid(axis="x", alpha=.15)

    axis_b.bar([0, 1], [int(summary.n_tumor_rois), int(summary.n_stroma_rois)], color=[ORANGE, BLUE], width=.62)
    axis_b.set_xticks([0, 1])
    axis_b.set_xticklabels(["Tumor\nROI", "Stroma\nROI"])
    for x, value in zip([0, 1], [int(summary.n_tumor_rois), int(summary.n_stroma_rois)]):
        axis_b.text(x, value + 4, str(value), ha="center", fontsize=8)
    axis_b.set_ylim(0, 135)
    axis_b.set_ylabel("ROIs")
    axis_b.set_title(
        f"B  ROI composition\n{int(summary.n_patients)} paired patients; {int(summary.n_matrix_rois)} ROIs",
        loc="left", fontsize=9.6, weight="bold", color=INK, linespacing=1.35,
    )
    axis_b.grid(axis="y", alpha=.15)

    axis_c.set_axis_off()
    axis_c.set_xlim(0, 1)
    axis_c.set_ylim(0, 1)
    axis_c.set_title("C  Interpretation boundaries", loc="left", fontsize=10.5, weight="bold", color=INK)
    add_box(axis_c, .08, .57, .84, .24, "BAX-only mitochondrial proxy", "Only BAX is available\n(1/3 final Mito3 genes); this is not\nMito3.", ORANGE)
    add_box(axis_c, .08, .27, .84, .22, "Myeloid proxy", "Only SPP1 is available\n(1/3 genes); it is not a complete\nmyeloid-ecology score.", BLUE)
    axis_c.text(.50, .09, "No geometric ROI coordinates or cell segmentation:\npaired-ROI differences do not prove spatial adjacency or cell contact.", ha="center", fontsize=7.0, color=MUTED)
    figure.suptitle("Supplementary Figure S4. GeoMx coverage audit and ROI-level spatial constraints", x=.01, ha="left", fontsize=13, weight="bold", color=INK)
    figure.tight_layout(rect=(0, .03, 1, .88))
    save_figure(figure, "Supplementary_FigureS4_GeoMx_coverage_and_roi_limits")


def make_figure_s5(data):
    figure5 = data["figure5"].loc[data["figure5"].evidence_layer.eq("TCGA_RPPA_chain")].copy()
    raw_and_combined = figure5.copy()
    raw_and_combined[["edge", "analysis"]] = raw_and_combined["feature"].str.split(
        " | ", n=1, expand=True, regex=False
    )
    extra_rows = []
    edge_map = {
        ("Akt_pS473", "Rictor_total"): "Rictor total → Akt pS473",
        ("p70S6K_pT389", "Akt_pS473"): "Akt pS473 → p70S6K pT389",
    }
    for label, frame in [
        ("Conditional adjustment", data["rppa_conditional"]),
        ("Composition-proxy adjustment", data["rppa_composition"]),
        ("Combined protein/composition adjustment", data["rppa_combined"]),
    ]:
        for _, row in frame.iterrows():
            edge = edge_map.get((row.target, row.predictor))
            if edge:
                extra_rows.append({"edge": edge, "analysis": label, "n": row.n, "effect": row.rho_partial, "fdr": row.fdr})
    sensitivity = pd.DataFrame(extra_rows)
    raw_and_combined.to_csv(SOURCE_DATA / "Supplementary_FigureS5_main_RPPA_edges.tsv", sep="\t", index=False)
    sensitivity.to_csv(SOURCE_DATA / "Supplementary_FigureS5_adjustment_sensitivity.tsv", sep="\t", index=False)
    figure = plt.figure(figsize=(13.8, 5.4), facecolor="white")
    grid = figure.add_gridspec(1, 3, width_ratios=[1.05, 1.20, .88])
    axis_a = figure.add_subplot(grid[0, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    axis_c = figure.add_subplot(grid[0, 2])
    main = raw_and_combined.loc[:, ["edge", "analysis", "effect", "fdr", "n"]].copy()
    edge_order = ["Rictor total → Akt pS473", "Akt pS473 → p70S6K pT389"]
    analysis_order = ["Unadjusted", "Protein + composition adjusted"]
    y_positions = {edge: position for position, edge in enumerate(edge_order)}
    for analysis, marker, color, offset in [(analysis_order[0], "o", BLUE, -.10), (analysis_order[1], "s", GREEN, .10)]:
        subset = main.loc[main.analysis.eq(analysis)]
        for _, row in subset.iterrows():
            y = y_positions[row.edge] + offset
            axis_a.plot([0, row.effect], [y, y], color=color, lw=1.6)
            axis_a.scatter(row.effect, y, marker=marker, s=48, color=color, edgecolor="white", zorder=3)
            axis_a.text(row.effect + .012, y, f"ρ={row.effect:.3f}\nFDR={row.fdr:.3g}", va="center", fontsize=6.8)
    axis_a.axvline(0, color=GRAY, lw=.8, ls="--")
    axis_a.set_yticks(range(len(edge_order)))
    axis_a.set_yticklabels(edge_order, fontsize=7.5)
    axis_a.set_xlim(-.05, .34)
    axis_a.set_xlabel("Association")
    axis_a.set_title("A  Raw versus combined adjustment", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_a.grid(axis="x", alpha=.15)

    for edge, color in zip(edge_order, [GRAY, GREEN]):
        subset = sensitivity.loc[sensitivity.edge.eq(edge)].copy()
        positions = np.arange(len(subset))
        axis_b.plot(positions, subset.effect, color=color, marker="o", lw=1.8, label=edge)
        for x, (_, row) in zip(positions, subset.iterrows()):
            axis_b.text(x, row.effect + .014, f"{row.effect:.2f}\nFDR={row.fdr:.2g}", ha="center", fontsize=6.6)
    axis_b.axhline(0, color=GRAY, lw=.8, ls="--")
    axis_b.set_xticks(np.arange(3))
    axis_b.set_xticklabels(["Conditional\n(n=122)", "Composition\n(n=115)", "Combined\n(n=113)"], fontsize=7.2)
    axis_b.set_ylim(0, .31)
    axis_b.set_ylabel("Partial Spearman rho")
    axis_b.set_title("B  Adjustment sensitivity", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_b.legend(frameon=False, fontsize=6.9, loc="upper right")
    axis_b.grid(axis="y", alpha=.15)

    axis_c.set_axis_off()
    axis_c.set_xlim(0, 1)
    axis_c.set_ylim(0, 1)
    axis_c.set_title("C  Assay and inference boundaries", loc="left", fontsize=10.5, weight="bold", color=INK)
    add_box(axis_c, .08, .58, .84, .22, "RPPA signal", "Akt_pS473 is a multiplexed\nAKT1/2/3 antibody feature, not\nAKT1-specific phosphosite quantification.", GOLD)
    add_box(axis_c, .08, .30, .84, .20, "Cohort context", "TCGA-PAAD RPPA lacks ICI treatment.\nAdjustment uses total-protein and RNA\ncomposition proxies, not direct purity.", BLUE)
    axis_c.text(.50, .10, "Only the downstream pS473–p70S6K edge is retained\nafter the combined adjustment. This is not a complete\nRICTOR activation chain.", ha="center", fontsize=6.9, color=MUTED)
    figure.suptitle("Supplementary Figure S5. Sensitivity of TCGA-PAAD RPPA pathway-edge associations to adjustment", x=.01, ha="left", fontsize=13, weight="bold", color=INK)
    figure.tight_layout(rect=(0, .03, 1, .88))
    save_figure(figure, "Supplementary_FigureS5_RPPA_sensitivity_and_assay_boundaries")


def make_figure_s6(data):
    screen = data["prism_screen"].copy()
    candidates = data["candidates"].copy()
    screen["neglog10_fdr"] = -np.log10(np.clip(screen.fdr.astype(float), 1e-300, 1))
    screen["passes_rule"] = (screen.delta_lfc_low_minus_high < -0.25) & (screen.fdr < 0.20)
    screen.to_csv(SOURCE_DATA / "Supplementary_FigureS6_PRISM_state_selectivity_screen.tsv", sep="\t", index=False)
    candidates.to_csv(SOURCE_DATA / "Supplementary_FigureS6_candidate_identity_and_ranking.tsv", sep="\t", index=False)
    figure = plt.figure(figsize=(13.8, 8.2), facecolor="white")
    grid = figure.add_gridspec(2, 2, width_ratios=[1.05, 1.05], height_ratios=[1, 1])
    axis_a = figure.add_subplot(grid[0, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    axis_c = figure.add_subplot(grid[1, 0])
    axis_d = figure.add_subplot(grid[1, 1])
    axis_a.scatter(screen.delta_lfc_low_minus_high, screen.neglog10_fdr, s=8, color="#9FB3BD", alpha=.45, linewidth=0)
    axis_a.axvline(-.25, color=ORANGE, ls="--", lw=1.0)
    axis_a.axhline(-np.log10(.20), color=ORANGE, ls="--", lw=1.0)
    axis_a.text(-.23, .66, "exploratory\nrule", color=ORANGE, fontsize=7.2, ha="left")
    axis_a.set_xlabel("Low-state minus high-state median PRISM LFC")
    axis_a.set_ylabel("−log10(BH FDR)")
    axis_a.set_title("A  Complete 4,686-compound state-selectivity screen", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_a.text(.97, .94, "0 compounds met both criteria", transform=axis_a.transAxes, ha="right", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=.25", fc="white", ec="#D7E0E4", lw=.8))
    axis_a.grid(alpha=.15)

    top = screen.nsmallest(12, "wilcox_p").sort_values("wilcox_p", ascending=True)
    positions = np.arange(len(top))[::-1]
    axis_b.hlines(positions, 0, -np.log10(top.wilcox_p), color=GRAY, lw=1.7)
    axis_b.scatter(-np.log10(top.wilcox_p), positions, color=ORANGE, s=28, zorder=3)
    axis_b.set_yticks(positions)
    axis_b.set_yticklabels(top.compound, fontsize=7.1)
    axis_b.set_xlabel("−log10(unadjusted Wilcoxon P)")
    axis_b.set_title("B  Nominal top hits remain FDR-null", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_b.text(.97, .04, "All displayed compounds: BH FDR=1.0", transform=axis_b.transAxes, ha="right", fontsize=7.2, color=MUTED)
    axis_b.grid(axis="x", alpha=.15)

    candidates["joint_rank"] = (candidates.unipert_rictor_rank + candidates.unipert_mtor_rank) / 2
    for _, row in candidates.iterrows():
        color = GRAY if row.compound == "SB-2343" else (ORANGE if row.compound == "GSK2126458" else GREEN)
        marker = "X" if row.compound == "SB-2343" else "o"
        axis_c.scatter(row.joint_rank, row.prism_median_lfc, s=72, marker=marker, color=color, edgecolor="white", zorder=3)
        axis_c.text(row.joint_rank + .3, row.prism_median_lfc + .06, row.compound, fontsize=7.1)
    axis_c.axvspan(0, 20, color=PALE_GREEN, alpha=.7, zorder=0)
    axis_c.set_xlim(0, 22)
    axis_c.set_ylim(-4.05, -1.0)
    axis_c.set_xlabel("Mean UniPert RICTOR/MTOR rank (lower is closer)")
    axis_c.set_ylabel("Median PRISM LFC across 37 PDAC lines")
    axis_c.set_title("C  Test-priority ranking is distinct from state selectivity", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis_c.grid(alpha=.15)

    axis_d.set_axis_off()
    axis_d.set_xlim(0, 1)
    axis_d.set_ylim(0, 1)
    axis_d.set_title("D  Translational boundary", loc="left", fontsize=10.5, weight="bold", color=INK)
    add_box(axis_d, .08, .58, .84, .20, "What the screen supports", "Broad PDAC viability and target-space\nevidence can prioritize a small\nexperimental test set.", GREEN)
    add_box(axis_d, .08, .31, .84, .20, "What the screen does not support", "No Mito-state-selective compound, target\nengagement, ICI synergy, or patient benefit\nis demonstrated.", ORANGE)
    axis_d.text(.5, .10, "SB-2343 is displayed in gray because its exact compound identity\nremains unresolved and it should not be purchased for validation.", ha="center", fontsize=6.9, color=MUTED)
    figure.suptitle("Supplementary Figure S6. Full PRISM state-selectivity screen and bounded PI3K/mTOR test prioritization", x=.01, ha="left", fontsize=13, weight="bold", color=INK)
    figure.tight_layout(rect=(0, .03, 1, .90))
    save_figure(figure, "Supplementary_FigureS6_PRISM_screen_and_drug_prioritization")


def write_legends():
    LEGENDS_PATH.write_text(
        """# Supplementary table and figure legends (v1)\n\n"
        "## Supplementary Table S1. Cohort and dataset inventory.\n\n"
        "Public datasets are grouped by evidentiary role, not treated as interchangeable validation cohorts. The table lists retained analytic units, treatment or endpoint, manuscript figures and the claim boundary for every dataset used in the main narrative.\n\n"
        "## Supplementary Figure S2. PRINCE selection flow and restricted association-transportability roles.\n\n"
        "(A) Sample-to-analysis mapping for the PRINCE source distribution. Of 45 nivolumab-exposed expression records, 38 baseline records had a binary best-overall-response mapping (18 responders and 20 nonresponders); seven were excluded because one was NOT EVALUABLE and six had missing response labels. A1 contributed 20 mapped records (10 responders and 10 nonresponders) and C2 contributed 18 (8 responders and 10 nonresponders). (B) PDAC ICI-containing external cohorts. GSE179351 is a six-patient directional-sensitivity cohort, whereas GSE248014 is a 23-patient response-association transportability audit with three responders. Neither cohort has a non-ICI comparator, so they do not test treatment-effect modification. (C) Melanoma cohorts are cross-cancer boundary analyses and are not PDAC replications.\n\n"
        "## Supplementary Figure S3. Single-cell doublet control and myeloid-state annotation.\n\n"
        "(A) scDblFinder classifications before the retained-myeloid display in GSE154778 and GSE156405. (B) Retained GSE154778 myeloid state-candidate cell counts. (C) Mean inflammatory-state scores. These results provide immune-context annotation only; they do not localize complete Mito3 or demonstrate a functional death program.\n\n"
        "## Supplementary Figure S4. GeoMx coverage audit and ROI-level spatial constraints.\n\n"
        "(A) Gene coverage for each pre-defined module in GSE240078. (B) ROI composition across 40 patients. (C) The mitochondrial and myeloid readouts are BAX-only and SPP1-only proxies, respectively. The release lacks geometric coordinates and cell segmentation; paired-ROI comparisons do not demonstrate cellular proximity or contact.\n\n"
        "## Supplementary Figure S5. Sensitivity of TCGA-PAAD RPPA pathway-edge associations to adjustment.\n\n"
        "(A) Unadjusted and combined protein/composition-adjusted pathway-edge associations. (B) Conditional, composition-proxy and combined adjustment sensitivity. (C) Assay boundaries: Akt_pS473 is a multiplexed AKT1/2/3 antibody feature and TCGA-PAAD lacks ICI treatment. Only the downstream Akt_pS473–p70S6K_pT389 association remains after the combined adjustment.\n\n"
        "## Supplementary Figure S6. Full PRISM state-selectivity screen and bounded PI3K/mTOR test prioritization.\n\n"
        "(A) The 4,686-compound screen across 31 PDAC models; no compound met the exploratory low-state-minus-high-state effect and BH-FDR criteria. (B) Nominal top compounds remain FDR-null. (C) UniPert-assisted target-space ranking combined with broad PRISM activity in 37 PDAC lines prioritizes experimental candidates but does not establish state-selective sensitivity. (D) Interpretation boundaries, including the unresolved exact identity of SB-2343.\n\n"
        "## Supplementary Figure S7. UMAP visualization and marker support for descriptive GSE154778 myeloid-state annotation.\n\n"
        "(A) UMAP of 1,868 doublet-filtered GSE154778 myeloid-candidate cells, generated from the locked PCA representation with a fixed random seed and without re-clustering. (B) Pre-specified marker support across retained state candidates. (C) UMAP feature maps for the same markers. These panels support annotation only, not lineage, Mito3 localization, death-program activity or tumor–myeloid interaction.\n\n"
        "## Supplementary Figure S8. Patient-level SecAct sender robustness and restricted myeloid-endpoint communication context.\n\n"
        "Six pre-specified SecAct pairs with monocyte, macrophage or dendritic endpoints are shown. (A) Sender-expression detection among patients with measurable sender-cell expression. (B) Overall strength and BH-adjusted P values; receiver activity and strength remain condition-level SecAct inferences. (C) CCL5-to-dendritic context occurs in both cohorts with different senders. The panels do not demonstrate direct ligand action, physical contact, conserved sender identity or causal Mito3 control.\n""",
        encoding="utf-8",
    )


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    TABLE_OUTPUT.mkdir(parents=True, exist_ok=True)
    configure_style()
    data = read_inputs()
    build_table_s1(data)
    make_figure_s2(data)
    make_figure_s3(data)
    make_figure_s4(data)
    make_figure_s5(data)
    make_figure_s6(data)
    write_legends()


if __name__ == "__main__":
    main()
