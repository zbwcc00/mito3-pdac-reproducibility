from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"D:/第九篇大论文")
RESULTS = ROOT / "data/03_results/feasibility"
FIGURES = ROOT / "manuscript/figures"
OUT = RESULTS / "figure4_mechanism_localization"

INK = "#20313B"
MUTED = "#667781"
BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
GOLD = "#E69F00"
PINK = "#CC79A7"
GRAY = "#8A99A3"
PALE = "#F7F9FA"


def configure_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
    })


def read_inputs():
    state_summary = pd.read_csv(
        RESULTS / "next_stage_nine_point_audit/myeloid_state_role_summary.tsv", sep="\t"
    )
    secact = pd.read_csv(
        RESULTS / "next_stage_nine_point_audit/secact_patient_level_sender_robustness.tsv", sep="\t"
    )
    region_scores = pd.read_csv(
        RESULTS / "gse240078/patient_region_module_scores.tsv", sep="\t"
    )
    paired_tests = pd.read_csv(
        RESULTS / "targeted_cptac_spatial/spatial_paired_tumor_stroma_tests.tsv", sep="\t"
    )
    coupling = pd.read_csv(
        RESULTS / "targeted_cptac_spatial/spatial_tumor_module_vs_stroma_myeloid.tsv", sep="\t"
    )
    return state_summary, secact, region_scores, paired_tests, coupling


def plot_myeloid_states(axis, state_summary):
    candidates = state_summary.loc[
        state_summary.annotation_role.eq("myeloid_state_candidate")
    ].copy()
    order = ["Inflammatory_SPP1_IL1B", "Macrophage_C1Q_APOE", "Dendritic_FCER1A"]
    candidates = candidates.loc[candidates.subcluster_label.isin(order)].copy()
    candidates["subcluster_label"] = pd.Categorical(
        candidates.subcluster_label, categories=order, ordered=True
    )
    candidates = candidates.sort_values("subcluster_label")
    columns = [
        "macrophage_C1Q_APOE_mean",
        "inflammatory_SPP1_IL1B_mean",
        "monocyte_FCN1_S100A8_mean",
        "dendritic_FCER1A_mean",
    ]
    display_names = ["C1Q/APOE\nmacrophage", "SPP1/IL1B\ninflammatory", "FCN1/S100A8\nmonocyte", "FCER1A\ndendritic"]
    values = candidates[columns].to_numpy(dtype=float)
    image = axis.imshow(values, cmap="YlOrRd", vmin=0, vmax=2.5, aspect="auto")
    for row_index, row_values in enumerate(values):
        for column_index, value in enumerate(row_values):
            color = "white" if value > 1.65 else INK
            axis.text(column_index, row_index, f"{value:.2f}", ha="center", va="center", fontsize=8, color=color)
    labels = []
    for _, row in candidates.iterrows():
        label = str(row.subcluster_label).replace("_", "\n")
        labels.append(f"{label}\n(n={int(row.n_cells)} cells)")
    axis.set_xticks(range(len(display_names)))
    axis.set_xticklabels(display_names, fontsize=7.7)
    axis.set_yticks(range(len(labels)))
    axis.set_yticklabels(labels, fontsize=7.6)
    axis.tick_params(length=0, pad=5)
    axis.set_title("A  Distinct myeloid states in GSE154778", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis.text(0, -0.28, "Scores are state annotations in doublet-filtered cells; they do not measure tumor-cell Mito3.",
              transform=axis.transAxes, ha="left", va="top", fontsize=7.5, color=MUTED)
    colorbar = plt.colorbar(image, ax=axis, fraction=.045, pad=.03)
    colorbar.set_label("Mean signature score", fontsize=7.5)
    colorbar.ax.tick_params(labelsize=7.5, length=2)
    return candidates


def plot_spatial_compartments(axis, region_scores, paired_tests):
    modules = [
        ("mitochondrial_oxidative", "BAX-only mitochondrial proxy\n(1/3 final Mito3 genes)", ORANGE),
        ("membrane_cytoskeleton", "Cytoskeleton\n(3/6 genes)", PINK),
        ("innate_input", "Innate input\n(6/8 genes)", BLUE),
    ]
    random = np.random.default_rng(20260828)
    tick_locations = []
    tick_labels = []
    for module_index, (module, module_label, color) in enumerate(modules):
        subset = region_scores.loc[region_scores.module.eq(module), ["patient_id", "region", "score"]].copy()
        paired = subset.pivot(index="patient_id", columns="region", values="score").dropna()
        x_stroma = module_index * 3.0
        x_tumor = x_stroma + .92
        for _, row in paired.iterrows():
            axis.plot([x_stroma, x_tumor], [row.stroma, row.tumor], color="#CBD5DB", lw=.55, alpha=.60, zorder=1)
        for x_value, region, marker in [(x_stroma, "stroma", "o"), (x_tumor, "tumor", "s")]:
            values = paired[region].to_numpy()
            jitter = random.normal(0, .048, len(values))
            axis.scatter(np.repeat(x_value, len(values)) + jitter, values, s=15, marker=marker,
                         color=color, alpha=.75, edgecolor="white", linewidth=.28, zorder=3)
            axis.scatter(x_value, np.median(values), s=52, marker="D", color=INK,
                         edgecolor="white", linewidth=.7, zorder=4)
        test = paired_tests.loc[paired_tests.module.eq(module)].iloc[0]
        axis.text((x_stroma + x_tumor) / 2, .97,
                  f"Δ={test.delta_tumor_minus_stroma:+.3f}\nFDR={test.fdr:.2g}",
                  transform=axis.get_xaxis_transform(),
                  ha="center", va="top", fontsize=7.5, color=color, weight="bold")
        tick_locations.extend([x_stroma, x_tumor])
        tick_labels.extend(["Stroma", "Tumor"])
        axis.text((x_stroma + x_tumor) / 2, -0.18, module_label, transform=axis.get_xaxis_transform(),
                  ha="center", va="top", fontsize=7.5, color=INK)
    axis.set_xticks(tick_locations)
    axis.set_xticklabels(tick_labels, fontsize=7.6)
    axis.set_ylabel("Patient-level ROI score")
    axis.set_title("B  Tumor–stroma compartmentalization in GSE240078", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis.text(0, -0.35, "Forty paired tumor/stroma ROI profiles; diamonds show medians. This GeoMx dataset has no geometric coordinates.",
              transform=axis.transAxes, ha="left", va="top", fontsize=7.5, color=MUTED)


def plot_spatial_coupling(axis, region_scores, coupling):
    tumor = region_scores.loc[
        region_scores.module.eq("innate_input") & region_scores.region.eq("tumor"),
        ["patient_id", "score"]
    ].rename(columns={"score": "tumor_innate"})
    stroma = region_scores.loc[
        region_scores.module.eq("myeloid_ecology") & region_scores.region.eq("stroma"),
        ["patient_id", "score"]
    ].rename(columns={"score": "stroma_myeloid"})
    paired = tumor.merge(stroma, on="patient_id", validate="one_to_one")
    axis.scatter(paired.tumor_innate, paired.stroma_myeloid, s=38, color=GOLD,
                 edgecolor="white", linewidth=.55, alpha=.85)
    coefficients = np.polyfit(paired.tumor_innate, paired.stroma_myeloid, 1)
    x_values = np.linspace(paired.tumor_innate.min(), paired.tumor_innate.max(), 100)
    axis.plot(x_values, np.polyval(coefficients, x_values), color=INK, lw=1.35)
    test = coupling.loc[coupling.tumor_module.eq("innate_input")].iloc[0]
    axis.text(.04, .96, f"ρ={test.spearman_r:.3f}\nFDR={test.fdr:.3g}\nn={int(test.n_patients)} paired patients",
              transform=axis.transAxes, ha="left", va="top", fontsize=8,
              bbox=dict(boxstyle="round,pad=.28", fc="white", ec="#D7E0E4", lw=.8))
    axis.set_xlabel("Tumor innate-input score (6/8 genes)")
    axis.set_ylabel("Stromal SPP1 proxy\n(myeloid ecology; 1/3 genes)")
    axis.set_title("C  Cross-compartment coupling is an ROI proxy", loc="left", fontsize=10.5, weight="bold", color=INK)
    axis.text(.5, -.30, "The inverse association is not evidence of direct cell contact or ligand–receptor signaling.",
              transform=axis.transAxes, ha="center", va="top", fontsize=7.5, color=MUTED)
    axis.grid(alpha=.16, lw=.6)


def plot_secact_context(axis, secact):
    receivers = ["Monocyte", "Macrophage", "Dendritic"]
    signals = secact.loc[secact.receiver.isin(receivers)].copy()
    signals["pair"] = signals.sender + " — " + signals.secreted_protein + " → " + signals.receiver
    signals["cohort_display"] = signals.cohort.replace({"GSE154778": "GSE154778", "GSE156405": "GSE156405"})
    signals = signals.sort_values(["cohort_display", "secact_overall_strength"], ascending=[True, False]).reset_index(drop=True)
    colors = signals.cohort.map({"GSE154778": BLUE, "GSE156405": GREEN})
    positions = np.arange(len(signals))[::-1]
    axis.hlines(positions, 0, np.log10(signals.secact_overall_strength), color="#D9E2E6", lw=1.0, zorder=1)
    axis.scatter(np.log10(signals.secact_overall_strength), positions, s=58, color=colors,
                 edgecolor="white", linewidth=.65, zorder=3)
    axis.set_yticks(positions)
    axis.set_yticklabels(signals.pair, fontsize=7.5)
    axis.tick_params(axis="y", length=0)
    axis.set_xlabel("SecAct overall strength (log10 scale)")
    axis.set_title("D  Patient-level SecAct communication context", loc="left", fontsize=10.5, weight="bold", color=INK)
    for y_value, (_, row) in zip(positions, signals.iterrows()):
        axis.text(1.02, y_value, f"{row.cohort}; n={int(row.n_patients_sender_measured)}",
                  transform=axis.get_yaxis_transform(), ha="left", va="center", fontsize=7.5,
                  color=BLUE if row.cohort == "GSE154778" else GREEN, clip_on=False)
    axis.text(.5, -.28, "Pre-specified pairs with a monocyte, macrophage, or dendritic endpoint; signals are cohort-specific support, not causal validation.",
              transform=axis.transAxes, ha="center", va="top", fontsize=7.5, color=MUTED, wrap=True)
    axis.grid(axis="x", alpha=.16, lw=.6)
    return signals


def write_evidence_summary(candidates, paired_tests, coupling, secact_signals):
    rows = []
    for _, row in candidates.iterrows():
        rows.append({
            "evidence_layer": "single_cell_myeloid_state",
            "cohort": "GSE154778",
            "feature": row.subcluster_label,
            "n": int(row.n_cells),
            "effect": row.inflammatory_SPP1_IL1B_mean,
            "p_value": np.nan,
            "fdr": np.nan,
            "interpretation": "state annotation only; not tumor-cell Mito3 localization",
        })
    for module in ["mitochondrial_oxidative", "membrane_cytoskeleton", "innate_input"]:
        row = paired_tests.loc[paired_tests.module.eq(module)].iloc[0]
        rows.append({
            "evidence_layer": "paired_spatial_roi",
            "cohort": "GSE240078",
            "feature": module,
            "n": int(row.n_patients),
            "effect": row.delta_tumor_minus_stroma,
            "p_value": row.paired_wilcox_p,
            "fdr": row.fdr,
            "interpretation": "paired tumor/stroma ROI compartmentalization; limited gene coverage",
        })
    row = coupling.loc[coupling.tumor_module.eq("innate_input")].iloc[0]
    rows.append({
        "evidence_layer": "cross_compartment_roi_coupling",
        "cohort": "GSE240078",
        "feature": "tumor_innate_input_vs_stromal_SPP1_proxy",
        "n": int(row.n_patients),
        "effect": row.spearman_r,
        "p_value": row.p,
        "fdr": row.fdr,
        "interpretation": "same-patient ROI coupling, not geometric proximity or direct contact",
    })
    for _, row in secact_signals.iterrows():
        rows.append({
            "evidence_layer": "patient_level_SecAct",
            "cohort": row.cohort,
            "feature": row.pair,
            "n": int(row.n_patients_sender_measured),
            "effect": row.secact_overall_strength,
            "p_value": np.nan,
            "fdr": row.secact_overall_padj,
            "interpretation": "cohort-specific communication context; not a causal tumor–myeloid mechanism",
        })
    pd.DataFrame(rows).to_csv(OUT / "Figure4_evidence_summary.tsv", sep="\t", index=False)


def write_results(candidates, paired_tests, coupling, secact_signals):
    inflammatory = candidates.loc[candidates.subcluster_label.eq("Inflammatory_SPP1_IL1B")].iloc[0]
    mito = paired_tests.loc[paired_tests.module.eq("mitochondrial_oxidative")].iloc[0]
    cytoskeleton = paired_tests.loc[paired_tests.module.eq("membrane_cytoskeleton")].iloc[0]
    innate = paired_tests.loc[paired_tests.module.eq("innate_input")].iloc[0]
    coupling_row = coupling.loc[coupling.tumor_module.eq("innate_input")].iloc[0]
    text = f"""# Results 4 and Figure 4 legend (v1)

## Results 4. Single-cell and paired-ROI analyses localize immune-context evidence but do not establish direct tumor–myeloid contact

Doublet-filtered GSE154778 cells resolved an `Inflammatory_SPP1_IL1B` myeloid state ({int(inflammatory.n_cells)} cells) alongside `Macrophage_C1Q_APOE` and `Dendritic_FCER1A` states. The inflammatory state had the highest mean inflammatory signature score among the retained state candidates ({inflammatory.inflammatory_SPP1_IL1B_mean:.2f}). This analysis identifies a distinct inflammatory myeloid context; it does not localize the three-gene `Mito3` score to tumor cells and does not imply mitoxyperilysis activation.

In GSE240078, 40 patients had paired tumor and stroma ROIs. Tumor ROI scores were higher than matched stroma for the BAX-only mitochondrial proxy (tumor-minus-stroma median difference, {mito.delta_tumor_minus_stroma:+.3f}; FDR={mito.fdr:.2g}) and for the partial cytoskeletal module ({cytoskeleton.delta_tumor_minus_stroma:+.3f}; FDR={cytoskeleton.fdr:.2g}), whereas the partial innate-input module was lower in tumor ROI ({innate.delta_tumor_minus_stroma:+.3f}; FDR={innate.fdr:.2g}). BAX is one of three final Mito3 genes, but this BAX-only measure is not Mito3. Because the GeoMx release has neither ROI coordinates nor full module coverage (BAX-only mitochondrial proxy, 1/3 final Mito3 genes; stromal myeloid ecology, 1/3 genes), these are compartmentalization observations rather than cell-level spatial localization.

Across the paired patients, tumor innate input was inversely associated with the matched stromal SPP1 myeloid proxy (Spearman rho={coupling_row.spearman_r:.3f}; FDR={coupling_row.fdr:.3g}). This same-patient cross-compartment relationship cannot establish physical proximity, ligand–receptor action, or a direction of signaling. Patient-level SecAct analysis yielded {len(secact_signals)} pre-specified communications with a myeloid or dendritic endpoint across GSE154778 and GSE156405. The sender identities differed between cohorts; these results are therefore presented as an auxiliary immune-context layer rather than a conserved tumor–myeloid mechanism.

## Figure 4 legend

### Figure 4. Single-cell and paired-ROI analyses localize immune-context evidence without establishing direct tumor–myeloid contact.

**A**, mean state-signature scores for retained doublet-filtered GSE154778 myeloid state candidates. **B**, paired tumor/stroma ROI scores from 40 GSE240078 patients; lines connect matched ROIs and diamonds denote medians. Module coverage is indicated explicitly. **C**, patient-level association between tumor innate input and matched stromal SPP1 myeloid proxy. **D**, patient-level SecAct overall strengths for pre-specified communications with monocyte, macrophage, or dendritic endpoints. All spatial readouts are ROI-based and lack geometric coordinate information; SecAct is an auxiliary communication-context analysis and not causal validation.
"""
    (ROOT / "manuscript/results_4_and_figure4_legend_v1.md").write_text(text, encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    configure_style()
    state_summary, secact, region_scores, paired_tests, coupling = read_inputs()
    figure = plt.figure(figsize=(15.2, 10.6), facecolor="white")
    grid = figure.add_gridspec(2, 2, width_ratios=[1.03, 1.30], height_ratios=[1, 1.05])
    axis_a = figure.add_subplot(grid[0, 0])
    axis_b = figure.add_subplot(grid[0, 1])
    axis_c = figure.add_subplot(grid[1, 0])
    axis_d = figure.add_subplot(grid[1, 1])
    candidates = plot_myeloid_states(axis_a, state_summary)
    plot_spatial_compartments(axis_b, region_scores, paired_tests)
    plot_spatial_coupling(axis_c, region_scores, coupling)
    secact_signals = plot_secact_context(axis_d, secact)
    figure.suptitle("Supplementary Figure S9. Immune-context annotations and paired-ROI associations", x=.01, ha="left", fontsize=14.0, weight="bold", color=INK)
    figure.subplots_adjust(left=.08, right=.89, top=.89, bottom=.10, hspace=.80, wspace=.72)
    for suffix in ("pdf", "png", "svg"):
        figure.savefig(FIGURES / f"Figure4_immune_context_localization.{suffix}", facecolor="white")
        figure.savefig(FIGURES / f"Supplementary_FigureS9_immune_context_localization.{suffix}", facecolor="white")
    plt.close(figure)
    write_evidence_summary(candidates, paired_tests, coupling, secact_signals)
    write_results(candidates, paired_tests, coupling, secact_signals)


if __name__ == "__main__":
    main()
