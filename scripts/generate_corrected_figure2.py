"""Replace Figure 2 with the audited three-gene Mito3 analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(r"D:/第九篇大论文")
INPUT = ROOT / "data/03_results/feasibility/frozen_mito3_sensitivity/PRINCE_Mito3_patient_scores.tsv"
OUT = ROOT / "manuscript/figures"
RESULTS = ROOT / "data/03_results/feasibility/frozen_mito3_sensitivity"
STATISTICS = RESULTS / "Mito3_primary_response_statistics.tsv"

ORANGE = "#D55E00"
GRAY = "#5C6770"
BLUE = "#0072B2"


def save(figure, name):
    for extension in ("pdf", "png", "svg"):
        figure.savefig(OUT / f"{name}.{extension}", bbox_inches="tight", dpi=400 if extension == "png" else None)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 9.5, "axes.titlesize": 10.5,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    data = pd.read_csv(INPUT, sep="\t")
    data["response"] = data.responder.map({1: "Responder", 0: "Nonresponder"})
    statistics = pd.read_csv(STATISTICS, sep="\t").iloc[0]
    figure, ax_a = plt.subplots(figsize=(5.3, 3.4))
    palette = {"Responder": "#F1BDAA", "Nonresponder": "#CBD5DB"}
    sns.boxplot(data=data, x="response", y="Mito3", hue="response", order=["Responder", "Nonresponder"],
                palette=palette, width=.55, showfliers=False, linewidth=1, legend=False, ax=ax_a)
    sns.stripplot(data=data, x="response", y="Mito3", hue="response", order=["Responder", "Nonresponder"],
                  palette={"Responder": ORANGE, "Nonresponder": GRAY}, jitter=.13, size=4.8, legend=False, ax=ax_a)
    ax_a.axhline(0, color="#CBD5E1", lw=.8, zorder=0)
    ax_a.set_xlabel("")
    ax_a.set_ylabel("Mito3 (BAX/BAK1/BID; within-cohort score)")
    ax_a.set_title("Lower baseline Mito3 in PRINCE responders", loc="left")
    ax_a.text(.5, .98, f"n={int(statistics.responders)} vs {int(statistics.nonresponders)}; Δmedian={statistics.median_difference:.3f}\nWilcoxon P={statistics.wilcoxon_p:.4f}; arm-permutation P={statistics.arm_stratified_permutation_p:.4f}",
              ha="center", va="top", transform=ax_a.transAxes, fontsize=8.1,
              bbox=dict(boxstyle="round,pad=.30", fc="white", ec="#CBD5E1"))
    figure.tight_layout()
    save(figure, "Figure2_prince_response_association")
    plt.close(figure)


if __name__ == "__main__":
    main()
