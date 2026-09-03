from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:/第九篇大论文")
EXPR = ROOT / "data/03_results/feasibility/integrated_pdac_stratification/pdac_cellline_expression_modules_dependency.tsv"
PRISM = ROOT / "data/03_results/feasibility/prism_19q4/pdac_compound_sensitivity_long.tsv"
OUT = ROOT / "data/03_results/feasibility/prism_state_selective_screen"
OUT.mkdir(parents=True, exist_ok=True)
MATRIX = ROOT / "public_data_audit/evidence_chain_output/evidence_chain_matrix.tsv"


def main():
    expr = pd.read_csv(EXPR, sep="\t")
    prism = pd.read_csv(PRISM, sep="\t")
    d = prism.merge(expr[["ModelID", "mitochondrial_oxidative", "innate_input", "myeloid_ecology"]],
                    left_on="cell_line", right_on="ModelID", how="inner")
    d["lfc_median"] = pd.to_numeric(d["lfc_median"], errors="coerce")
    d = d.dropna(subset=["lfc_median", "mitochondrial_oxidative"])
    threshold = d["mitochondrial_oxidative"].median()
    d["state"] = np.where(d["mitochondrial_oxidative"] <= threshold, "mito_low", "mito_high")
    rows = []
    for (broad_id, name, moa, target), block in d.groupby(["broad_id", "name", "moa", "target"], dropna=False):
        low = block.loc[block.state == "mito_low", "lfc_median"]
        high = block.loc[block.state == "mito_high", "lfc_median"]
        if len(low) < 10 or len(high) < 10:
            continue
        test = mannwhitneyu(low, high, alternative="two-sided")
        rows.append({
            "broad_id": broad_id,
            "compound": name,
            "moa": moa,
            "target": target,
            "n_mito_low": len(low),
            "n_mito_high": len(high),
            "median_lfc_mito_low": low.median(),
            "median_lfc_mito_high": high.median(),
            "delta_lfc_low_minus_high": low.median() - high.median(),
            "wilcox_p": test.pvalue,
        })
    result = pd.DataFrame(rows)
    result["fdr"] = multipletests(result.wilcox_p, method="fdr_bh")[1]
    result["target_text"] = (result["moa"].fillna("").astype(str) + " " + result["target"].fillna("").astype(str)).str.lower()
    result["pathway_relevance"] = np.where(result.target_text.str.contains("mtor|pi3k|akt|rps6kb1|bcl2|autophagy|rho|rock", regex=True), "prespecified_pathway_related", "other")
    result.sort_values(["fdr", "delta_lfc_low_minus_high"]).to_csv(OUT / "all_compounds_state_selectivity.tsv", sep="\t", index=False)
    candidates = result[(result.delta_lfc_low_minus_high < -0.25) & (result.fdr < 0.20)].copy()
    candidates.sort_values(["fdr", "delta_lfc_low_minus_high"]).to_csv(OUT / "candidate_state_selective_compounds.tsv", sep="\t", index=False)
    related = result[result.pathway_relevance == "prespecified_pathway_related"].sort_values(["fdr", "delta_lfc_low_minus_high"])
    related.to_csv(OUT / "prespecified_pathway_related_compounds.tsv", sep="\t", index=False)
    top = result.sort_values(["fdr", "delta_lfc_low_minus_high"]).head(10)
    lines = [
        "# PRISM mitochondrial-state-selective drug screen",
        "",
        f"The screen matched {d.ModelID.nunique()} PDAC cell models and {result.broad_id.nunique()} compounds. Cell lines were split at the median tumor `mitochondrial_oxidative` expression score; lower PRISM log-fold-change indicates stronger growth inhibition.",
        "",
        "## Result",
        "",
        f"{len(candidates)} compounds met the exploratory ranking rule (low-state minus high-state median LFC < -0.25 and FDR < 0.20). The top candidates and all prespecified-pathway-related compounds are reported in machine-readable tables.",
        "",
        "No compound from this screen can be called RICTOR-specific, ICI-synergistic, or clinically effective. The analysis is a cell-line state-association screen with a median split, not a randomized treatment comparison. FDR is reported across the complete tested compound set.",
        "",
        "## Top-ranked rows",
        "",
        top[["compound", "moa", "target", "delta_lfc_low_minus_high", "fdr"]].to_string(index=False),
        "",
        "## Outputs",
        "",
        "- `all_compounds_state_selectivity.tsv`",
        "- `candidate_state_selective_compounds.tsv`",
        "- `prespecified_pathway_related_compounds.tsv`",
    ]
    (OUT / "PRISM_state_selective_screen.md").write_text("\n".join(lines), encoding="utf-8")
    if MATRIX.exists():
        matrix = pd.read_csv(MATRIX, sep="\t", dtype=str).fillna("")
        candidate_names = "; ".join(candidates.sort_values("fdr").head(5)["compound"].astype(str)) or "none"
        row = {
            "evidence_domain": "functional_in_silico",
            "dataset": "PRISM 19Q4 expression+drug",
            "unit": "PDAC_cell_line",
            "n": str(int(d.ModelID.nunique())),
            "events": "NA",
            "candidate_or_module": "mitochondrial_oxidative state-selective screen",
            "contrast": "mito-low vs mito-high cell lines",
            "effect_direction": "low_state_more_sensitive_exploratory",
            "effect_size": f"n_candidates={len(candidates)}; top={candidate_names}",
            "ci_low": "NA",
            "ci_high": "NA",
            "p_value": "NA",
            "fdr": "see_output_table",
            "endpoint": "PRISM log-fold-change",
            "supports_chain_component": "candidate drug prioritization",
            "independent_or_reused": "existing project PRISM data",
            "confounder_adjustment": "median state split; compound-wise Wilcoxon with global BH FDR",
            "interpretation": "Exploratory compounds show preferential inhibition in low mitochondrial-oxidative PDAC cell lines.",
            "limitation": "Cell-line association; no RICTOR specificity, ICI combination, patient benefit, or causal mechanism.",
            "status": "exploratory_drug_screen",
            "responders": "NA",
            "nonresponders": "NA",
        }
        duplicate = (matrix.evidence_domain == row["evidence_domain"]) & (matrix.dataset == row["dataset"]) & (matrix.candidate_or_module == row["candidate_or_module"])
        if not duplicate.any():
            matrix = pd.concat([matrix, pd.DataFrame([row])], ignore_index=True)
            matrix.to_csv(MATRIX, sep="\t", index=False)
    print(f"models={d.ModelID.nunique()} compounds={result.broad_id.nunique()} candidates={len(candidates)}")
    print(top[["compound", "moa", "delta_lfc_low_minus_high", "fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()
