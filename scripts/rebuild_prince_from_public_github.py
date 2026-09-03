"""Rebuild PRINCE Mito3 analyses from the official public GitHub release."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
from statsmodels.formula.api import logit


ROOT = Path(r"D:/第九篇大论文")
COMMIT = "c01cf276b1cef27e61f7349bccfac37c0c1d6ab7"
SOURCE_ARCHIVE = ROOT / f"data/01_source_data/PRINCE_ParkerICI_{COMMIT}.zip"
SOURCE = ROOT / f"data/01_source_data/prince-trial-data-{COMMIT}"
OUTPUT = ROOT / "data/03_results/official_prince_public_reproduction_v1"
GENE_AUDIT = ROOT / "data/03_results/feasibility/pdac_mechanism_gene_set_audit.tsv"
CURRENT_PATIENT_SCORES = ROOT / "data/03_results/feasibility/frozen_mito3_sensitivity/PRINCE_Mito3_patient_scores.tsv"
CURRENT_STATISTICS = ROOT / "data/03_results/feasibility/frozen_mito3_sensitivity/Mito3_primary_response_statistics.tsv"
SEED = 20260828
BOOTSTRAPS = 10_000
PERMUTATIONS = 20_000
MITO3_GENES = ("BAX", "BAK1", "BID")
RESPONSE_MAP = {
    "COMPLETE RESPONSE": 1,
    "PARTIAL RESPONSE": 1,
    "STABLE DISEASE": 0,
    "PROGRESSIVE DISEASE": 0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rank_auc_lower_score_is_response(responders: np.ndarray, nonresponders: np.ndarray) -> float:
    return float(
        np.less.outer(responders, nonresponders).mean()
        + 0.5 * np.equal.outer(responders, nonresponders).mean()
    )


def responder_minus_nonresponder_median(frame: pd.DataFrame, score_column: str) -> float:
    responder_scores = frame.loc[frame["responder"].eq(1), score_column]
    nonresponder_scores = frame.loc[frame["responder"].eq(0), score_column]
    return float(responder_scores.median() - nonresponder_scores.median())


def bootstrap_intervals(frame: pd.DataFrame, score_column: str) -> dict[str, float]:
    responder_scores = frame.loc[frame["responder"].eq(1), score_column].to_numpy()
    nonresponder_scores = frame.loc[frame["responder"].eq(0), score_column].to_numpy()
    generator = np.random.default_rng(SEED + 1)
    differences = np.empty(BOOTSTRAPS, dtype=float)
    auc_values = np.empty(BOOTSTRAPS, dtype=float)
    for index in range(BOOTSTRAPS):
        responder_resample = generator.choice(responder_scores, size=len(responder_scores), replace=True)
        nonresponder_resample = generator.choice(nonresponder_scores, size=len(nonresponder_scores), replace=True)
        differences[index] = float(np.median(responder_resample) - np.median(nonresponder_resample))
        auc_values[index] = rank_auc_lower_score_is_response(responder_resample, nonresponder_resample)
    difference_low, difference_high = np.quantile(differences, [0.025, 0.975])
    auc_low, auc_high = np.quantile(auc_values, [0.025, 0.975])
    return {
        "bootstrap_median_difference_percentile_interval_low": float(difference_low),
        "bootstrap_median_difference_percentile_interval_high": float(difference_high),
        "bootstrap_rank_auc_percentile_interval_low": float(auc_low),
        "bootstrap_rank_auc_percentile_interval_high": float(auc_high),
        "bootstrap_direction_probability_lower_in_responders": float(np.mean(differences < 0)),
    }


def arm_stratified_permutation_pvalue(frame: pd.DataFrame, score_column: str) -> float:
    observed = responder_minus_nonresponder_median(frame, score_column)
    generator = np.random.default_rng(SEED)
    labels = frame["responder"].to_numpy(dtype=int)
    arms = frame["arm"].to_numpy(dtype=str)
    values = frame[score_column].to_numpy(dtype=float)
    null_statistics = np.empty(PERMUTATIONS, dtype=float)
    arm_levels = frame["arm"].astype(str).unique()
    for index in range(PERMUTATIONS):
        shuffled_labels = labels.copy()
        for arm in arm_levels:
            arm_indices = np.flatnonzero(arms == arm)
            shuffled_labels[arm_indices] = generator.permutation(shuffled_labels[arm_indices])
        null_statistics[index] = float(np.median(values[shuffled_labels == 1]) - np.median(values[shuffled_labels == 0]))
    return float((np.sum(np.abs(null_statistics) >= abs(observed)) + 1) / (PERMUTATIONS + 1))


def benjamini_hochberg(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    return np.minimum(adjusted, 1.0)


def read_gene_tpm_values(metadata: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    records: list[dict[str, float | str]] = []
    for metadata_row in metadata.itertuples(index=False):
        report_path = SOURCE / "RNAseq/Data" / getattr(metadata_row, "Filename")
        report = pd.read_csv(report_path, sep="\t", usecols=["Gene Symbol", "TPM"])
        report["Gene Symbol"] = report["Gene Symbol"].str.upper().str.strip()
        selected = report.set_index("Gene Symbol")["TPM"].reindex(genes)
        if selected.isna().any():
            missing = ";".join(selected.index[selected.isna()])
            raise ValueError(f"Missing required gene(s) in {report_path.name}: {missing}")
        record: dict[str, float | str] = {"deidentified_id": str(getattr(metadata_row, "_0"))}
        record.update(selected.astype(float).to_dict())
        records.append(record)
    return pd.DataFrame(records)


def main() -> None:
    if not SOURCE_ARCHIVE.is_file() or not SOURCE.is_dir():
        raise FileNotFoundError("The fixed official PRINCE archive and extracted source directory are required.")
    OUTPUT.mkdir(parents=True, exist_ok=True)

    clinical = pd.read_csv(SOURCE / "Clinical/PICI0002_ph2_clinical.csv", dtype={"Deidentified.ID": str})
    metadata = pd.read_csv(SOURCE / "RNAseq/NatureMed_GX_ph2_metadata.csv", dtype={"Deidentified.ID": str})
    if len(metadata) != 65 or metadata["Filename"].nunique() != 65:
        raise ValueError("Unexpected PRINCE RNA metadata cardinality.")

    gene_audit = pd.read_csv(GENE_AUDIT, sep="\t")
    gene_audit["gene_or_feature"] = gene_audit["gene_or_feature"].str.upper().str.strip()
    seven_modules = [
        "death_pathway_exclusion",
        "innate_inflammatory_context",
        "innate_input",
        "mTORC2_membrane_dynamics",
        "membrane_cytoskeleton",
        "mitochondrial_oxidative",
        "myeloid_ecology",
    ]
    historical_gene_sets = {
        module: gene_audit.loc[gene_audit["module"].eq(module), "gene_or_feature"].drop_duplicates().tolist()
        for module in seven_modules
    }
    corrected_gene_sets = historical_gene_sets.copy()
    corrected_gene_sets["mitochondrial_oxidative"] = list(MITO3_GENES)
    all_required_genes = sorted(set(gene for genes in historical_gene_sets.values() for gene in genes))
    expression = read_gene_tpm_values(metadata, all_required_genes)
    log2_tpm = np.log2(expression[all_required_genes].to_numpy(dtype=float) + 1)
    standardized = (log2_tpm - log2_tpm.mean(axis=0)) / log2_tpm.std(axis=0, ddof=1)
    standardized_expression = pd.DataFrame(standardized, columns=all_required_genes)
    standardized_expression.insert(0, "deidentified_id", expression["deidentified_id"])

    for module, genes in historical_gene_sets.items():
        standardized_expression[f"historical_{module}"] = standardized_expression[genes].mean(axis=1)
    for module, genes in corrected_gene_sets.items():
        standardized_expression[f"corrected_{module}"] = standardized_expression[genes].mean(axis=1)
    standardized_expression["Mito3"] = standardized_expression[list(MITO3_GENES)].mean(axis=1)

    metadata_renamed = metadata.rename(columns={"Deidentified.ID": "deidentified_id"})
    clinical_renamed = clinical.rename(columns={"Deidentified.ID": "deidentified_id"})
    scored = standardized_expression.merge(metadata_renamed, on="deidentified_id", how="left", validate="one_to_one")
    scored = scored.merge(clinical_renamed, on="deidentified_id", how="left", validate="one_to_one")
    scored["responder"] = scored["Best Overall Response"].map(RESPONSE_MAP)
    scored["arm"] = scored["Actual Arm"].astype(str)
    scored["received_nivolumab"] = scored["Received Nivolumab"].astype(str)
    retained = scored.loc[scored["received_nivolumab"].eq("Y") & scored["responder"].notna()].copy()
    retained["responder"] = retained["responder"].astype(int)
    if (len(retained), int(retained["responder"].sum())) != (38, 18):
        raise ValueError("Official public source did not reproduce the expected PRINCE response cohort.")

    historical_rows: list[dict[str, float | int | str]] = []
    corrected_rows: list[dict[str, float | int | str]] = []
    for module in seven_modules:
        for score_name, score_column, result_rows in (
            ("historical", f"historical_{module}", historical_rows),
            ("corrected", f"corrected_{module}", corrected_rows),
        ):
            responder_scores = retained.loc[retained["responder"].eq(1), score_column].to_numpy()
            nonresponder_scores = retained.loc[retained["responder"].eq(0), score_column].to_numpy()
            result_rows.append({
                "screen": score_name,
                "module": module,
                "gene_set": ";".join(corrected_gene_sets[module] if score_name == "corrected" else historical_gene_sets[module]),
                "n": len(retained),
                "responders": int(retained["responder"].sum()),
                "nonresponders": int((1 - retained["responder"]).sum()),
                "median_responder": float(np.median(responder_scores)),
                "median_nonresponder": float(np.median(nonresponder_scores)),
                "median_difference": float(np.median(responder_scores) - np.median(nonresponder_scores)),
                "wilcoxon_p": float(mannwhitneyu(responder_scores, nonresponder_scores, alternative="two-sided", method="asymptotic").pvalue),
                "rank_auc_lower_is_response": rank_auc_lower_score_is_response(responder_scores, nonresponder_scores),
            })
    historical_screen = pd.DataFrame(historical_rows)
    corrected_screen = pd.DataFrame(corrected_rows)
    historical_screen["BH_FDR"] = benjamini_hochberg(historical_screen["wilcoxon_p"])
    corrected_screen["BH_FDR"] = benjamini_hochberg(corrected_screen["wilcoxon_p"])
    historical_screen.to_csv(OUTPUT / "historical_seven_module_screen.tsv", sep="\t", index=False)
    corrected_screen.to_csv(OUTPUT / "corrected_seven_module_screen.tsv", sep="\t", index=False)

    responder_scores = retained.loc[retained["responder"].eq(1), "Mito3"].to_numpy()
    nonresponder_scores = retained.loc[retained["responder"].eq(0), "Mito3"].to_numpy()
    standard_deviation = retained["Mito3"].std(ddof=0)
    retained["Mito3_z"] = (retained["Mito3"] - retained["Mito3"].mean()) / standard_deviation
    model = logit("responder ~ Mito3_z + C(arm)", data=retained).fit(disp=False)
    coefficient_interval = model.conf_int().loc["Mito3_z"]
    primary = {
        "n": len(retained),
        "responders": int(retained["responder"].sum()),
        "nonresponders": int((1 - retained["responder"]).sum()),
        "median_difference": responder_minus_nonresponder_median(retained, "Mito3"),
        "wilcoxon_p": float(mannwhitneyu(responder_scores, nonresponder_scores, alternative="two-sided", method="asymptotic").pvalue),
        "rank_auc_lower_is_response": rank_auc_lower_score_is_response(responder_scores, nonresponder_scores),
        "rank_biserial_correlation": float(2 * rank_auc_lower_score_is_response(responder_scores, nonresponder_scores) - 1),
        "bootstrap_iterations": BOOTSTRAPS,
        "arm_stratified_permutation_p": arm_stratified_permutation_pvalue(retained, "Mito3"),
        "permutations": PERMUTATIONS,
        "random_seed": SEED,
        "adjusted_odds_ratio_per_sd": float(np.exp(model.params["Mito3_z"])),
        "adjusted_odds_ratio_ci_low": float(np.exp(coefficient_interval.iloc[0])),
        "adjusted_odds_ratio_ci_high": float(np.exp(coefficient_interval.iloc[1])),
        "adjusted_logit_p": float(model.pvalues["Mito3_z"]),
    }
    primary.update(bootstrap_intervals(retained, "Mito3"))
    primary_statistics = pd.DataFrame([primary])
    primary_statistics.to_csv(OUTPUT / "official_public_Mito3_primary_statistics.tsv", sep="\t", index=False)

    current_scores = pd.read_csv(CURRENT_PATIENT_SCORES, sep="\t", dtype={"deidentified_id": str})
    current_scores["deidentified_id"] = current_scores["deidentified_id"].astype(str)
    score_comparison = retained[["deidentified_id", "responder", "arm", "Mito3"]].merge(
        current_scores[["deidentified_id", "responder", "arm", "Mito3"]],
        on="deidentified_id",
        how="outer",
        suffixes=("_official", "_current"),
        indicator=True,
    )
    score_comparison["Mito3_difference"] = score_comparison["Mito3_official"] - score_comparison["Mito3_current"]
    score_comparison.to_csv(OUTPUT / "official_vs_current_patient_score_comparison.tsv", sep="\t", index=False)
    current_statistics = pd.read_csv(CURRENT_STATISTICS, sep="\t").iloc[0].to_dict()
    comparable_metrics = [
        "n", "responders", "nonresponders", "median_difference", "wilcoxon_p",
        "rank_auc_lower_is_response", "rank_biserial_correlation",
        "bootstrap_median_difference_percentile_interval_low",
        "bootstrap_median_difference_percentile_interval_high",
        "bootstrap_rank_auc_percentile_interval_low",
        "bootstrap_rank_auc_percentile_interval_high",
        "bootstrap_direction_probability_lower_in_responders",
        "arm_stratified_permutation_p", "adjusted_odds_ratio_per_sd",
        "adjusted_odds_ratio_ci_low", "adjusted_odds_ratio_ci_high", "adjusted_logit_p",
    ]
    metric_comparison = pd.DataFrame([
        {
            "metric": metric,
            "official_public": primary[metric],
            "current": current_statistics[metric],
            "difference": float(primary[metric] - current_statistics[metric]),
        }
        for metric in comparable_metrics
    ])
    metric_comparison.to_csv(OUTPUT / "official_vs_current_metric_comparison.tsv", sep="\t", index=False)

    source_manifest = pd.DataFrame([
        {
            "source": "GitHub archive",
            "repository": "https://github.com/ParkerICI/prince-trial-data",
            "commit": COMMIT,
            "relative_path": SOURCE_ARCHIVE.relative_to(ROOT).as_posix(),
            "bytes": SOURCE_ARCHIVE.stat().st_size,
            "sha256": sha256(SOURCE_ARCHIVE),
        },
        *[
            {
                "source": "GitHub source file",
                "repository": "https://github.com/ParkerICI/prince-trial-data",
                "commit": COMMIT,
                "relative_path": path.relative_to(SOURCE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in [SOURCE / "Clinical/PICI0002_ph2_clinical.csv", SOURCE / "RNAseq/NatureMed_GX_ph2_metadata.csv"]
        ],
    ])
    source_manifest.to_csv(OUTPUT / "official_source_manifest.tsv", sep="\t", index=False)

    retained.to_csv(OUTPUT / "official_public_Mito3_patient_scores.tsv", sep="\t", index=False)
    max_score_difference = float(score_comparison["Mito3_difference"].abs().max())
    complete_mapping_match = bool(score_comparison["_merge"].eq("both").all())
    label_match = bool(
        score_comparison["responder_official"].eq(score_comparison["responder_current"]).all()
        and score_comparison["arm_official"].eq(score_comparison["arm_current"]).all()
    )
    summary = {
        "official_repository": "https://github.com/ParkerICI/prince-trial-data",
        "commit": COMMIT,
        "archive_sha256": sha256(SOURCE_ARCHIVE),
        "source_population": {
            "clinical_records": int(len(clinical)),
            "baseline_RNA_records": int(len(metadata)),
            "nivolumab_exposed_RNA_records": int(scored["received_nivolumab"].eq("Y").sum()),
            "retained_response_records": int(len(retained)),
            "responders": int(retained["responder"].sum()),
            "nonresponders": int((1 - retained["responder"]).sum()),
        },
        "Mito3_definition": "mean of gene-wise z-scored log2(TPM + 1) for BAX, BAK1, and BID across 65 baseline RNA records",
        "primary_statistics": primary,
        "comparison_with_current": {
            "complete_sample_mapping_match": complete_mapping_match,
            "response_and_arm_match": label_match,
            "maximum_absolute_Mito3_difference": max_score_difference,
            "maximum_absolute_primary_metric_difference": float(metric_comparison["difference"].abs().max()),
        },
        "interpretation": "The public GitHub release reproduces the current PRINCE Mito3 cohort and statistics from independent per-sample RNA reports. This does not authorize rehosting the source files because the repository has no explicit data license.",
    }
    (OUTPUT / "official_public_reproduction_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUTPUT / "README.md").write_text(
        "# Official public PRINCE reproduction\n\n"
        "This directory rebuilds the PRINCE Mito3 analysis from the official ParkerICI GitHub release at the fixed commit recorded in `official_source_manifest.tsv`. "
        "It does not use the prior SecAct-distributed expression matrix. The official source is public, but it has no explicit repository license; this project therefore records provenance and aggregate results without rehosting the source files.\n\n"
        "Outputs include historical and corrected seven-module screens, the primary Mito3 statistics, patient-level comparison against the pre-existing local analysis, and input hashes.\n",
        encoding="utf-8",
    )
    print(primary_statistics.to_string(index=False))
    print(metric_comparison.to_string(index=False))


if __name__ == "__main__":
    main()
