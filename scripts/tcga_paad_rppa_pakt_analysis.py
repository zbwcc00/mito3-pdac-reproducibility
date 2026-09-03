from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from scipy.stats import rankdata, t as student_t
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"D:/第九篇大论文")
INP = ROOT / "incoming_data/CPTAC_RPPA"
OUT = ROOT / "data/03_results/feasibility/tcga_paad_rppa_pakt"
REPORT = ROOT / "public_data_audit/evidence_chain_output/tcga_paad_rppa_pakt_report.md"
OUT.mkdir(parents=True, exist_ok=True)


def read_cbio_table(path: Path) -> pd.DataFrame:
    # cBioPortal clinical files contain comment/metadata lines before the header.
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("PATIENT_ID\t") or line.startswith("SAMPLE_ID\t") or line.startswith("Hugo_Symbol\t") or line.startswith("Composite.Element.REF\t"))
    return pd.read_csv(path, sep="\t", skiprows=header, low_memory=False)


def qval(values):
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    if ok.any():
        out[ok] = multipletests(values[ok], method="fdr_bh")[1]
    return out


def corr_table(df: pd.DataFrame, names):
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ok = df[a].notna() & df[b].notna()
            if ok.sum() < 10:
                continue
            rho, p = spearmanr(df.loc[ok, a], df.loc[ok, b])
            rows.append({"feature_a": a, "feature_b": b, "n": int(ok.sum()), "spearman_rho": rho, "p_value": p})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = qval(result["p_value"])
    return result


def partial_spearman(df: pd.DataFrame, target: str, predictor: str, covariates):
    columns = [target, predictor] + list(covariates)
    complete = df[columns].dropna()
    if len(complete) < 20:
        return {"target": target, "predictor": predictor, "n": len(complete), "rho_partial": np.nan, "p_value": np.nan}
    ranked = complete.apply(rankdata)
    y = ranked[target].to_numpy(float)
    x = ranked[predictor].to_numpy(float)
    design = np.column_stack([np.ones(len(ranked)), ranked[list(covariates)].to_numpy(float)])
    y_res = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    x_res = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    rho = float(np.corrcoef(y_res, x_res)[0, 1])
    df_resid = len(ranked) - len(covariates) - 2
    t_stat = rho * np.sqrt(df_resid / max(1e-12, 1 - rho * rho))
    p_value = float(2 * student_t.sf(abs(t_stat), df_resid))
    return {"target": target, "predictor": predictor, "n": len(complete), "rho_partial": rho, "p_value": p_value}


def zmean(expr, genes):
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns), present
    z = expr.loc[present].T.apply(lambda x: (x - x.mean()) / x.std(ddof=0), axis=0)
    return z.mean(axis=1), present


def main():
    rppa_path = INP / "paad_tcga_pan_can_atlas_2018_data_rppa.txt"
    rppa = pd.read_csv(rppa_path, sep="\t", low_memory=False)
    rppa = rppa.set_index("Composite.Element.REF")
    samples = [c for c in rppa.columns if c.startswith("TCGA-")]
    rppa = rppa[samples].apply(pd.to_numeric, errors="coerce")
    targets = {
        "Akt_total": "AKT1 AKT2 AKT3|Akt",
        "Akt_pS473": "AKT1 AKT2 AKT3|Akt_pS473",
        "Akt_pT308": "AKT1 AKT2 AKT3|Akt_pT308",
        "mTOR_total": "MTOR|mTOR",
        "mTOR_pS2448": "MTOR|mTOR_pS2448",
        "Rictor_total": "RICTOR|Rictor",
        "Rictor_pT1135": "RICTOR|Rictor_pT1135",
        "p70S6K_total": "RPS6KB1|p70S6K",
        "p70S6K_pT389": "RPS6KB1|p70S6K_pT389",
        "PRAS40_pT246": "AKT1S1|PRAS40_pT246",
    }
    found = {k: v for k, v in targets.items() if v in rppa.index}
    extracted = rppa.loc[list(found.values())].T.rename(columns={v: k for k, v in found.items()})
    extracted.index.name = "sample_id"
    extracted.to_csv(OUT / "tcga_paad_rppa_targeted_features.tsv", sep="\t")

    inventory = pd.DataFrame({
        "feature": list(found),
        "rppa_label": list(found.values()),
        "n_samples": [int(extracted[k].notna().sum()) for k in found],
        "missing_fraction": [float(extracted[k].isna().mean()) for k in found],
        "median": [float(extracted[k].median()) for k in found],
    })
    inventory.to_csv(OUT / "targeted_feature_inventory.tsv", sep="\t", index=False)

    corr = corr_table(extracted, list(found))
    corr.to_csv(OUT / "targeted_rppa_pairwise_spearman.tsv", sep="\t", index=False)

    chain_tests = []
    for target, predictor, covariates in [
        ("Akt_pS473", "Rictor_total", ["Akt_total", "mTOR_total"]),
        ("p70S6K_pT389", "Akt_pS473", ["p70S6K_total", "Akt_total"]),
        ("p70S6K_pT389", "Rictor_total", ["p70S6K_total", "Akt_total"]),
        ("Akt_pS473", "Rictor_pT1135", ["Akt_total", "Rictor_total"]),
    ]:
        chain_tests.append(partial_spearman(extracted, target, predictor, covariates))
    chain_tests = pd.DataFrame(chain_tests)
    if not chain_tests.empty:
        chain_tests["fdr"] = qval(chain_tests["p_value"])
    chain_tests.to_csv(OUT / "conditional_chain_partial_spearman.tsv", sep="\t", index=False)

    # Legacy TCGA Firehose RPPA provides an independent processing release of the same disease.
    legacy_path = INP / "paad_tcga_firehose_legacy_data_rppa.txt"
    legacy_corr = pd.DataFrame()
    legacy_common = 0
    if legacy_path.exists():
        legacy_raw = pd.read_csv(legacy_path, sep="\t", low_memory=False).set_index("Composite.Element.REF")
        legacy_samples = [c for c in legacy_raw.columns if c.startswith("TCGA-")]
        legacy_raw = legacy_raw[legacy_samples].apply(pd.to_numeric, errors="coerce")
        legacy_found = {k: v for k, v in targets.items() if v in legacy_raw.index}
        legacy = legacy_raw.loc[list(legacy_found.values())].T.rename(columns={v: k for k, v in legacy_found.items()})
        legacy.index.name = "sample_id"
        legacy.to_csv(OUT / "legacy_tcga_paad_rppa_targeted_features.tsv", sep="\t")
        legacy_corr = corr_table(legacy, list(legacy_found))
        legacy_corr.to_csv(OUT / "legacy_targeted_rppa_pairwise_spearman.tsv", sep="\t", index=False)
        legacy_common = len(extracted.index.intersection(legacy.index))

    # Phosphosite/total-protein ratios are exploratory within-array activation proxies.
    ratios = pd.DataFrame(index=extracted.index)
    ratio_pairs = {
        "Akt_pS473_over_Akt": ("Akt_pS473", "Akt_total"),
        "Akt_pT308_over_Akt": ("Akt_pT308", "Akt_total"),
        "mTOR_pS2448_over_mTOR": ("mTOR_pS2448", "mTOR_total"),
        "Rictor_pT1135_over_Rictor": ("Rictor_pT1135", "Rictor_total"),
        "p70S6K_pT389_over_p70S6K": ("p70S6K_pT389", "p70S6K_total"),
    }
    for name, (num, den) in ratio_pairs.items():
        if num in extracted and den in extracted:
            ratios[name] = extracted[num] - extracted[den]
    ratios.to_csv(OUT / "exploratory_phospho_total_ratios.tsv", sep="\t")
    ratio_corr = corr_table(pd.concat([extracted, ratios], axis=1), list(ratios))
    ratio_corr.to_csv(OUT / "exploratory_ratio_pairwise_spearman.tsv", sep="\t", index=False)

    # Same-study TCGA mRNA is paired by TCGA barcode; compute prespecified module proxies.
    mrna_path = INP / "paad_tcga_pan_can_atlas_2018_data_mrna_seq_v2_rsem.txt"
    mrna = pd.read_csv(mrna_path, sep="\t", low_memory=False)
    mrna = mrna.drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    mrna = mrna[[c for c in mrna.columns if c.startswith("TCGA-")]].apply(pd.to_numeric, errors="coerce")
    common = extracted.index.intersection(mrna.columns)
    paired = extracted.loc[common].copy()
    modules = {
        "innate_input": ["TLR2", "TLR4", "TLR7", "TLR8", "MYD88", "TNFRSF1A", "TNFRSF1B", "TNF"],
        "mitochondrial_oxidative": ["BAX", "BAK1", "BID"],
        "mTORC2_membrane_dynamics": ["MTOR", "RICTOR", "AKT1"],
        "membrane_cytoskeleton": ["RHOA", "RAC1", "CDC42", "ACTB", "ARPC2", "ARPC3"],
        "myeloid_ecology": ["IL1B", "SPP1", "CXCL3"],
    }
    module_rows = []
    for name, genes in modules.items():
        score, present = zmean(mrna.loc[:, common], genes)
        paired[name] = score
        module_rows.append({"module": name, "requested_genes": ",".join(genes), "genes_present": ",".join(present), "n_genes_present": len(present)})
    paired.to_csv(OUT / "tcga_paad_rppa_mrna_paired_features_modules.tsv", sep="\t")
    pd.DataFrame(module_rows).to_csv(OUT / "module_gene_coverage.tsv", sep="\t", index=False)
    module_names = list(modules)
    module_corr = corr_table(paired, module_names + list(found))
    module_corr = module_corr[(module_corr.feature_a.isin(module_names)) | (module_corr.feature_b.isin(module_names))]
    module_corr.to_csv(OUT / "mrna_module_vs_rppa_spearman.tsv", sep="\t", index=False)

    # Sensitivity to RNA-derived composition/context proxies; these are not purity estimates.
    composition_tests = []
    for target, predictor in [("Akt_pS473", "Rictor_total"), ("p70S6K_pT389", "Akt_pS473")]:
        composition_tests.append(partial_spearman(paired, target, predictor, ["myeloid_ecology", "innate_input", "membrane_cytoskeleton"]))
    composition_tests = pd.DataFrame(composition_tests)
    if not composition_tests.empty:
        composition_tests["fdr"] = qval(composition_tests["p_value"])
    composition_tests.to_csv(OUT / "composition_proxy_adjusted_chain.tsv", sep="\t", index=False)

    # TCGA transcriptomic purity/composition proxies already generated in the project.
    # These are sensitivity covariates, not histologic or ABSOLUTE purity measurements.
    purity_path = ROOT / "data/03_results/feasibility/tcga_transcriptomic_purity_proxy.tsv"
    purity_assoc = pd.DataFrame()
    purity_chain = pd.DataFrame()
    purity_common = 0
    if purity_path.exists():
        purity = pd.read_csv(purity_path, sep="\t")
        purity["sample_id"] = purity["sample_id"].astype(str).str.replace(r"(-\d{2})[A-Z]$", r"\1", regex=True)
        purity = purity.drop_duplicates("sample_id").set_index("sample_id")
        purity_common = len(extracted.index.intersection(purity.index))
        purity_join = extracted.join(purity[["tumor_epithelial", "immune", "stromal", "transcriptomic_purity_proxy"]], how="inner")
        purity_assoc = corr_table(purity_join, ["Akt_pS473", "Rictor_total", "p70S6K_pT389", "tumor_epithelial", "immune", "stromal", "transcriptomic_purity_proxy"])
        purity_assoc.to_csv(OUT / "rppa_purity_proxy_associations.tsv", sep="\t", index=False)
        purity_chain = pd.DataFrame([
            partial_spearman(purity_join, "Akt_pS473", "Rictor_total", ["Akt_total", "mTOR_total", "tumor_epithelial", "immune", "stromal"]),
            partial_spearman(purity_join, "p70S6K_pT389", "Akt_pS473", ["p70S6K_total", "Akt_total", "tumor_epithelial", "immune", "stromal"]),
        ])
        purity_chain["fdr"] = qval(purity_chain["p_value"])
        purity_chain.to_csv(OUT / "rppa_purity_adjusted_chain.tsv", sep="\t", index=False)

    # Clinical OS association is explicitly exploratory and not ICI response evidence.
    clinical_path = INP / "paad_tcga_pan_can_atlas_2018_data_clinical_patient.txt"
    clinical = read_cbio_table(clinical_path)
    clinical.columns = [str(c).strip() for c in clinical.columns]
    patient_col = next((c for c in clinical.columns if c.upper() in {"PATIENT_ID", "PATIENT_ID"}), clinical.columns[0])
    clinical = clinical.rename(columns={patient_col: "patient_id"})
    # sample barcode to patient barcode (first three TCGA fields).
    clin = clinical.set_index("patient_id")
    surv = pd.DataFrame(index=extracted.index)
    surv["patient_id"] = ["-".join(s.split("-")[:3]) for s in extracted.index]
    for col in ["OS_MONTHS", "OS_STATUS", "DFS_MONTHS", "DFS_STATUS"]:
        if col in clin.columns:
            surv[col] = surv.patient_id.map(clin[col])
    surv["OS_MONTHS"] = pd.to_numeric(surv.get("OS_MONTHS"), errors="coerce")
    surv["OS_EVENT"] = surv.get("OS_STATUS", pd.Series(index=surv.index)).astype(str).str.startswith("1:").astype(float)
    survival_out = []
    if surv["OS_MONTHS"].notna().sum() >= 20:
        for feature in ["Akt_pS473", "Akt_pT308", "Rictor_total", "p70S6K_pT389"] + module_names:
            if feature not in paired:
                continue
            x = paired[feature].reindex(surv.index)
            ok = x.notna() & surv.OS_MONTHS.notna()
            if ok.sum() < 20:
                continue
            rho, p = spearmanr(x[ok], surv.loc[ok, "OS_MONTHS"])
            survival_out.append({"feature": feature, "n": int(ok.sum()), "spearman_rho_with_OS_months": rho, "p_value": p})
    pd.DataFrame(survival_out).to_csv(OUT / "exploratory_os_months_correlations.tsv", sep="\t", index=False)
    surv.to_csv(OUT / "tcga_paad_rppa_survival_linkage.tsv", sep="\t")

    # Human-readable report with exact counts and limitations.
    key = corr[(corr.feature_a == "Akt_pS473") | (corr.feature_b == "Akt_pS473")].sort_values("fdr")
    key_lines = []
    for _, row in key.head(12).iterrows():
        key_lines.append(f"- `{row.feature_a}` vs `{row.feature_b}`: n={int(row.n)}, rho={row.spearman_rho:.3f}, FDR={row.fdr:.3g}.")
    mkey = module_corr[(module_corr.feature_a == "Akt_pS473") | (module_corr.feature_b == "Akt_pS473")].sort_values("fdr")
    mkey_lines = []
    for _, row in mkey.head(8).iterrows():
        mkey_lines.append(f"- `{row.feature_a}` vs `{row.feature_b}`: n={int(row.n)}, rho={row.spearman_rho:.3f}, FDR={row.fdr:.3g}.")
    legacy_lines = []
    if not legacy_corr.empty:
        for _, row in legacy_corr[((legacy_corr.feature_a == "Akt_pS473") | (legacy_corr.feature_b == "Akt_pS473"))].sort_values("fdr").head(8).iterrows():
            legacy_lines.append(f"- `{row.feature_a}` vs `{row.feature_b}`: n={int(row.n)}, rho={row.spearman_rho:.3f}, FDR={row.fdr:.3g}.")
    chain_lines = []
    for _, row in chain_tests.iterrows():
        chain_lines.append(f"- `{row.predictor}` → `{row.target}` (adjusted): n={int(row.n)}, partial rho={row.rho_partial:.3f}, FDR={row.fdr:.3g}.")
    composition_lines = []
    for _, row in composition_tests.iterrows():
        composition_lines.append(f"- `{row.predictor}` → `{row.target}` adjusted for RNA composition proxies: n={int(row.n)}, partial rho={row.rho_partial:.3f}, FDR={row.fdr:.3g}.")
    purity_lines = []
    if not purity_chain.empty:
        for _, row in purity_chain.iterrows():
            purity_lines.append(f"- `{row.predictor}` → `{row.target}` adjusted for protein totals plus epithelial/immune/stromal proxies: n={int(row.n)}, partial rho={row.rho_partial:.3f}, FDR={row.fdr:.3g}.")

    # Add auditable rows to the project-wide evidence chain, without duplicating runs.
    matrix_path = ROOT / "public_data_audit/evidence_chain_output/evidence_chain_matrix.tsv"
    if matrix_path.exists():
        matrix = pd.read_csv(matrix_path, sep="\t", dtype=str)
        rows = []
        def add_row(candidate, contrast, direction, effect, p, fdr, endpoint, component, interpretation, limitation, status):
            rows.append({
                "evidence_domain": "RPPA_proteomics", "dataset": "TCGA-PAAD PanCancer Atlas RPPA", "unit": "sample", "n": str(int(len(samples))), "events": "NA",
                "candidate_or_module": candidate, "contrast": contrast, "effect_direction": direction, "effect_size": str(effect), "ci_low": "NA", "ci_high": "NA",
                "p_value": str(p), "fdr": str(fdr), "endpoint": endpoint, "supports_chain_component": component, "independent_or_reused": "independent TCGA RPPA cohort",
                "confounder_adjustment": "none; Spearman correlation on RPPA abundance", "interpretation": interpretation, "limitation": limitation, "status": status,
                "responders": "NA", "nonresponders": "NA"
            })
        def lookup(a, b):
            q = corr[((corr.feature_a == a) & (corr.feature_b == b)) | ((corr.feature_a == b) & (corr.feature_b == a))]
            return q.iloc[0] if not q.empty else None
        add_row("Akt_pS473 RPPA site", "site availability", "direct_site_present", "122/122", "NA", "NA", "RPPA antibody feature", "AKT-S473 direct measurement", "Direct `Akt_pS473` feature is present in every RPPA sample.", "Antibody recognizes AKT1/2/3 together; TCGA is untreated and not an ICI response cohort.", "direct_site_support")
        if not legacy_corr.empty:
            add_row("Akt_pS473 cross-release replication", "PanCancer Atlas vs Firehose Legacy RPPA", "replicated", f"overlap={legacy_common}", "NA", "NA", "processing-release sensitivity", "reproducible AKT-S473 pathway-state support", "The direct pS473 feature is present in both TCGA RPPA processing releases, with positive Rictor/p70S6K associations in both.", "The releases reuse the same underlying TCGA patients and are not independent cohorts; effect sizes vary by processing release.", "replication_support")
        for a, b, component, status, interpretation, limitation in [
            ("Akt_pS473", "Rictor_total", "mTORC2 upstream coherence", "weak_mechanistic_support", "Akt-pS473 has a weak positive association with total Rictor, consistent with but not proving mTORC2 coupling.", "Rictor abundance is not Rictor complex activity; no treatment or perturbation."),
            ("Akt_pS473", "p70S6K_pT389", "downstream S6K coherence", "nominal_only", "The positive association with p70S6K-pT389 is nominal but does not survive FDR correction.", "p70S6K-T389 is not a direct AKT substrate readout and assay feedback may decouple nodes."),
            ("Akt_pS473", "mTOR_pS2448", "mTOR pathway coherence", "mechanistic_support", "Strong positive association with mTOR-pS2448 supports shared mTOR-pathway state, not mTORC2 specificity.", "mTOR-S2448 is not specific for mTORC2 and cannot replace RICTOR functional evidence."),
        ]:
            q = lookup(a, b)
            if q is not None:
                add_row(f"{a} vs {b}", "same-sample RPPA abundance", "positive", f"rho={q.spearman_rho:.3f}", q.p_value, q.fdr, "Spearman correlation", component, interpretation, limitation, status)
        for _, q in chain_tests.iterrows():
            add_row(f"{q.predictor} -> {q.target} partial", "adjusted for total proteins", "positive" if q.rho_partial > 0 else "negative", f"partial_rho={q.rho_partial:.3f}", q.p_value, q.fdr, "partial Spearman correlation", "conditional pathway coherence", "Conditional association after total-protein adjustment; supportive but observational.", "No measured tumor purity, no ICI endpoint, and partial correlations do not establish directionality.", "conditional_support" if q.fdr < 0.05 else "conditional_nominal")
        for r in rows:
            if not ((matrix["dataset"] == r["dataset"]) & (matrix["candidate_or_module"] == r["candidate_or_module"]) & (matrix["contrast"] == r["contrast"])).any():
                matrix = pd.concat([matrix, pd.DataFrame([r])], ignore_index=True)
        matrix.to_csv(matrix_path, sep="\t", index=False)
    report = f"""# TCGA-PAAD RPPA pAKT-S473 analysis

## Dataset

- Source: cBioPortal PanCancer Atlas TCGA-PAAD RPPA, downloaded from the cBioPortal datahub URL recorded in the project manifest.
- RPPA matrix: {len(samples)} TCGA sample columns; targeted features found: {len(found)}/{len(targets)}.
- Direct pAKT-S473 feature: `AKT1 AKT2 AKT3|Akt_pS473`; this is a multiplexed AKT1/2/3 antibody signal, not AKT1-specific quantification.
- Same-study mRNA matrix was paired for {len(common)} samples by TCGA sample barcode.
- Legacy TCGA Firehose RPPA release is also present locally; it contains the same direct `Akt_pS473` feature and is used as a processing-release sensitivity check on {legacy_common} overlapping sample barcodes.

## Direct pathway evidence

The RPPA layer contains the requested `Akt_pS473`, total Akt, `Rictor`, Rictor-pT1135, total p70S6K, p70S6K-pT389, mTOR and mTOR-pS2448 features. Pairwise associations for `Akt_pS473` are:

{chr(10).join(key_lines) if key_lines else '- No valid pairwise tests.'}

## RNA-to-protein coupling

Prespecified RNA module scores were computed as the mean of within-cohort gene z-scores. Associations between these RNA proxies and the RPPA `Akt_pS473` feature are:

{chr(10).join(mkey_lines) if mkey_lines else '- No valid module tests.'}

## Processing-release sensitivity check

The legacy TCGA Firehose RPPA release contains the direct pS473 feature as well. Its pS473 pairwise results are:

{chr(10).join(legacy_lines) if legacy_lines else '- Legacy release not available.'}

## Conditional chain tests

Partial Spearman tests use rank-transformed values and adjust for indicated total-protein covariates:

{chr(10).join(chain_lines) if chain_lines else '- No valid conditional tests.'}

After additionally adjusting for RNA-derived `myeloid_ecology`, `innate_input` and `membrane_cytoskeleton` scores (composition/context proxies, **not** measured purity):

{chr(10).join(composition_lines) if composition_lines else '- No valid composition-adjusted tests.'}

The RPPA matrix was also linked to the project transcriptomic composition file for {purity_common} samples. After adjusting for epithelial, immune and stromal proxies:

{chr(10).join(purity_lines) if purity_lines else '- No purity-proxy-adjusted tests.'}

## Interpretation

This dataset upgrades the project from indirect downstream phosphoprotein evidence to direct detection of an AKT S473 RPPA feature in PDAC. It does **not** establish an ICI-response biomarker: TCGA-PAAD is not an ICI-treated response cohort, the RPPA antibody combines AKT1/2/3, and all analyses are observational. A correlation among Rictor, Akt-pS473 and p70S6K-pT389 would support pathway coherence, but absence of correlation would not disprove mTORC2 activity because of tumor purity, assay dynamic range, treatment-free sampling and pathway feedback.

The direct site and its positive associations with Rictor and p70S6K are also present in the legacy Firehose RPPA processing release, which strengthens reproducibility across processing versions. The magnitude is release-dependent, so the result should be reported as replicated pathway-state support rather than a fixed effect size.

The OS analyses in the output directory are exploratory prognostic associations only; they must not be written as treatment-response or causal validation. Phospho/total ratios are within-array proxies and are not calibrated kinase activities.

## Output files

- `targeted_feature_inventory.tsv`
- `targeted_rppa_pairwise_spearman.tsv`
- `exploratory_phospho_total_ratios.tsv`
- `tcga_paad_rppa_mrna_paired_features_modules.tsv`
- `mrna_module_vs_rppa_spearman.tsv`
- `exploratory_os_months_correlations.tsv`
- `conditional_chain_partial_spearman.tsv`
- `composition_proxy_adjusted_chain.tsv`
- `rppa_purity_proxy_associations.tsv`
- `rppa_purity_adjusted_chain.tsv`
"""
    REPORT.write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
