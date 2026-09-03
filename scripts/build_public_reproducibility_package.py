"""Build a public, non-participant-level reproducibility package for this study.

The package deliberately excludes PRINCE source files and participant-level
outputs. It links to the official public source, and preserves exact source
paths and SHA-256 hashes for every released file.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path


ROOT = Path(r"D:/第九篇大论文")
PACKAGE = ROOT / "public_reproducibility_package_v8"

INCLUDED_FILES = {
    "README.md": "_generated_readme_",
    "data/PUBLIC_DATASET_MANIFEST.tsv": "PDAC_public_dataset_manifest.tsv",
    "documentation/mito3_provenance_and_lock_audit.md": "manuscript/mito3_provenance_timeline_audit_v1.md",
    "documentation/official_prince_public_reproduction_audit_v1.md": "manuscript/prince_public_source_reproduction_audit_v1.md",
    "documentation/citation_and_source_verification.md": "manuscript/citation_and_data_source_verification_report_v1.md",
    "documentation/evidence_ledger.md": "PDAC_mitoxyperilysis_evidence_ledger_2026-08-26.md",
    "documentation/evidence_chain_matrix.tsv": "public_data_audit/evidence_chain_output/evidence_chain_matrix.tsv",
    "scripts/gse240078_spatial_audit.R": "scripts/gse240078_spatial_audit.R",
    "scripts/tcga_paad_rppa_pakt_analysis.py": "scripts/tcga_paad_rppa_pakt_analysis.py",
    "scripts/prism_state_selective_drug_screen.py": "scripts/prism_state_selective_drug_screen.py",
    "scripts/generate_figure4_mechanism_localization.py": "manuscript/generate_figure4_mechanism_localization.py",
    "scripts/rebuild_prince_from_public_github.py": "scripts/rebuild_prince_from_public_github.py",
    "scripts/build_public_reproducibility_package.py": "submission_reproducibility_staging_v1/scripts/build_public_reproducibility_package.py",
    "results/GSE248014_locked_Mito3_summary_v1.csv": "data/03_results/external_ici/GSE248014/GSE248014_locked_Mito3_summary_v1.csv",
    "results/PRINCE_official_source_manifest.tsv": "data/03_results/official_prince_public_reproduction_v1/official_source_manifest.tsv",
    "results/PRINCE_historical_seven_module_screen.tsv": "data/03_results/official_prince_public_reproduction_v1/historical_seven_module_screen.tsv",
    "results/PRINCE_corrected_seven_module_screen.tsv": "data/03_results/official_prince_public_reproduction_v1/corrected_seven_module_screen.tsv",
    "results/PRINCE_official_public_Mito3_primary_statistics.tsv": "data/03_results/official_prince_public_reproduction_v1/official_public_Mito3_primary_statistics.tsv",
    "results/PRINCE_official_vs_current_metric_comparison.tsv": "data/03_results/official_prince_public_reproduction_v1/official_vs_current_metric_comparison.tsv",
    "results/GSE240078_module_coverage.tsv": "data/03_results/feasibility/gse240078/module_coverage.tsv",
    "results/GSE240078_paired_tumor_stroma_tests.tsv": "data/03_results/feasibility/targeted_cptac_spatial/spatial_paired_tumor_stroma_tests.tsv",
    "results/GSE240078_tumor_stroma_myeloid_association.tsv": "data/03_results/feasibility/targeted_cptac_spatial/spatial_tumor_module_vs_stroma_myeloid.tsv",
    "results/TCGA_PAAD_RPPA_targeted_spearman.tsv": "data/03_results/feasibility/tcga_paad_rppa_pakt/targeted_rppa_pairwise_spearman.tsv",
    "results/TCGA_PAAD_RPPA_module_gene_coverage.tsv": "data/03_results/feasibility/tcga_paad_rppa_pakt/module_gene_coverage.tsv",
    "results/PRISM_complete_state_selectivity.tsv": "data/03_results/feasibility/prism_state_selective_screen/all_compounds_state_selectivity.tsv",
    "results/PRISM_state_selectivity_report.md": "data/03_results/feasibility/prism_state_selective_screen/PRISM_state_selective_screen.md",
    "tables/Supplementary_Table_S1_cohort_and_dataset_inventory.tsv": "manuscript/tables/Supplementary_Table_S1_cohort_and_dataset_inventory.tsv",
    "tables/Supplementary_Table_S4_GeoMx_compartmentalization.tsv": "manuscript/tables/Supplementary_Table_S4_GeoMx_compartmentalization.tsv",
    "tables/Supplementary_Table_S5_TCGA_RPPA_targeted_associations.tsv": "manuscript/tables/Supplementary_Table_S5_TCGA_RPPA_targeted_associations.tsv",
    "tables/Supplementary_Table_S6_PRISM_complete_state_selectivity_screen.tsv": "manuscript/tables/Supplementary_Table_S6_PRISM_complete_state_selectivity_screen.tsv",
    "figures/Figure3_external_ICI_audit_v2.pdf": "manuscript/figures/Figure3_external_ICI_audit_v2.pdf",
    "figures/Figure4_immune_context_localization.pdf": "manuscript/figures/Figure4_immune_context_localization.pdf",
    "figures/Figure5_protein_drug_prioritization.pdf": "manuscript/figures/Figure5_protein_drug_prioritization.pdf",
    "figures/Supplementary_FigureS1_Mito3_frozen_sensitivity.pdf": "manuscript/figures/Supplementary_FigureS1_Mito3_frozen_sensitivity.pdf",
    "figures/Supplementary_FigureS4_GeoMx_coverage_and_roi_limits.pdf": "manuscript/figures/Supplementary_FigureS4_GeoMx_coverage_and_roi_limits.pdf",
    "figures/Supplementary_FigureS9_immune_context_localization.pdf": "manuscript/figures/Supplementary_FigureS9_immune_context_localization.pdf",
    "figures/Supplementary_FigureS10_protein_drug_prioritization.pdf": "manuscript/figures/Supplementary_FigureS10_protein_drug_prioritization.pdf",
}

FORBIDDEN_TERMS = ("clinical_response", "patient_scores", "participant_level", "deidentified_id")


README = """# Public reproducibility package (v8 candidate)

## Scope

This local release candidate contains public-resource provenance,
non-participant-level derived results, publication figures, and historical
scripts used for selected public-data analyses in the pancreatic cancer Mito3
study. It is an evidence and result-release package, not a complete executable
workflow archive or a public repository.

PRINCE deidentified limited clinical data and processed RNA-seq files are
publicly released in the official ParkerICI `prince-trial-data` repository at
commit `c01cf276b1cef27e61f7349bccfac37c0c1d6ab7` (Padrón *et al.*, 2022;
DOI `10.1038/s41591-022-01829-9`). This package does not redistribute PRINCE
source files because that repository has no explicit repository-level data
license. It contains a public-source reconstruction script, a fixed-source
manifest, and non-participant-level aggregate outputs only. The provenance audit distinguishes the candidate module
dictionary, PRINCE seven-module selection, final Mito3 definition, and the
scoring rule frozen before the GSE248014 response-association transportability
audit; it does not represent prospective preregistration.

According to institutional policy, analysis of these public, deidentified
secondary datasets does not constitute human-subjects research requiring ethics
review. No new ethics approval or exemption is claimed.

## Included materials

- `data/`: accession-level public-resource manifest.
- `documentation/`: source verification, evidence mapping, and provenance audit.
- `scripts/`: historical public-data analysis and figure-generation scripts.
- `results/`, `tables/`, and `figures/`: non-participant-level result snapshots
  and figures supporting the manuscript.
- `MANIFEST.tsv`: release-relative path, original project-relative path, byte
  count, and SHA-256 digest for each included artifact.

## Reuse and rerun guidance

The public accessions and URLs are listed in `data/PUBLIC_DATASET_MANIFEST.tsv`.
Scripts retain their original project-relative input layout and should be
treated as provenance-preserving analysis code. Before executing them in a new
environment, obtain source data from the cited public repositories, set a local
project root, and verify file versions against the accompanying source
verification report. For PRINCE, use the fixed official commit above and the
hashes in `results/PRINCE_official_source_manifest.tsv`. The package does not
distribute third-party raw data, software environments, or external-tool archives.

## Interpretation boundary

Results here support exploratory molecular-context and pharmacologic-priority
analyses. They do not establish a causal mechanism, drug efficacy, a
deployable clinical response biomarker, calibration, clinical utility, or ICI
treatment-effect modification. The GSE248014 analysis is an
response-association transportability audit in a cohort with three responders;
its performance interval is not stably interpretable and its direction did not
reproduce the PRINCE association.

## Release exclusions

- PRINCE clinical and expression source files, participant-level outputs, and response tables.
- Third-party raw data and any access-controlled materials.
- Local caches, virtual environments, downloaded third-party software, and raw
  data whose redistribution terms are not verified here.
- Reproducible environment lock files and the archival DOI, which remain
  required before public release.

## Integrity check

Run `python verify_package.py` from this directory. It verifies that every
released file still matches `MANIFEST.tsv` and that no release path contains a
restricted-cohort marker.
"""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_safe_path(relative_path: str, source_relative_path: str) -> None:
    combined = f"{relative_path} {source_relative_path}".lower()
    if any(term in combined for term in FORBIDDEN_TERMS):
        raise ValueError(f"Refusing restricted-cohort material: {source_relative_path}")


def write_verifier() -> None:
    verifier = '''from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FORBIDDEN_TERMS = ("clinical_response", "patient_scores", "participant_level", "deidentified_id")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


with (ROOT / "MANIFEST.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\\t"))

errors = []
for row in rows:
    relative_path = row["release_relative_path"]
    if any(term in relative_path.lower() for term in FORBIDDEN_TERMS):
        errors.append(f"restricted marker in path: {relative_path}")
        continue
    path = ROOT / relative_path
    if not path.is_file():
        errors.append(f"missing: {relative_path}")
    elif path.stat().st_size != int(row["bytes"]):
        errors.append(f"size mismatch: {relative_path}")
    elif sha256(path) != row["sha256"]:
        errors.append(f"SHA-256 mismatch: {relative_path}")

if errors:
    raise SystemExit("Package verification failed:\\n- " + "\\n- ".join(errors))
print(f"Verified {len(rows)} released files; no restricted path markers found.")
'''
    (PACKAGE / "verify_package.py").write_text(verifier, encoding="utf-8")


def main() -> None:
    if PACKAGE.exists():
        raise FileExistsError(f"Refusing to overwrite existing package: {PACKAGE}")
    PACKAGE.mkdir()
    (PACKAGE / "README.md").write_text(README, encoding="utf-8")

    manifest_rows = []
    for release_relative_path, source_relative_path in INCLUDED_FILES.items():
        if source_relative_path == "_generated_readme_":
            continue
        verify_safe_path(release_relative_path, source_relative_path)
        source = ROOT / source_relative_path
        if not source.is_file():
            raise FileNotFoundError(source)
        destination = PACKAGE / release_relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        manifest_rows.append((release_relative_path, source_relative_path, destination.stat().st_size, sha256(destination)))

    manifest = PACKAGE / "MANIFEST.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        handle.write("release_relative_path\toriginal_project_relative_path\tbytes\tsha256\n")
        for row in manifest_rows:
            handle.write("\t".join(map(str, row)) + "\n")
    write_verifier()
    print(f"Built {PACKAGE} with {len(manifest_rows)} released files.")


if __name__ == "__main__":
    main()
