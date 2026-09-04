# Mito3 reproducibility release (v1.0.2)

## Scope

This release contains public-resource provenance,
non-participant-level derived results, publication figures, and historical
scripts used for selected public-data analyses in the pancreatic cancer Mito3
study. It is a versioned public code and result release; source data are
retrieved from the cited providers and are not bundled.

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
- `scripts/analyze_prince_mito3_survival.py`: post hoc OS/PFS sensitivity analysis and Kaplan–Meier export.
- `environment/`: locked Python and R package versions plus OS requirements.
- `results/`, `tables/`, and `figures/`: non-participant-level result snapshots
  and figures supporting the manuscript.
- `CITATION.cff` and `zenodo.json`: citation and archival-metadata templates.
- `documentation/third_party_license_and_checksum_audit.md`: redistribution
  decisions and source-verification boundaries.
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
distribute third-party raw data, virtual environments, or external-tool
archives. Run `scripts/check_source_urls.py` to refresh the date-stamped URL
reachability record.

## Interpretation boundary

Results here support exploratory molecular-context and pharmacologic-priority
analyses. They do not establish a causal mechanism, drug efficacy, a
deployable clinical response biomarker, calibration, clinical utility, or ICI
treatment-effect modification. The PRINCE OS/PFS outputs are post hoc
survival sensitivity analyses and do not establish prognostic or predictive
utility. The GSE248014 analysis is an
response-association transportability audit in a cohort with three responders;
its performance interval is not stably interpretable and its direction did not
reproduce the PRINCE association.

## Release exclusions

- PRINCE clinical and expression source files, participant-level outputs, and response tables.
- Third-party raw data and any access-controlled materials.
- Local caches, virtual environments, downloaded third-party software, and raw
  data whose redistribution terms are not verified here.
- Release `v1.0.1` remains archived at DOI `10.5281/zenodo.22282838`.
- Release `v1.0.2` adds the PRINCE OS/PFS survival sensitivity analysis. Its
  version-specific DOI will be recorded after the Zenodo archive is published;
  the concept DOI for all versions is `10.5281/zenodo.22281618`.

## Integrity check

Run `python verify_package.py` from this directory. It verifies that every
released file still matches `MANIFEST.tsv` and that no release path contains a
restricted-cohort marker.
