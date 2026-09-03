# PRINCE official-public-source reproduction audit (v1)

Date: 2026-09-02

## Scope

This audit independently rebuilt the PRINCE discovery analysis from the
official ParkerICI public release rather than the previously used
SecAct-distributed expression matrix. It evaluates reproducibility of the
existing Mito3 result; it does not change the post-discovery status of the
PRINCE screen or convert the association into biomarker validation.

## Fixed source

- Study: Padrón *et al.*, *Nature Medicine* 2022, DOI
  `10.1038/s41591-022-01829-9`.
- Official release: `https://github.com/ParkerICI/prince-trial-data` at commit
  `c01cf276b1cef27e61f7349bccfac37c0c1d6ab7`.
- Downloaded archive SHA-256:
  `3f7e268858214ef126eeef2b0004ea80ed66864eecd2419471714068fbf2002a`.
- The original article's Data availability statement identifies this release
  as containing a deidentified limited clinical dataset and processed RNA-seq
  files. It states that the complete clinical dataset is commercially sensitive
  and that confidential or identifiable patient-level raw data cannot be
  shared.

The official clinical file and RNA metadata file are byte-identical to the
copies previously retained under `incoming_data/prince-trial-data_minimal/`.

## Independent reconstruction

The rebuild read 65 per-sample public RNA reports, transformed expression as
`log2(TPM + 1)`, standardized each gene over the 65 baseline RNA records, and
calculated Mito3 as the mean standardized expression of `BAX`, `BAK1`, and
`BID`. It linked RNA metadata to the public clinical table, retained 45
nivolumab-exposed RNA records, then retained 38 with a binary best-overall
response (18 responders and 20 nonresponders).

## Reproduction result

All 17 pre-specified numerical checks match the existing primary result to
floating-point precision (maximum absolute difference `3.33e-16`):

- responder-minus-nonresponder median difference: `-0.832975`;
- two-sided Wilcoxon P value: `0.0316509`;
- lower-Mito3 rank AUC: `0.705556`; rank-biserial correlation: `0.411111`;
- primary 10,000-resample bootstrap percentile interval for the median
  difference: `-1.423514` to `-0.218976`;
- bootstrap percentile interval for rank AUC: `0.519444` to `0.872222`;
- lower-in-responder bootstrap direction probability: `0.9953`;
- 20,000-permutation arm-stratified P value: `0.0174991`;
- arm-adjusted odds ratio per standard deviation: `0.442741` (95% CI
  `0.206982` to `0.947037`; P=`0.0357083`).

The historical PPIF-inclusive mitochondrial score reproduces raw P=`0.0186007`
and BH-FDR=`0.0651023`; the corrected final Mito3 (`BAX`/`BAK1`/`BID`) has raw
P=`0.0316509` and seven-module BH-FDR=`0.1107781`. These are distinct score
definitions and must remain separately labelled.

## Files

- Rebuild script: `scripts/rebuild_prince_from_public_github.py`.
- Result directory:
  `data/03_results/official_prince_public_reproduction_v1/`.
- The result directory includes the fixed-source manifest, historical and
  corrected seven-module screens, primary statistics, and exact comparisons
  against the prior local result.

## Release implication

The manuscript can accurately identify the PRINCE limited clinical and
processed RNA-seq source as publicly released at the official repository and
can cite its fixed commit. This finding does **not** authorize rehosting the
source files: the GitHub repository has no explicit repository-level data
license. The project should publish code, source manifests and derived
aggregate outputs, while linking readers to the original repository for PRINCE
source data.
