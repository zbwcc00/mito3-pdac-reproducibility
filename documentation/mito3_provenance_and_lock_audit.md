# Mito3 provenance and analysis-lock audit (v1)

Audit date: 2026-08-29

## Audit question

This audit asks four separate questions that must not be conflated:

1. Were `BAX`, `BAK1`, and `BID` identified as a mechanistic mitochondrial core before access to PRINCE response labels?
2. Was the three-gene `Mito3` score designated as the sole primary PRINCE response score before examining PRINCE outcomes?
3. Was the final three-gene scoring rule and direction fixed before the independent GSE248014 audit?
4. Do the retained records constitute prospective preregistration?

## Executive finding

The retained local record supports a **two-stage provenance**, not a prospectively preregistered primary biomarker analysis.

- The `BAX`/`BAK1`/`BID` mitochondrial-effector core was mechanistically specified before the PRINCE clinical files entered the project.
- The exact focal name `Mito3` and its elevation to the manuscript's primary response-associated score occurred after the initial PRINCE seven-module response screen. PRINCE must therefore remain a discovery/internal-association cohort.
- The corrected three-gene rule, observed direction, and intended analysis were documented before the GSE248014 expression file entered the project. GSE248014 can therefore be described as a locally documented frozen transportability audit.
- The record is not prospective preregistration. A Git repository and Zenodo
  release were created after the analysis and therefore do not provide a
  pre-analysis public timestamp, OSF registration, or other prospective record.
  NTFS creation and modification times remain supportive provenance only.

## Reconstructed timeline

| Local time | Evidence | What it establishes | What it does not establish |
|---|---|---|---|
| 2026-08-24 22:48 | `第九篇大论文_线粒体外周氧化应激外科生信课题设计_v1.md` | The proposed mechanism already included `BAX/BAK1/BID` as the mitochondrial damage layer. | The primary clinical aim was postoperative recurrence/ERFS, not PRINCE ICI response. |
| 2026-08-24 23:11 | `PDAC_公共数据库数据清单与v2研究方案_2026-08-24.md` | The project again specified a `BAX/BAK1/BID` mitochondrial oxidative layer and continuous module scoring before PRINCE clinical-data access. | It did not designate a three-gene score as the sole primary ICI-response hypothesis. |
| 2026-08-24 23:23 | `PDAC_mitoxyperilysis_gene_mechanism_evidence.tsv` | `BAX`, `BAK1`, and `BID` were each recorded as core mitochondrial-oxidative features; `PPIF` was separately recorded as a negative control. | This evidence table was not an immutable registration and did not define a PRINCE response analysis. |
| 2026-08-25 23:58 | `pdac_mechanism_gene_set_audit.tsv` and its generating script | The audit explicitly assigned three genes to the mitochondrial `core_state_score` and `PPIF` to `negative_control`. | It did not yet name `Mito3` or identify it as the sole primary clinical score. |
| 2026-08-26 20:43 | PRINCE clinical files were created under `incoming_data/prince-trial-data_minimal/`. | The retained local record places clinical response-label availability after the mechanistic gene audit. | A local file time cannot prove that no investigator had prior conceptual knowledge of the trial results. |
| 2026-08-26 20:48 | First PRINCE `response_tests.tsv` | Seven module-level response associations were examined. The implementation grouped genes by the broad `module` column, causing `PPIF` to be included in the initial mitochondrial score despite its negative-control role. | This was not a valid final three-gene Mito3 analysis and was not a preregistered single-primary test. |
| 2026-08-27 22:04-22:42 | PRINCE robustness and locked-model scripts | The mitochondrial module was identified as the strongest response-associated candidate and increasingly treated as the focal model. | These files were created after the initial PRINCE results and therefore cannot establish pre-PRINCE primary-score designation. |
| 2026-08-28 16:06 | `locked_treatment_stratification_specification.md` | The name `Mito3`, genes `BAX/BAK1/BID`, lower-in-responders direction, continuous scoring, and external-validation rule were explicitly fixed. The file itself states that the score was frozen after existing PRINCE and external GEO analyses. | It does not retroactively convert PRINCE into a confirmatory cohort. |
| 2026-08-28 18:52 | `mito3_definition_correction_v1.md` | The implementation error was documented, `PPIF` was removed, and `Mito3 = mean(z(BAX), z(BAK1), z(BID))` was made operational. The correction retained the direction and nominal PRINCE association. | Because this occurred after PRINCE inspection, the corrected PRINCE result remains discovery evidence. |
| 2026-08-28 20:19 | `GSE248014_table_RNAseq.csv.gz` entered the project. | The raw external expression file arrived after the explicit three-gene lock and correction record. | This is a local chronology, not a public preregistration. |
| 2026-08-28 20:29 | `GSE248014_locked_Mito3_summary_v1.csv` | The frozen external audit produced an opposite, imprecise point estimate. | It does not prove biological reversal, but it fails to reproduce the PRINCE direction. |

## Adjudication

### 1. Was the three-gene biological core specified before PRINCE labels?

**Supported, with moderate confidence.** Multiple files created on 2026-08-24 and 2026-08-25 identify `BAX`, `BAK1`, and `BID` as the mitochondrial core, while the retained PRINCE clinical files were created on 2026-08-26.

### 2. Was Mito3 the prespecified sole primary PRINCE response score?

**Not supported.** The initial project aim concerned postoperative recurrence/ERFS, the first PRINCE analysis tested seven modules, and the first exact `Mito3` records were created after the initial PRINCE signal was known. Calling PRINCE a confirmatory validation cohort would be misleading.

### 3. Was the final score frozen before GSE248014?

**Supported as a local retrospective lock.** The explicit score definition and direction were documented at 16:06, the correction record was created at 18:52, and the GSE248014 expression file entered the project at 20:19 on 2026-08-28.

### 4. Is this prospective preregistration?

**No.** The lock was created after discovery analyses, was stored locally, and lacks immutable third-party timestamping. It can support transparent sequencing of discovery and external audit, not prospective registration.

## Required manuscript language

Recommended wording:

> The BAX/BAK1/BID mitochondrial-effector module was mechanistically specified before local access to the PRINCE clinical files. After exploratory comparison of seven prespecified modules in PRINCE, the three-gene Mito3 score was selected as the focal response-associated candidate. Its gene membership, scoring rule, and lower-in-responders direction were then frozen before the GSE248014 transportability audit. This was a retrospective analysis lock, not prospective preregistration.

Acceptable short label:

> retrospectively defined and frozen before external transportability audit

Avoid:

- `prospectively prespecified Mito3 biomarker`
- `preregistered primary score`
- `externally validated analysis-locked biomarker`
- any wording implying that the sole Mito3 hypothesis preceded the first PRINCE response screen

## Consequences for statistical interpretation

1. PRINCE should be labelled discovery or internal association evidence.
2. The seven-module BH FDR is decision-relevant multiplicity context, not a disposable secondary statistic.
3. PRINCE leave-one-out, bootstrap, component-deletion, and arm-adjusted analyses are internal robustness checks, not validation.
4. GSE248014 is the cleanest available frozen external audit because the score and direction predate local acquisition of its expression file.
5. The opposite GSE248014 estimate prevents a stable or transportable biomarker claim but, with three responders and a different regimen, does not prove biological effect reversal.

## Provenance limitations

- The GitHub repository and Zenodo release were created after this analysis;
  they document the final release but do not retroactively establish
  prospective preregistration.
- No OSF/AsPredicted/public protocol timestamp was identified.
- Windows creation and modification times can be altered by copying or manual changes.
- Current SHA-256 hashes preserve the state observed on 2026-08-29 but do not retroactively make earlier files immutable.
- The available record supports transparent retrospective chronology, not proof of investigator blinding.

## Audit conclusion

The strongest defensible provenance claim is narrower than the current title may imply: **the biological three-gene core predates PRINCE label access, but its selection as the focal clinical score is post-discovery; the final score was then frozen before GSE248014.** This distinction should be carried consistently through the title, Abstract, Methods, Figure 1, and Discussion during the next manuscript revision.
