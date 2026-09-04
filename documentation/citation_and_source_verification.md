# Citation and data-source verification report (v1)

**Audit date:** 2026-09-04 (post-release status update; source checks below were performed on 2026-08-29)
**Scope frozen before verification:** `citation_and_data_source_claim_inventory_v1.tsv` (13 literature claims, 13 data resources, and 6 software/portal resources).

## Overall conclusion

All 13 pre-existing references were verified against Crossref, and their
titles, journals, bibliographic fields, and available PubMed/Europe PMC
abstracts support the bounded narrative roles assigned in the manuscript.
No DOI/author/title mismatch was found. The audit identified a completeness
gap—not a false citation: the original draft cited no source studies for the
GEO cohorts and did not cite the software or portal resources used in Methods.
References 14–28 now close that gap.

## Claim-level outcome

| Scope | Outcome | Audit conclusion |
|---|---:|---|
| Existing literature claims L01–L13 | 13/13 pass | Each DOI resolved to the cited title and supported the stated context; none was used to claim a causal Mito3 mechanism or clinical efficacy. |
| GEO/CPTAC/TCGA/DepMap/PRISM source claims D02–D12 | 11/11 pass | GEO accession pages and linked original papers support cohort identity; cBioPortal study IDs and exact DepMap/PRISM file hashes were confirmed. |
| Software/portal claims S01–S04, S06 | 5/5 pass | Versions and primary citations were confirmed from local package metadata, the runtime, or official records. |
| Provenance items D13, S05 | 2/2 pass | D13 has a dated API capture with request/response hashes; S05 has an immutable UniPert commit, source archive checksum, and pinned model URL/local checksum. D01 is closed by the fixed ParkerICI commit, archive checksum, and public release link. |

## Corrected timeline wording

GSE179351 is now described as applying “the same post-discovery three-gene
score and response grouping.” It is never called a locked or frozen validation
cohort. “Retrospectively defined and frozen” is retained only for the final
score definition before the GSE248014 transportability audit.

## Material source-verification findings

1. **GSE179351:** the GEO record and Parikh *et al.* confirm radiation plus
   dual checkpoint blockade in microsatellite-stable colorectal and pancreatic
   adenocarcinoma. The paper-level cohort is larger than the six mapped PDAC
   baseline samples used here; the manuscript correctly limits this to a
   directional sensitivity comparison.
2. **GSE248014:** the GEO record and Baretti *et al.* confirm the 14-day
   entinostat lead-in and entinostat–nivolumab regimen. This validates the
   reported treatment context, not external reproducibility of the Mito3
   association.
3. **GSE240078:** GEO lists 242 records, whereas the quantile-normalized
   matrix reconciles to 223 analyzable ROIs from 40 paired patients. The
   manuscript should retain these two different counts and state that it uses
   the reconciled ROI matrix. No geometric ROI coordinates or cell segmentation
   were found in the release.
4. **CPTAC/TCGA:** cBioPortal metadata confirms `paad_cptac_2021`,
   `paad_tcga_pan_can_atlas_2018`, and legacy `paad_tcga`. The legacy RPPA
   release is appropriately a processing sensitivity analysis, not an
   independent cohort.
5. **DepMap/PRISM:** local MD5 hashes exactly match the official DepMap 24Q4
   `CRISPRGeneDependency.csv` and PRISM primary-screen files. The
   reported data type and directionality are therefore traceable.

## Remaining provenance limitations

These are reproducibility requirements, not analytic defects:

- Link the PRINCE analysis to the official ParkerICI `prince-trial-data`
  repository at commit `c01cf276b1cef27e61f7349bccfac37c0c1d6ab7` and retain
  its archive and public-file checksums. Do not rehost PRINCE source files:
  the official repository has no explicit repository-level data license.
- The dated L1000FWD API capture is retained in
  `documentation/l1000fwd_provenance_2026-09-04/`; because the service exposes
  no release identifier, the snapshot is pinned by endpoint, UTC timestamps,
  request JSON, result IDs, raw response hashes, and the two input-file hashes.
- The UniPert source revision and archive checksum are retained in
  `documentation/unipert_provenance_2026-09-04.json`. The analysis remains
  labelled **UniPert-assisted target-space prioritization**, not UniPert-G2CP
  perturbation-phenotype validation.
- The previous archival DOI remains available at GitHub tag `v1.0.2` and DOI
  `10.5281/zenodo.22297472`. The `v1.0.3` tag is a separate immutable release
  candidate awaiting Zenodo archival; these new provenance records are
  post-tag evidence and must not be silently retrofitted into that tag.

## Updated manuscript boundary

Citation completion does not change the study inference: PRINCE remains a
discovery/internal association; GSE179351 remains directional only; GSE248014
remains the frozen external transportability audit whose opposite, imprecise
estimate leaves external transportability unsupported and unresolved.
