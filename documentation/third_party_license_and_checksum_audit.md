# Third-party licence and checksum audit

**Audit date:** 2026-09-03  
**Release candidate:** `public_reproducibility_repo_v1`

## Scope and method

Canonical accessions, versions, local evidence, and checksums are recorded in
`data/CITATION_AND_DATA_SOURCE_VERIFICATION_MANIFEST.tsv`. The URL checker in
`scripts/check_source_urls.py` records HTTP status and content type in
`documentation/source_url_status.tsv`. The audit does not infer a licence from
the absence of a licence file.

## Redistribution decisions

| Resource | Evidence | Decision |
|---|---|---|
| ParkerICI `prince-trial-data` | GitHub API returned no repository licence; fixed commit and archive SHA-256 are in `results/PRINCE_official_source_manifest.tsv`. | Do not rehost source CSV/TSV files; link to the official commit and release only aggregate results and reconstruction code. |
| DepMap 24Q4 Public | Figshare record `10.25452/figshare.plus.27993248.v1` declares CC BY 4.0; local MD5 values are in the citation manifest. | Do not bundle raw files; cite the DOI and retain local checksums. |
| PRISM Repurposing 19Q3/19Q4 primary screen | Figshare record `10.6084/m9.figshare.9393293.v4` declares CC BY 4.0; local MD5 values are in the citation manifest. | Do not bundle raw files; cite the DOI and retain local checksums. |
| GEO accessions | NCBI GEO public records and the accession-specific source files are cited; GEO terms remain applicable. | Release accession links and derived non-participant-level summaries only. |
| cBioPortal Datahub | GitHub API returned no repository licence. | Release processing code and aggregate outputs; do not assume redistribution rights for downloaded study files. |
| SecAct and UniPert | Source repositories and model provenance are cited; UniPert source revision remains partial. | Keep third-party source trees/checkpoints out of the release and retain the recorded model checksum. |

## Checksums

The release `MANIFEST.tsv` contains SHA-256 values for every file shipped here.
PRINCE fixed-source checksums are retained in
`results/PRINCE_official_source_manifest.tsv`; DepMap and PRISM local checksums
are retained in the citation manifest. The L1000FWD retrieval log and a stable
UniPert source archive are still required before claiming a fully closed source
audit.

## Interpretation

An HTTP 200 response confirms that a source URL was reachable on the audit date;
it does not grant redistribution permission. A missing repository licence is
treated conservatively as “link only”.
