# Archival DOI handoff

The `v1.0.0` release is published on GitHub and archived in Zenodo.

- Repository: https://github.com/zbwcc00/mito3-pdac-reproducibility
- Release tag: `v1.0.0`
- Version DOI: `10.5281/zenodo.22281619`
- Concept DOI (all versions): `10.5281/zenodo.22281618`

1. Keep the published `v1.0.0` tag immutable.
2. Record any post-publication metadata correction as a new commit or release;
   do not rewrite the published tag.
3. Keep the DOI, repository URL, release tag, and final commit in the
   submission audit, then rerun `verify_package.py`.

Do not upload PRINCE source CSV/TSV files, participant-level outputs, raw GEO
files, virtual environments, or third-party software archives. The source
licence and checksum boundaries in
`documentation/third_party_license_and_checksum_audit.md` apply to the upload.
