# Archival DOI handoff

The release is prepared locally, but no GitHub or Zenodo credentials are
available in this environment, so an archival DOI has not been minted.

1. Create the public GitHub repository and push branch `main` from this
   directory. The first release commit is recorded in the local Git history.
2. Create a GitHub release tagged `v1.0.0` from the release commit.
3. Enable the repository in Zenodo, import the GitHub release, review the
   metadata in `zenodo.json`, and publish the Zenodo record.
4. Replace the DOI placeholder in the manuscript Data and code availability
   statement and `CITATION.cff` with the DOI returned by Zenodo.
5. Record the DOI, repository URL, release tag, and final commit in the
   submission audit, then rerun `verify_package.py`.

Do not upload PRINCE source CSV/TSV files, participant-level outputs, raw GEO
files, virtual environments, or third-party software archives. The source
licence and checksum boundaries in
`documentation/third_party_license_and_checksum_audit.md` apply to the upload.
