# Reproduction environment

The locked environments were captured on 2026-09-03 from the runtimes used to
generate the release artifacts.

- Operating system: Windows 11 x86_64 (UCRT); the analysis scripts use POSIX
  paths only after replacing the project-root variable.
- Python: 3.12.13; install the packages in `python-requirements.lock`.
- R: 4.4.3 (ucrt); Bioconductor packages were installed from the Bioconductor
  3.20 repository. The exact package set is in `R-session-lock.txt`.
- Source data are not bundled. Retrieve each accession from the URL in
  `data/PUBLIC_DATASET_MANIFEST.tsv` and verify the fixed-source records before
  running a script.

The lock files record package versions, not operating-system binaries. A new
reproduction should record its OS, interpreter versions, and any compiler or
BLAS differences in the run log.
