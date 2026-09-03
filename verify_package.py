from __future__ import annotations

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
    rows = list(csv.DictReader(handle, delimiter="\t"))

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
    raise SystemExit("Package verification failed:\n- " + "\n- ".join(errors))
print(f"Verified {len(rows)} released files; no restricted path markers found.")
