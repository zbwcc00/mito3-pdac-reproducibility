from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    paths = sorted(path for path in ROOT.rglob("*") if path.is_file() and path != MANIFEST)
    rows = [(path.relative_to(ROOT).as_posix(), "staging-package", path.stat().st_size, sha256(path)) for path in paths]
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        handle.write("release_relative_path\toriginal_project_relative_path\tbytes\tsha256\n")
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Wrote {MANIFEST} with {len(rows)} files.")


if __name__ == "__main__":
    main()
