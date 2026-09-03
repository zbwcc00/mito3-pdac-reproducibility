"""Check canonical source URLs recorded in the release manifest."""

from __future__ import annotations

import csv
import re
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "CITATION_AND_DATA_SOURCE_VERIFICATION_MANIFEST.tsv"
OUTPUT = ROOT / "documentation" / "source_url_status.tsv"
URL_PATTERN = re.compile(r"https?://[^\s;]+")


def check(url: str) -> tuple[str, str, str]:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Mito3-reproducibility-audit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return str(response.status), response.headers.get_content_type(), ""
    except urllib.error.HTTPError as error:
        if error.code not in {403, 405}:
            return str(error.code), "", str(error.reason)
        request = urllib.request.Request(url, method="GET", headers={"User-Agent": "Mito3-reproducibility-audit/1.0", "Range": "bytes=0-1023"})
        try:
            with urllib.request.urlopen(request, timeout=20) as response:
                return str(response.status), response.headers.get_content_type(), "GET fallback"
        except Exception as fallback_error:  # pragma: no cover - network dependent
            return "ERROR", "", str(fallback_error)
    except Exception as error:  # pragma: no cover - network dependent
        return "ERROR", "", str(error)


def main() -> None:
    urls: set[str] = set()
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            urls.update(URL_PATTERN.findall(row.get("canonical_source", "")))
    checked_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    rows = []
    for url in sorted(urls):
        status, content_type, note = check(url)
        rows.append({"url": url, "status": status, "content_type": content_type, "checked_at_utc": checked_at, "note": note})
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Checked {len(rows)} canonical URLs; wrote {OUTPUT}")


if __name__ == "__main__":
    main()
