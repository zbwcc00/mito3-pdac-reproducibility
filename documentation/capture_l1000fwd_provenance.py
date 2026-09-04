"""Capture auditable L1000FWD API provenance for the three locked queries."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import requests


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "documentation" / "l1000fwd_provenance_2026-09-04"
BASE_URL = "https://maayanlab.cloud/l1000fwd/"
QUERIES = {
    "high_module_state": {
        "up_genes": ["BAX", "BAK1", "BID", "TLR2", "TLR4", "TLR7", "TLR8", "MYD88", "TNFRSF1A", "TNFRSF1B", "TNF", "IL1B", "SPP1"],
        "down_genes": [],
    },
    "high_inflammatory_low_mito": {
        "up_genes": ["TLR2", "TLR4", "TLR7", "TLR8", "MYD88", "TNFRSF1A", "TNFRSF1B", "TNF", "IL1B", "SPP1"],
        "down_genes": ["BAX", "BAK1", "BID"],
    },
    "high_mito_low_inflammatory": {
        "up_genes": ["BAX", "BAK1", "BID"],
        "down_genes": ["TLR2", "TLR4", "TLR7", "TLR8", "MYD88", "TNFRSF1A", "TNFRSF1B", "TNF", "IL1B", "SPP1"],
    },
}


def utc_now():
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    records = []
    for name, request_body in QUERIES.items():
        search_url = BASE_URL + "sig_search"
        search_started = utc_now()
        search = session.post(search_url, json=request_body, timeout=90)
        search_body = search.content
        search_path = OUTPUT / f"{name}.sig_search.response.json"
        search_path.write_bytes(search_body)
        search_json = search.json()
        result_id = search_json["result_id"]

        result_url = BASE_URL + "result/topn/" + result_id
        result_started = utc_now()
        result = session.get(result_url, timeout=90)
        result_body = result.content
        result_path = OUTPUT / f"{name}.result_topn.response.json"
        result_path.write_bytes(result_body)
        result_json = result.json()

        records.append(
            {
                "query_name": name,
                "search": {
                    "retrieval_started_utc": search_started,
                    "retrieval_completed_utc": utc_now(),
                    "method": "POST",
                    "endpoint": search_url,
                    "request_json": request_body,
                    "http_status": search.status_code,
                    "result_id": result_id,
                    "response_file": search_path.name,
                    "response_sha256": sha256_bytes(search_body),
                },
                "result_topn": {
                    "retrieval_started_utc": result_started,
                    "retrieval_completed_utc": utc_now(),
                    "method": "GET",
                    "endpoint": result_url,
                    "http_status": result.status_code,
                    "response_file": result_path.name,
                    "response_sha256": sha256_bytes(result_body),
                    "top_level_keys": sorted(result_json),
                },
            }
        )

    input_files = [
        ROOT.parent / "incoming_data" / "LINCS_L1000FWD" / "CD_signature_metadata.csv",
        ROOT.parent / "incoming_data" / "LINCS_L1000FWD" / "CD_signatures_binary_42809.gmt",
    ]
    inputs = [
        {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in input_files
    ]
    try:
        script_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        script_commit = None

    log = {
        "provenance_schema": "l1000fwd-api-retrieval-v1",
        "captured_at_utc": utc_now(),
        "service": "L1000FWD",
        "base_url": BASE_URL,
        "source_release": None,
        "source_release_note": "The public API does not expose a release identifier in the captured response; endpoint and response bytes are therefore pinned by URL, timestamp, result_id, and SHA-256.",
        "script_commit_at_capture": script_commit,
        "queries": records,
        "input_files": inputs,
        "interpretation_limitation": "This is an exploratory top-50 signature-ranking screen. It does not establish clinical efficacy, and cell line, dose, time, and batch heterogeneity remain unresolved.",
    }
    (OUTPUT / "retrieval_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "queries": len(records), "script_commit": script_commit}, indent=2))


if __name__ == "__main__":
    main()
