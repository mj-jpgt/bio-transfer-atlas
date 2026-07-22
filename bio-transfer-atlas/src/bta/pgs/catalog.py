"""
PGS Catalog API helpers.
"""

import time
from typing import Any

import requests


API_BASE = "https://www.pgscatalog.org/rest"


def search_scores(trait: str, page_size: int = 100) -> list[dict[str, Any]]:
    results = []
    page = 1
    while True:
        r = requests.get(
            f"{API_BASE}/score/search",
            params={"trait_search": trait, "limit": page_size, "offset": (page - 1) * page_size},
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
        batch = data.get("results", [])
        results.extend(batch)
        if not data.get("next"):
            break
        page += 1
        time.sleep(0.25)
    return results


def get_score_metadata(pgs_id: str) -> dict[str, Any]:
    r = requests.get(f"{API_BASE}/score/{pgs_id}", timeout=30)
    r.raise_for_status()
    return r.json()


def harmonized_score_url(pgs_id: str, build: str = "GRCh37") -> str:
    return (
        f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/"
        f"{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_{build}.txt.gz"
    )
