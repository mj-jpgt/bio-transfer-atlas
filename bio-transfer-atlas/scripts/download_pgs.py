"""
Download PGS Catalog metadata and harmonized scoring files for MVP traits.

Steps:
  1. Query PGS Catalog REST API for each MVP trait.
  2. Save raw JSON responses to data/raw/pgs_catalog/metadata/api/.
  3. Parse responses to build a candidate score table.
  4. Save candidate table to data/raw/pgs_catalog/metadata/candidate_scores.tsv.
  5. Download harmonized scoring files (GRCh37) for all candidates.

Usage:
    python scripts/download_pgs.py
    TRAITS="coronary artery disease,type 2 diabetes" python scripts/download_pgs.py
"""

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from loguru import logger
from tqdm import tqdm

API_BASE = "https://www.pgscatalog.org/rest"
FTP_SCORE_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"

MVP_TRAITS = [
    "coronary artery disease",
    "type 2 diabetes",
    "body mass index",
]

TRAITS = [
    t.strip()
    for t in os.environ.get("TRAITS", ",".join(MVP_TRAITS)).split(",")
    if t.strip()
]

BUILD = "GRCh37"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


TRAIT_KEYWORDS = {
    "coronary artery disease": ["coronary artery disease", "coronary heart disease", "cad", "myocardial infarction"],
    "type 2 diabetes": ["type 2 diabetes", "t2d", "type ii diabetes"],
    "body mass index": ["body mass index", "bmi"],
}


def matches_trait(reported: str, query: str) -> bool:
    r = reported.lower()
    return any(kw in r for kw in TRAIT_KEYWORDS.get(query, [query.lower()]))


def fetch_all_scores_filtered(trait_queries: list[str]) -> dict[str, list]:
    """Page through /rest/score/all and return scores matching any trait query."""
    matched: dict[str, list] = {t: [] for t in trait_queries}
    url = f"{API_BASE}/score/all"
    params: dict = {"limit": 100, "offset": 0}
    total = None
    seen = 0

    while True:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
        if total is None:
            total = data.get("count", 0)
            logger.info(f"Total scores in catalog: {total}")

        batch = data.get("results", [])
        if not batch:
            break
        seen += len(batch)

        for score in batch:
            reported = score.get("trait_reported", "") or ""
            for tq in trait_queries:
                if matches_trait(reported, tq):
                    matched[tq].append(score)

        logger.info(f"  Scanned {seen}/{total} scores ...")
        if not data.get("next"):
            break
        params["offset"] += params["limit"]
        time.sleep(0.15)

    return matched


def download_stream(url: str, dest: Path, desc: str = "") -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        logger.info(f"Already exists, skipping: {dest}")
        return
    logger.info(f"Downloading {desc or url} -> {dest}")
    r = requests.get(url, stream=True, timeout=60)
    r.raise_for_status()
    total = int(r.headers.get("content-length", 0))
    with open(dest, "wb") as f, tqdm(
        total=total, unit="B", unit_scale=True, desc=dest.name
    ) as bar:
        for chunk in r.iter_content(chunk_size=131072):
            f.write(chunk)
            bar.update(len(chunk))
    logger.success(f"Saved: {dest}  sha256={sha256_file(dest)[:12]}...")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    api_dir = root / "data" / "raw" / "pgs_catalog" / "metadata" / "api"
    scores_dir = root / "data" / "raw" / "pgs_catalog" / "scores"
    api_dir.mkdir(parents=True, exist_ok=True)

    cache_all = api_dir / "all_scores_filtered.json"

    if cache_all.exists():
        logger.info(f"Loading cached filtered results: {cache_all}")
        with open(cache_all) as f:
            matched = json.load(f)
    else:
        logger.info("Paging through PGS Catalog /score/all ...")
        matched = fetch_all_scores_filtered(TRAITS)
        with open(cache_all, "w") as f:
            json.dump(matched, f, indent=2)

    for trait, scores in matched.items():
        logger.success(f"  {trait}: {len(scores)} scores matched")

    all_candidates = []
    for trait, scores in matched.items():
        for score in scores:
            pgs_id = score.get("id", "")
            hm = (score.get("ftp_harmonized_scoring_files") or {})
            has_grch37 = BUILD in hm
            all_candidates.append({
                "pgs_id": pgs_id,
                "trait_query": trait,
                "reported_trait": score.get("trait_reported", ""),
                "development_ancestry": str(
                    (score.get("ancestry_distribution") or {}).get("dev", {})
                )[:120],
                "n_variants_original": score.get("variants_number"),
                "genome_build": score.get("variants_genomebuild", ""),
                "publication": (score.get("publication") or {}).get("doi", ""),
                "license": score.get("license", ""),
                "has_grch37_harmonized": has_grch37,
                "harmonized_url": hm.get(BUILD, {}).get("positions", "") if has_grch37 else "",
            })

    if not all_candidates:
        logger.error("No candidates found.")
        return

    df = pd.DataFrame(all_candidates)
    out_tsv = root / "data" / "raw" / "pgs_catalog" / "metadata" / "candidate_scores.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.success(f"Candidate scores table: {out_tsv}  ({len(df)} rows)")

    # Only download scores that have a GRCh37 harmonized file
    to_download = df[df["has_grch37_harmonized"] == True]
    logger.info(f"Downloading {len(to_download)} harmonized scoring files (build={BUILD}) ...")

    for _, row in tqdm(to_download.iterrows(), total=len(to_download), desc="Scoring files"):
        pgs_id = row["pgs_id"]
        url = row["harmonized_url"]
        if not url:
            continue
        fname = f"{pgs_id}_hmPOS_{BUILD}.txt.gz"
        dest = scores_dir / pgs_id / fname
        try:
            download_stream(url, dest, f"{pgs_id}")
        except requests.HTTPError as e:
            logger.warning(f"  Skipped {pgs_id}: {e}")
        time.sleep(0.05)

    logger.success("PGS Catalog downloads complete.")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    main()
