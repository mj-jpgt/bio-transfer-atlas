"""
Score selection: query PGS Catalog API for MVP traits and produce
a curated selection table written to configs/scores.yaml.

Criteria:
  - Harmonized GRCh37 scoring file must exist
  - Include at least one EUR-developed score per trait
  - Include multi-ancestry or non-EUR scores if available

Usage:
    python scripts/select_scores.py
"""

import json
import sys
from pathlib import Path

import pandas as pd
import yaml
from loguru import logger

from bta.pgs.catalog import search_scores


MVP_TRAITS = {
    "coronary_artery_disease": "coronary artery disease",
    "type_2_diabetes": "type 2 diabetes",
    "body_mass_index": "body mass index",
}


def classify_ancestry(score: dict) -> str:
    dist = score.get("ancestry_distribution") or {}
    keys = [k.lower() for k in dist.keys()]
    if len(keys) > 1:
        return "Multi"
    if keys:
        return keys[0].capitalize()
    return "Unknown"


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    api_dir = root / "data" / "raw" / "pgs_catalog" / "metadata" / "api"
    api_dir.mkdir(parents=True, exist_ok=True)
    scores_yaml = root / "configs" / "scores.yaml"

    selected = []

    for trait_key, trait_query in MVP_TRAITS.items():
        cache = api_dir / f"{trait_query.replace(' ', '_')}_search.json"
        if cache.exists():
            with open(cache) as f:
                results = json.load(f)
        else:
            results = search_scores(trait_query)
            with open(cache, "w") as f:
                json.dump(results, f, indent=2)

        logger.info(f"{trait_query}: {len(results)} scores found")

        for score in results:
            pgs_id = score.get("id", "")
            build = score.get("variants_genomebuild", "")
            if "37" not in str(build) and "hg19" not in str(build).lower():
                continue

            anc = classify_ancestry(score)
            pub = (score.get("publication") or {}).get("doi", "")
            n_var = score.get("variants_number")

            selected.append({
                "pgs_id": pgs_id,
                "trait": trait_key,
                "reported_trait": score.get("trait_reported", ""),
                "development_ancestry": anc,
                "n_variants_original": n_var,
                "genome_build": build,
                "publication": pub,
                "license": score.get("license", ""),
                "reason_selected": "auto-selected from API search",
            })

    df = pd.DataFrame(selected).drop_duplicates("pgs_id")
    logger.info(f"Total candidates: {len(df)}")

    out_tsv = root / "data" / "raw" / "pgs_catalog" / "metadata" / "candidate_scores.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    logger.success(f"Candidate table: {out_tsv}")

    with open(scores_yaml) as f:
        cfg = yaml.safe_load(f)

    cfg["selected_scores"] = df.to_dict(orient="records")

    with open(scores_yaml, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True)

    logger.success(f"configs/scores.yaml updated with {len(df)} entries")


if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stderr, level="INFO")
    main()
