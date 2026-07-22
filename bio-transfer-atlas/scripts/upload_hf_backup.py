"""
Upload archivable pipeline data to Hugging Face (cold backup).

Selective plan (default): only bulky artifacts safe to delete locally.
See logs/hf_archive_plan.json for the exact manifest.

Usage:
  set HF_TOKEN=...
  set HF_XET_HIGH_PERFORMANCE=1
  set HF_XET_CACHE=%TEMP%\\hf_xet_cache

  python scripts/upload_hf_backup.py --dry-run
  python scripts/upload_hf_backup.py --target features
  python scripts/upload_hf_backup.py --target vcf
  python scripts/upload_hf_backup.py --target interim
  python scripts/upload_hf_backup.py --verify --target features
  python scripts/upload_hf_backup.py --delete-local --target features
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
REPO_ID = "mj0jpgg/fairgen"

# Score pgens kept locally for chr8-21 pipeline + chr22 intervention
KEEP_SCORE_CHROMS = [str(c) for c in list(range(8, 22)) + [22]]
INTERIM_IGNORE = [f"**/chr{c}.score.*" for c in KEEP_SCORE_CHROMS]

ARCHIVE_TARGETS = {
    "features": {
        "local": DATA / "features",
        "hub": "backup/features",
        # .vcor LD caches are 285GB, OneDrive-unfriendly, regenerable — skip upload
        "ignore_patterns": ["**/*.vcor", "**/*.log"],
        "delete": "tree",
    },
    "vcf": {
        "local": DATA / "raw" / "1000g" / "vcf_grch38",
        "hub": "backup/raw/1000g/vcf_grch38",
        "ignore_patterns": None,
        "delete": "tree",
    },
    "interim": {
        "local": DATA / "interim",
        "hub": "backup/interim",
        "ignore_patterns": INTERIM_IGNORE,
        "delete": "interim_selective",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Backup archivable data to Hugging Face.")
    p.add_argument(
        "--target",
        choices=[*ARCHIVE_TARGETS.keys(), "all"],
        default="all",
    )
    p.add_argument("--repo-id", default=REPO_ID)
    p.add_argument("--verify", action="store_true")
    p.add_argument("--delete-local", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def iter_local_files(local: Path, ignore_patterns: list[str] | None) -> list[Path]:
    import fnmatch

    out: list[Path] = []
    for f in local.rglob("*"):
        if not f.is_file():
            continue
        rel = f.relative_to(local).as_posix()
        if ignore_patterns and any(fnmatch.fnmatch(rel, pat) for pat in ignore_patterns):
            continue
        out.append(f)
    return out


def target_stats(name: str) -> dict:
    cfg = ARCHIVE_TARGETS[name]
    local: Path = cfg["local"]
    if not local.exists():
        return {"target": name, "n_files": 0, "bytes": 0, "gb": 0.0}
    files = iter_local_files(local, cfg.get("ignore_patterns"))
    total = sum(f.stat().st_size for f in files)
    return {
        "target": name,
        "local": str(local.relative_to(ROOT)),
        "hub": cfg["hub"],
        "n_files": len(files),
        "bytes": total,
        "gb": round(total / (1024**3), 3),
        "ignore_patterns": cfg.get("ignore_patterns"),
    }


def write_plan(targets: list[str]) -> Path:
    plan = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "repo_id": REPO_ID,
        "keep_local": {
            "data/modeling": "master tables, models",
            "data/processed": "PGS weights, score matrices",
            "data/labels": "concordance labels",
            "data/annotations": "gene/pathway maps",
            "data/raw/panukbb,bbj,finngen": "GWAS chr8-21",
            "data/raw/ensembl,reference,reactome,gnomad": "reference",
            f"data/interim/1000g_grch38/chr{{{','.join(KEEP_SCORE_CHROMS)}}}.score.*": "active pgens",
        },
        "targets": [target_stats(t) for t in targets],
    }
    out = ROOT / "logs" / "hf_archive_plan.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(plan, indent=2), encoding="utf-8")
    return out


def get_api():
    if not os.environ.get("HF_TOKEN") and not os.environ.get("HUGGING_FACE_HUB_TOKEN"):
        print("ERROR: Set HF_TOKEN in the environment.", file=sys.stderr)
        sys.exit(1)
    from huggingface_hub import HfApi

    return HfApi()


def upload_target(api, name: str, repo_id: str, dry_run: bool) -> None:
    cfg = ARCHIVE_TARGETS[name]
    local: Path = cfg["local"]
    if not local.exists():
        print(f"SKIP {name}: {local} missing")
        return
    st = target_stats(name)
    print(f"\n=== Upload {name} ===")
    print(f"  {st['n_files']:,} files, {st['gb']:.2f} GB")
    print(f"  {st['local']} -> {st['hub']}/")
    if cfg.get("ignore_patterns"):
        print(f"  ignore patterns: {cfg['ignore_patterns']}")
    if dry_run:
        return
    # Upload features by subfolder to avoid huge single-commit failures
    if name == "features":
        for sub in sorted(p for p in local.iterdir() if p.is_dir()):
            print(f"\n  -- subfolder {sub.name} --")
            api.upload_folder(
                repo_id=repo_id,
                folder_path=str(sub),
                path_in_repo=f"{cfg['hub']}/{sub.name}",
                repo_type="dataset",
                commit_message=f"Backup features/{sub.name}",
                ignore_patterns=cfg.get("ignore_patterns"),
            )
        print(f"  DONE: {cfg['hub']}")
        return
    kwargs = dict(
        repo_id=repo_id,
        folder_path=str(local),
        path_in_repo=cfg["hub"],
        repo_type="dataset",
        commit_message=f"Backup {name}",
    )
    if cfg.get("ignore_patterns"):
        kwargs["ignore_patterns"] = cfg["ignore_patterns"]
    api.upload_folder(**kwargs)
    print(f"  DONE: {cfg['hub']}")


def list_remote(api, repo_id: str, prefix: str) -> set[str]:
    from huggingface_hub import list_repo_files

    return {p for p in list_repo_files(repo_id, repo_type="dataset") if p.startswith(prefix)}


def verify_target(api, name: str, repo_id: str) -> bool:
    cfg = ARCHIVE_TARGETS[name]
    local: Path = cfg["local"]
    hub: str = cfg["hub"]
    prefix = f"{hub}/"
    local_files = iter_local_files(local, cfg.get("ignore_patterns"))
    remote = list_remote(api, repo_id, prefix)
    missing = []
    for f in local_files:
        rel = f.relative_to(local).as_posix()
        if f"{prefix}{rel}" not in remote:
            missing.append(rel)
    print(f"\n=== Verify {name} ===")
    print(f"  local:  {len(local_files):,}")
    print(f"  remote: {len(remote):,}")
    if missing:
        print(f"  MISSING: {len(missing):,}")
        for m in missing[:15]:
            print(f"    - {m}")
        return False
    print("  PASS")
    return True


def delete_interim_selective() -> None:
    import fnmatch

    interim = DATA / "interim"
    removed = 0
    for f in interim.rglob("*"):
        if not f.is_file():
            continue
        rel = f.relative_to(interim).as_posix()
        if any(fnmatch.fnmatch(rel, pat) for pat in INTERIM_IGNORE):
            continue
        f.unlink()
        removed += 1
    print(f"Deleted {removed:,} interim files (kept chr{','.join(KEEP_SCORE_CHROMS)} score pgens)")


def delete_target(name: str) -> None:
    import shutil

    cfg = ARCHIVE_TARGETS[name]
    mode = cfg["delete"]
    if mode == "tree":
        local: Path = cfg["local"]
        if local.exists():
            shutil.rmtree(local)
            local.mkdir(parents=True, exist_ok=True)
            print(f"Deleted {local}")
    elif mode == "interim_selective":
        delete_interim_selective()


def main() -> None:
    args = parse_args()
    targets = list(ARCHIVE_TARGETS) if args.target == "all" else [args.target]
    plan = write_plan(targets)
    print(f"Plan -> {plan}")

    if args.dry_run:
        for t in targets:
            upload_target(None, t, args.repo_id, dry_run=True)
        return

    api = get_api()
    if args.verify or args.delete_local:
        if not all(verify_target(api, t, args.repo_id) for t in targets):
            sys.exit(1)
        if args.delete_local:
            for t in targets:
                delete_target(t)
        return

    for t in targets:
        upload_target(api, t, args.repo_id, dry_run=False)


if __name__ == "__main__":
    main()
