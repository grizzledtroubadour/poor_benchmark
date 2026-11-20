#!/usr/bin/env python3
"""
Extract DIPS-Plus archives and inspect their directory layout.

Usage:
    python scripts/extract_and_inspect.py [--core-only] [--dry-run]

--core-only: extract only core contents (pairs, split lists, metadata) from final_raw_dips.tar.gz.
--dry-run: list archive contents without extracting.
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/extract_and_inspect.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"

CORE_ARCHIVE = DATA_DIR / "final_raw_dips.tar.gz"
MSA_ARCHIVE = DATA_DIR / "interim_external_feats_dips_msas.tar.gz"


def run_cmd(cmd, cwd=None):
    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("Command failed: %s", result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result.stdout


def list_archive(archive_path):
    """List tar.gz archive contents."""
    if not archive_path.exists():
        logger.warning("Archive not found: %s", archive_path)
        return []
    stdout = run_cmd(["tar", "-tzf", str(archive_path)])
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    return lines


def extract_archive(archive_path, dest_dir, patterns=None):
    """Extract a tar.gz archive, optionally filtering by pattern."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["tar", "-xzf", str(archive_path), "-C", str(dest_dir)]
    if patterns:
        cmd.extend(["--wildcards"] + patterns)
    run_cmd(cmd)
    logger.info("Extracted %s to %s", archive_path.name, dest_dir)


def build_directory_tree(paths):
    """Build a directory-tree summary from an archive listing."""
    tree = {}
    for p in paths:
        parts = p.split("/")
        node = tree
        for part in parts:
            if not part:
                continue
            if part not in node:
                node[part] = {}
            node = node[part]
    return tree


def tree_summary(tree, prefix=""):
    """Convert a tree to a printable string."""
    lines = []
    keys = sorted(tree.keys())
    for i, key in enumerate(keys):
        is_last = i == len(keys) - 1
        lines.append(f"{prefix}{'└── ' if is_last else '├── '}{key}")
        if tree[key]:
            lines.extend(tree_summary(tree[key], prefix + ("    " if is_last else "│   ")))
    return lines


def main():
    parser = argparse.ArgumentParser(description="Extract and inspect DIPS-Plus archives")
    parser.add_argument("--core-only", action="store_true", help="Extract core archive only")
    parser.add_argument("--dry-run", action="store_true", help="List archive contents only")
    parser.add_argument("--patterns", nargs="+", help="Extract only matching file patterns (e.g. *.txt)")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summaries = {}

    # Core archive
    logger.info("Inspecting core archive: %s", CORE_ARCHIVE)
    core_paths = list_archive(CORE_ARCHIVE)
    summaries["core"] = {
        "archive": CORE_ARCHIVE.name,
        "total_entries": len(core_paths),
        "top_level_dirs": sorted(set(p.split("/")[0] for p in core_paths if "/" in p)),
    }
    logger.info("Core archive entries: %d", len(core_paths))

    if args.dry_run:
        tree = build_directory_tree(core_paths[:5000])
        tree_text = "\n".join(tree_summary(tree))
        (OUTPUT_DIR / "dips_plus_structure_tree.txt").write_text(tree_text)
        logger.info("Dry-run tree saved to %s", OUTPUT_DIR / "dips_plus_structure_tree.txt")
    else:
        patterns = args.patterns if args.patterns else None
        extract_archive(CORE_ARCHIVE, DATA_DIR, patterns=patterns)

    # Optional MSA archive
    if not args.core_only and MSA_ARCHIVE.exists():
        logger.info("Inspecting MSA archive: %s", MSA_ARCHIVE)
        msa_paths = list_archive(MSA_ARCHIVE)
        summaries["msa"] = {
            "archive": MSA_ARCHIVE.name,
            "total_entries": len(msa_paths),
            "top_level_dirs": sorted(set(p.split("/")[0] for p in msa_paths if "/" in p)),
        }
        logger.info("MSA archive entries: %d", len(msa_paths))
        if not args.dry_run:
            extract_archive(MSA_ARCHIVE, DATA_DIR)

    (OUTPUT_DIR / "extraction_summary.json").write_text(json.dumps(summaries, indent=2))
    logger.info("Summary saved to %s", OUTPUT_DIR / "extraction_summary.json")


if __name__ == "__main__":
    main()
