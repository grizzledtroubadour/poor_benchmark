#!/usr/bin/env python3
"""
下载 DIPS-Plus v1.3.0 数据集归档。

用法：
    python scripts/download_dips_plus.py [--verify] [--msa]

默认仅下载核心归档 final_raw_dips.tar.gz；加 --msa 同时下载外部特征与 MSA 归档。
"""

import argparse
import hashlib
import json
import logging
import os
import subprocess
import sys
from pathlib import Path
from urllib.request import urlopen

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/download_dips_plus.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ZENODO_RECORD = "8140981"
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"

FILES = {
    "final_raw_dips.tar.gz": {
        "description": "核心数据集归档（pair 文件、HDF5、划分列表、元数据）",
        "required": True,
    },
    "interim_external_feats_dips_msas.tar.gz": {
        "description": "外部特征与 MSA 归档（可选）",
        "required": False,
    },
}


def fetch_record_metadata():
    """从 Zenodo API 获取记录元数据。"""
    logger.info("Fetching Zenodo record metadata from %s", ZENODO_API)
    with urlopen(ZENODO_API, timeout=60) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    files = {f["key"]: f for f in data.get("files", [])}
    logger.info("Found %d files in Zenodo record", len(files))
    return files


def download_file(filename, files_meta):
    """使用 aria2c 下载单个文件，支持断点续传。"""
    if filename not in files_meta:
        raise ValueError(f"{filename} not found in Zenodo record")

    meta = files_meta[filename]
    url = meta["links"]["self"]
    expected_size = meta["size"]
    expected_md5 = meta["checksum"].replace("md5:", "")

    out_path = DATA_DIR / filename
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if out_path.exists():
        actual_size = out_path.stat().st_size
        if actual_size == expected_size:
            logger.info("%s already exists with correct size (%d bytes)", filename, expected_size)
            return out_path, expected_md5
        logger.info("%s exists but size mismatch (%d vs %d); resuming", filename, actual_size, expected_size)

    logger.info("Downloading %s from %s", filename, url)
    cmd = [
        "aria2c",
        "-x", "4",
        "-s", "4",
        "-c",
        "--max-tries=10",
        "--retry-wait=30",
        "-d", str(DATA_DIR),
        "-o", filename,
        url,
    ]
    subprocess.run(cmd, check=True)
    return out_path, expected_md5


def md5sum(path):
    """计算文件 MD5。"""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_file(path, expected_md5, expected_size):
    """校验文件大小与 MD5。"""
    actual_size = path.stat().st_size
    if actual_size != expected_size:
        raise ValueError(f"Size mismatch for {path.name}: {actual_size} != {expected_size}")

    actual_md5 = md5sum(path)
    if actual_md5 != expected_md5:
        raise ValueError(f"MD5 mismatch for {path.name}: {actual_md5} != {expected_md5}")

    logger.info("Verification passed for %s", path.name)


def main():
    parser = argparse.ArgumentParser(description="Download DIPS-Plus v1.3.0")
    parser.add_argument("--verify", action="store_true", help="下载完成后校验 MD5")
    parser.add_argument("--msa", action="store_true", help="同时下载外部特征与 MSA 归档")
    args = parser.parse_args()

    files_meta = fetch_record_metadata()

    selected = ["final_raw_dips.tar.gz"]
    if args.msa:
        selected.append("interim_external_feats_dips_msas.tar.gz")

    for filename in selected:
        desc = FILES[filename]["description"]
        logger.info("=" * 60)
        logger.info("File: %s", filename)
        logger.info("Description: %s", desc)

        path, expected_md5 = download_file(filename, files_meta)

        if args.verify:
            expected_size = files_meta[filename]["size"]
            verify_file(path, expected_md5, expected_size)

    logger.info("=" * 60)
    logger.info("Download stage completed")


if __name__ == "__main__":
    main()
