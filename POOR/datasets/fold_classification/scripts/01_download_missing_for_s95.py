#!/usr/bin/env python3
"""
下载S95结构域提取所需的缺失PDB文件。
优先 .pdb.gz，回退 .cif.gz

缺失清单由本脚本自动生成：S95 列表的唯一 PDB ID 与 data/pdb/ 本地索引对比，
写入 output/missing_for_s95.txt 后执行下载。
"""

import os
import time
import gzip
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

TASK_DIR = Path(__file__).resolve().parents[1]  # datasets/fold_classification
ROOT = Path(__file__).resolve().parents[3]      # 仓库根（POOR）
DATA_PDB = ROOT / "data" / "pdb"
CATH_DOMAIN_LIST_S95 = (
    ROOT / "data" / "CATHv44" / "cath-classification-data"
    / "cath-domain-list-S95.txt"
)
# 缺失 PDB ID 清单（每行一个 4 字符 pdb id）；由本脚本自动生成
MISSING_IDS_FILE = TASK_DIR / "output" / "missing_for_s95.txt"
LOG_FILE = TASK_DIR / "output" / "download_missing_s95.log"

PDB_URL = "https://files.rcsb.org/download/{}.pdb.gz"
CIF_URL = "https://files.rcsb.org/download/{}.cif.gz"
MAX_RETRIES = 3
REQUEST_DELAY = 0.2


def log(msg):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


def generate_missing_ids_file():
    """对比 S95 列表的唯一 PDB ID 与 data/pdb/ 本地索引，生成缺失清单。"""
    s95_pdb_ids = set()
    with open(CATH_DOMAIN_LIST_S95) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tokens = line.split()
            if len(tokens) < 12:
                continue
            s95_pdb_ids.add(tokens[0][:4].lower())

    local_ids = set()
    for f in DATA_PDB.iterdir():
        if f.is_file() and f.suffix == '.gz':
            local_ids.add(f.name.split('.')[0].lower()[:4])

    missing = sorted(s95_pdb_ids - local_ids)
    MISSING_IDS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(MISSING_IDS_FILE, "w") as f:
        for pid in missing:
            f.write(pid + "\n")
    log(f"S95 唯一 PDB ID: {len(s95_pdb_ids)} | 本地已有: "
        f"{len(s95_pdb_ids & local_ids)} | 缺失: {len(missing)}")
    log(f"缺失清单已生成: {MISSING_IDS_FILE}")


def download_one(pdb_id: str):
    pdb_id = pdb_id.lower().strip()
    
    pdb_path = DATA_PDB / f"{pdb_id}.pdb.gz"
    cif_path = DATA_PDB / f"{pdb_id}.cif.gz"
    
    if pdb_path.exists() or cif_path.exists():
        return pdb_id, "already_exists", None
    
    for url_template, out_path, ext in [(PDB_URL, pdb_path, ".pdb.gz"), (CIF_URL, cif_path, ".cif.gz")]:
        url = url_template.format(pdb_id)
        for attempt in range(MAX_RETRIES):
            try:
                req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urlopen(req, timeout=60) as response:
                    data = response.read()
                    try:
                        gzip.decompress(data)
                    except gzip.BadGzipFile:
                        continue
                    with open(out_path, "wb") as f:
                        f.write(data)
                    return pdb_id, "success", ext
            except HTTPError as e:
                if e.code == 404:
                    break
                time.sleep(REQUEST_DELAY * (attempt + 1))
            except Exception:
                time.sleep(REQUEST_DELAY * (attempt + 1))
    
    return pdb_id, "failed", None


def main():
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    generate_missing_ids_file()

    with open(MISSING_IDS_FILE) as f:
        missing_ids = [line.strip() for line in f if line.strip()]
    
    total = len(missing_ids)
    log(f"开始下载 {total} 个缺失的PDB文件...")
    
    success = 0
    failed = 0
    skipped = 0
    pdb_count = 0
    cif_count = 0
    failed_ids = []
    
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(download_one, pid): pid for pid in missing_ids}
        for future in as_completed(futures):
            pdb_id, status, ext = future.result()
            if status == "success":
                success += 1
                if ext == ".pdb.gz":
                    pdb_count += 1
                else:
                    cif_count += 1
            elif status == "already_exists":
                skipped += 1
            else:
                failed += 1
                failed_ids.append(pdb_id)
            
            done = success + failed + skipped
            if done % 10 == 0 or done == total:
                log(f"进度: {done}/{total} | 成功: {success} (pdb:{pdb_count}, cif:{cif_count}) | 失败: {failed} | 跳过: {skipped}")
    
    if failed_ids:
        failed_file = TASK_DIR / "output" / "download_failed_s95.txt"
        with open(failed_file, "w") as f:
            for fid in failed_ids:
                f.write(fid + "\n")
        log(f"失败列表已保存: {failed_file}")
    
    log(f"完成! 总计: {total} | 成功: {success} | 失败: {failed} | 跳过: {skipped}")


if __name__ == "__main__":
    main()
