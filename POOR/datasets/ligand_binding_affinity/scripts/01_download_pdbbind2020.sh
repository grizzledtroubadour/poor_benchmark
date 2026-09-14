#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 脚本：01_download_pdbbind2020.sh
# 用途：检查、解压 PDBbind v2020 / v2020.R1 源数据包，并生成 source_manifest.json
# 使用：
#   cd binding_affinity/
#   bash scripts/01_download_pdbbind2020.sh
#
# 支持两种数据包格式：
#   A. 旧版 PDBbind-CN（3 个独立包）：
#        - PDBbind_v2020_plain_text_index.tar.gz
#        - PDBbind_v2020_refined.tar.gz
#        - PDBbind_v2020_other_PL.tar.gz
#   B. 新版 PDBbind+ v2020.R1 免费数据包（2 个包）：
#        - index.tar.gz
#        - P-L.tar.gz   （蛋白-配体复合物，按 release year 分区）
#
# 注意：
#   原 PDBbind-CN 站点（pdbbind.org.cn）已无法访问。
#   请从 PDBbind+ 手动下载：https://www.pdbbind-plus.org.cn/download
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data"
MANIFEST="${DATA_DIR}/source_manifest.json"

mkdir -p "${DATA_DIR}"

echo "=========================================="
echo "PDBbind v2020 源数据检查 / 解压脚本"
echo "=========================================="
echo "数据目录: ${DATA_DIR}"
echo

# ---------------------------------------------------------------------------
# 检查/解压新版 PDBbind+ v2020.R1 数据包
# ---------------------------------------------------------------------------
NEW_PACKAGES=("index.tar.gz" "P-L.tar.gz")
new_found=0
for pkg in "${NEW_PACKAGES[@]}"; do
    if [[ -f "${DATA_DIR}/${pkg}" ]]; then
        new_found=$((new_found + 1))
        size=$(du -h "${DATA_DIR}/${pkg}" | cut -f1)
        echo "[已找到] ${pkg} (${size})"
    fi
done

if [[ ${new_found} -eq ${#NEW_PACKAGES[@]} ]]; then
    echo
    echo "检测到 PDBbind+ v2020.R1 数据包。"
    for pkg in "${NEW_PACKAGES[@]}"; do
        src="${DATA_DIR}/${pkg}"
        top_dir=$(tar -tzf "${src}" | head -1 | cut -d/ -f1) || true
        dst="${DATA_DIR}/${top_dir}"
        if [[ -d "${dst}" ]]; then
            echo "[跳过] ${dst}/ 已存在"
        else
            echo "[解压] ${pkg} -> ${DATA_DIR}/"
            tar -xzf "${src}" -C "${DATA_DIR}"
        fi
    done
    echo
    echo "生成 source_manifest.json..."
    python3 - <<'PY' "${DATA_DIR}" "${MANIFEST}"
import os, sys, json, hashlib
from datetime import datetime

data_dir, manifest_path = sys.argv[1], sys.argv[2]

def md5(path, chunk=8192):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            c = f.read(chunk)
            if not c: break
            h.update(c)
    return h.hexdigest()

manifest = {
    "created_at": datetime.now().isoformat(),
    "source": "PDBbind v2020.R1 (PDBbind+ demo package)",
    "download_url": "https://www.pdbbind-plus.org.cn/download",
    "note": "PDBbind-CN (pdbbind.org.cn) is no longer accessible. This is the PDBbind+ v2020.R1 free demo package, with structures re-processed by the v2024 workflow.",
    "packages": [],
}

for tar_name in ["index.tar.gz", "P-L.tar.gz"]:
    tar_path = os.path.join(data_dir, tar_name)
    pkg = {"tarball": tar_name, "tarball_exists": os.path.isfile(tar_path)}
    if pkg["tarball_exists"]:
        pkg["tarball_size_bytes"] = os.path.getsize(tar_path)
        pkg["tarball_md5"] = md5(tar_path)
    manifest["packages"].append(pkg)

# 索引文件
index_dir = os.path.join(data_dir, "index")
expected_indices = ["INDEX_general_PL.2020R1.lst", "INDEX_general_PP.2020R1.lst",
                    "INDEX_general_PN.2020R1.lst", "INDEX_general_NL.2020R1.lst"]
manifest["index_files"] = {name: os.path.isfile(os.path.join(index_dir, name)) for name in expected_indices}

# P-L 结构统计
pl_dir = os.path.join(data_dir, "P-L")
if os.path.isdir(pl_dir):
    ranges = [d for d in os.listdir(pl_dir) if os.path.isdir(os.path.join(pl_dir, d))]
    manifest["P-L"] = {
        "year_ranges": sorted(ranges),
        "complex_dirs": sum(1 for r in ranges for _ in os.listdir(os.path.join(pl_dir, r)) if os.path.isdir(os.path.join(pl_dir, r, _))),
    }

with open(manifest_path, "w") as f:
    json.dump(manifest, f, indent=2)

print(f"Manifest written to: {manifest_path}")
PY

    echo
    echo "=========================================="
    echo "源数据准备完成"
    echo "=========================================="
    cat "${MANIFEST}"
    echo
    exit 0
fi

# ---------------------------------------------------------------------------
# 检查/解压旧版 PDBbind-CN 数据包
# ---------------------------------------------------------------------------
OLD_PACKAGES=(
    "PDBbind_v2020_plain_text_index.tar.gz"
    "PDBbind_v2020_refined.tar.gz"
    "PDBbind_v2020_other_PL.tar.gz"
)
OLD_EXTRACT_DIRS=(
    "PDBbind_v2020_plain_text_index"
    "PDBbind_v2020_refined"
    "PDBbind_v2020_other_PL"
)

old_found=0
missing=()
for pkg in "${OLD_PACKAGES[@]}"; do
    if [[ -f "${DATA_DIR}/${pkg}" ]]; then
        old_found=$((old_found + 1))
        size=$(du -h "${DATA_DIR}/${pkg}" | cut -f1)
        echo "[已找到] ${pkg} (${size})"
    else
        echo "[缺失]   ${pkg}"
        missing+=("${pkg}")
    fi
done

if [[ ${old_found} -eq ${#OLD_PACKAGES[@]} ]]; then
    echo
    echo "检测到旧版 PDBbind-CN 数据包。"
    for i in "${!OLD_PACKAGES[@]}"; do
        pkg="${OLD_PACKAGES[$i]}"
        dst="${DATA_DIR}/${OLD_EXTRACT_DIRS[$i]}"
        src="${DATA_DIR}/${pkg}"
        if [[ -d "${dst}" ]]; then
            echo "[跳过] ${dst}/ 已存在"
        else
            echo "[解压] ${pkg} -> ${dst}/"
            tar -xzf "${src}" -C "${DATA_DIR}"
            top_dir=$(tar -tzf "${src}" | head -1 | cut -d/ -f1)
            if [[ "${top_dir}" != "${OLD_EXTRACT_DIRS[$i]}" && -d "${DATA_DIR}/${top_dir}" ]]; then
                mv "${DATA_DIR}/${top_dir}" "${dst}"
            fi
        fi
    done
    exit 0
fi

# ---------------------------------------------------------------------------
# 未找到任何完整数据包
# ---------------------------------------------------------------------------
echo
echo "ERROR: 未在 ${DATA_DIR} 中找到完整的数据包。"
echo
echo "请从 PDBbind+ 手动下载（推荐 v2020.R1 免费数据包）："
echo "  1. 访问 https://www.pdbbind-plus.org.cn/download"
echo "  2. 注册并登录免费账号"
echo "  3. 下载以下文件并放入 ${DATA_DIR}/："
echo "       - index.tar.gz"
echo "       - P-L.tar.gz"
echo "  4. 重新运行本脚本"
echo
exit 1
