#!/usr/bin/env bash
set -euo pipefail

# ---- Config (you can override via env vars) -------------------
MTK_REPO="${MTK_REPO:-https://github.com/kahypar/mt-kahypar.git}"
MTK_COMMIT="${MTK_COMMIT:-0ef674a}"    # short/long SHA ok
BUILD_TYPE="${BUILD_TYPE:-Release}"

# If your system doesn't have suitable TBB, let Mt-KaHyPar download it:
KAHYPAR_DOWNLOAD_TBB="${KAHYPAR_DOWNLOAD_TBB:-OFF}"
# Old commit *requires* hwloc (no easy OFF switch). Ensure it's installed system-wide.

# ---- Paths ----------------------------------------------------
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

WORK_DIR="${WORK_DIR:-${SCRIPT_DIR}/build-mtkahypar}"
SRC_DIR="${WORK_DIR}/mt-kahypar"
BLD_DIR="${SRC_DIR}/build"
STAGE_DIR="${WORK_DIR}/_install"   # cmake --install prefix (for stable path)
DEST_DIR="${SCRIPT_DIR}/extern/mt-kahypar-library"

echo ">>> Mt-KaHyPar local install"
echo "    repo   : ${MTK_REPO}"
echo "    commit : ${MTK_COMMIT}"
echo "    build  : ${BUILD_TYPE}"
echo "    work   : ${WORK_DIR}"
echo "    dest   : ${DEST_DIR}"
echo

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}" "${DEST_DIR}"

# ---- Clone at specific commit (with submodules) ---------------
git clone "${MTK_REPO}" "${SRC_DIR}"
git -C "${SRC_DIR}" checkout "${MTK_COMMIT}"
git -C "${SRC_DIR}" submodule update --init --recursive

# Quick sanity (old tree uses external_tools/**)
test -d "${SRC_DIR}/external_tools" || {
  echo "!!! external_tools/ missing. Submodules not fetched correctly."
  exit 1
}

# Newer compilers reject this old growt revision because
# migration_table_mapped_reference uses ref without declaring it.
GROWT_ITERATOR="${SRC_DIR}/external_tools/growt/data-structures/migration_table_iterator.hpp"
if [[ -f "${GROWT_ITERATOR}" ]] && grep -q "sref.ref.refresh()" "${GROWT_ITERATOR}" && ! grep -q "base_mapped_reference& ref;" "${GROWT_ITERATOR}"; then
  echo ">>> Patching growt migration_table_iterator.hpp for newer compilers"
  sed -i 's/: _tab(table), _version(ver), _mref(mref)/: _tab(table), _version(ver), _mref(mref), ref(_mref)/' "${GROWT_ITERATOR}"
  sed -i '/base_mapped_reference _mref;/a\    base_mapped_reference\& ref;' "${GROWT_ITERATOR}"
fi

# ---- Configure & build ----------------------------------------
mkdir -p "${BLD_DIR}"

# NOTE: This old CMakeLists enforces TBB & HWLOC presence; install them or set KAHYPAR_DOWNLOAD_TBB=ON
cmake -S "${SRC_DIR}" -B "${BLD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DKAHYPAR_DOWNLOAD_TBB="${KAHYPAR_DOWNLOAD_TBB}" \
  -DKAHYPAR_DOWNLOAD_BOOST=OFF \
  -DKAHYPAR_ENFORCE_MINIMUM_TBB_VERSION=OFF \
  -DKAHYPAR_PYTHON=OFF \
  -DMT_KAHYPAR_DISABLE_BOOST=OFF

# Build the C interface library target (name is 'mtkahypar' in this tree)
cmake --build "${BLD_DIR}" --target mtkahypar -j"$(nproc)"

# ---- Install to staging, then copy .so ------------------------
# Some old versions don’t provide install rules; ignore install errors.
cmake --install "${BLD_DIR}" --prefix "${STAGE_DIR}" || true

# Candidate locations to search for the shared lib
CANDIDATES=(
  "${STAGE_DIR}/lib/libmtkahypar.so"
  "${STAGE_DIR}/lib64/libmtkahypar.so"
  "${STAGE_DIR}/lib/x86_64-linux-gnu/libmtkahypar.so"
  "${BLD_DIR}/libmtkahypar.so"
  "${BLD_DIR}/lib/libmtkahypar.so"
  "${BLD_DIR}/mt-kahypar/libmtkahypar.so"
)

FOUND=""
for f in "${CANDIDATES[@]}"; do
  [[ -f "$f" ]] && { FOUND="$f"; break; }
done

if [[ -z "${FOUND}" ]]; then
  FOUND="$(find "${BLD_DIR}" -maxdepth 6 -type f -name 'libmtkahypar.so' | head -n1 || true)"
fi

if [[ -z "${FOUND}" ]]; then
  echo "!!! Could not locate libmtkahypar.so after build."
  echo "    Tips:"
  echo "     - Ensure hwloc & TBB dev packages are installed."
  echo "     - Inspect build logs under: ${BLD_DIR}"
  exit 1
fi

install -m 644 "${FOUND}" "${DEST_DIR}/libmtkahypar.so"
echo ">>> Installed: ${DEST_DIR}/libmtkahypar.so"
