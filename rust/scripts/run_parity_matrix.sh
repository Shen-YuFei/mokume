#!/usr/bin/env bash
# Run the features2proteins parity matrix: for each dataset x method combo, run
# the Python mokume baseline and the Rust port, then diff the protein matrices
# with compare_protein_matrices.py. Reports are written under
# ${OUT_ROOT}/<PXD>/ and never committed to the repo.
#
# Usage:
#   bash scripts/run_parity_matrix.sh [PXD ...]
# With no arguments, runs the four Phase-1 datasets. Override behaviour via env:
#   DATA_ROOT, OUT_ROOT, PY_MOKUME, RUST_BIN, THREADS, COMBOS
#
# All multi-threaded steps use 24 threads per project convention.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${DATA_ROOT:-/home/shenyufei/Git-repository/Bigbio/Bigbio_data}"
OUT_ROOT="${OUT_ROOT:-/tmp/mokume-parity}"
THREADS="${THREADS:-24}"
PY_MOKUME="${PY_MOKUME:-conda run -n Bigbio python -m mokume.mokume_cli}"
RUST_BIN="${RUST_BIN:-${REPO_ROOT}/target/release/mokume}"
COMPARE="${REPO_ROOT}/scripts/compare_protein_matrices.py"

# Datasets: CLI args override the default Phase-1 set.
if [ "$#" -gt 0 ]; then
  DATASETS=("$@")
else
  DATASETS=(PXD003539 PXD004701 PXD041421 PXD030304)
fi

# Method combos: "name|quant|run_norm|sample_norm". Phase-1 default set.
DEFAULT_COMBOS=(
  "sum_none_none|sum|none|none"
  "topn_none_none|topn|none|none"
  "median_none_none|median|none|none"
  "sum_runmedian_none|sum|median|none"
  "sum_none_globalmedian|sum|none|globalmedian"
  "sum_none_conditionmedian|sum|none|conditionmedian"
)
if [ -n "${COMBOS:-}" ]; then
  read -r -a COMBO_LIST <<<"${COMBOS}"
else
  COMBO_LIST=("${DEFAULT_COMBOS[@]}")
fi

echo "[parity] building Rust release binary (one-off, for real-data speed)"
if ! cargo build --release -p mokume-cli --jobs "${THREADS}"; then
  echo "[parity] FATAL: cargo build failed" >&2
  exit 1
fi

run_one() {
  local pxd="$1" combo="$2"
  local name quant run_norm sample_norm
  IFS='|' read -r name quant run_norm sample_norm <<<"${combo}"

  local parquet="${DATA_ROOT}/cell-lines/${pxd}/qpx/${pxd}.feature.parquet"
  local sdrf="${DATA_ROOT}/cell_lines_mokume/${pxd}/sdrf/${pxd}.sdrf.tsv"
  if [ ! -f "${parquet}" ]; then
    echo "[parity][${pxd}] SKIP ${name}: missing parquet ${parquet}"
    return 0
  fi

  local out_dir="${OUT_ROOT}/${pxd}"
  mkdir -p "${out_dir}/python" "${out_dir}/rust"
  local py_csv="${out_dir}/python/${name}.csv"
  local rs_csv="${out_dir}/rust/${name}.csv"
  local report="${out_dir}/${name}_report.json"

  local sdrf_py=() sdrf_rs=()
  if [ -f "${sdrf}" ]; then
    sdrf_py=(--sdrf "${sdrf}")
    sdrf_rs=(--sdrf "${sdrf}")
  fi

  echo "[parity][${pxd}] ${name}: python baseline"
  # shellcheck disable=SC2086
  ${PY_MOKUME} features2proteins \
    --parquet "${parquet}" "${sdrf_py[@]}" --output "${py_csv}" \
    --quant-method "${quant}" \
    --run-normalization "${run_norm}" \
    --sample-normalization "${sample_norm}" \
    --duckdb-threads "${THREADS}" \
    || { echo "[parity][${pxd}] ${name}: python FAILED" >&2; return 1; }

  echo "[parity][${pxd}] ${name}: rust"
  "${RUST_BIN}" features2proteins \
    --parquet "${parquet}" "${sdrf_rs[@]}" --output "${rs_csv}" \
    --quant-method "${quant}" \
    --run-normalization "${run_norm}" \
    --sample-normalization "${sample_norm}" \
    --threads "${THREADS}" \
    || { echo "[parity][${pxd}] ${name}: rust FAILED" >&2; return 1; }

  echo "[parity][${pxd}] ${name}: compare -> ${report}"
  conda run -n Bigbio python "${COMPARE}" \
    --python "${py_csv}" --rust "${rs_csv}" --output "${report}" \
    || { echo "[parity][${pxd}] ${name}: compare FAILED" >&2; return 1; }
}

for pxd in "${DATASETS[@]}"; do
  for combo in "${COMBO_LIST[@]}"; do
    run_one "${pxd}" "${combo}"
  done
done

echo "[parity] done. Reports under ${OUT_ROOT}/<PXD>/<combo>_report.json"
