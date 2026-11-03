#!/usr/bin/env bash
# paper-datasets.sh - the simulated datasets of the QRAFT paper, one family at a time.
#
#   simulation/paper-datasets.sh FAMILY [--simphy-bin PATH] [--data-dir DIR] [--threads N] [--dry-run]
#
# FAMILY:
#   scalability  ILS-L2, five replicates: 1,000-300,000 taxa with 1,000 gene trees, and 1,000 taxa with
#                1,000-300,000 gene trees (the STELAR-X / ASTRAL-X datasets)
#   million      ILS-L2, five replicates: 100,000-1,000,000 taxa with 1,000 gene trees, and 1,000 taxa with
#                100,000-1,000,000 gene trees, in steps of 100,000
#   vary-ils     ILS-L1 to ILS-L4, five replicates: 1,000-40,000 taxa with 1,000 gene trees, and 1,000 taxa with
#                5,000-50,000 gene trees
#
# Every dataset uses seed 42, so the same command gives the same files as the paper's. --dry-run prints the
# simulate.sh commands without running them. Datasets already present (all R<n>/all_gt.tre) are skipped.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
family="${1:-}"
[[ -n "$family" ]] && shift || true
simphy_bin=""
data_dir="data"
threads="all"
dry_run=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --simphy-bin) simphy_bin="$2"; shift 2 ;;
    --data-dir) data_dir="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "Unexpected argument: $1" >&2; exit 1 ;;
  esac
done

replicates=5
spmin=100000

one() {   # one TAXA GENES SPMAX
  local t=$1 g=$2 hi=$3
  local dir="${data_dir%/}/t_${t}_g_${g}_sb_0.000001_spmin_${spmin}_spmax_${hi}"
  local cmd=("$SCRIPT_DIR/simulate.sh" -t "$t" -g "$g" --replicates "$replicates" --spmin "$spmin" --spmax "$hi"
             --seed 42 --data-dir "$data_dir" --threads "$threads")
  [[ -n "$simphy_bin" ]] && cmd+=(--simphy-bin "$simphy_bin")
  if (( dry_run )); then
    echo "${cmd[*]}"
    return
  fi
  local done_reps=0
  for r in $(seq 1 "$replicates"); do [[ -s "$dir/R$r/all_gt.tre" ]] && done_reps=$((done_reps + 1)); done
  if (( done_reps == replicates )); then
    echo "skip $dir (present)"
    return
  fi
  "${cmd[@]}"
}

case "$family" in
  scalability)
    for n in 1000 2500 5000 7500 10000 25000 50000 75000 100000 125000 150000 175000 200000 225000 250000 \
             275000 300000; do
      one "$n" 1000 200000
      [[ "$n" = 1000 ]] || one 1000 "$n" 200000
    done
    ;;
  million)
    for n in $(seq 100000 100000 1000000); do
      one "$n" 1000 200000
      one 1000 "$n" 200000
    done
    ;;
  vary-ils)
    for hi in 150000 200000 250000 300000; do   # ILS-L1, L2, L3, L4
      for t in 1000 5000 10000 20000 30000 40000; do one "$t" 1000 "$hi"; done
      for g in 5000 10000 25000 50000; do one 1000 "$g" "$hi"; done
    done
    ;;
  *)
    echo "Usage: $0 scalability|million|vary-ils [--simphy-bin PATH] [--data-dir DIR] [--threads N] [--dry-run]" >&2
    exit 1
    ;;
esac
