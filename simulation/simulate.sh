#!/usr/bin/env bash
# simulate.sh - simulate one dataset of the QRAFT paper with SimPhy-X (or SimPhy 1.0.2).
#
# Adapted from run_simulator.sh by Anik Saha (the benchmark harness of the STELAR-X and ASTRAL-X
# studies); the SimPhy command line is unchanged.
#
# Writes <data_dir>/t_T_g_G_sb_SB_spmin_MIN_spmax_MAX/R<n>/{all_gt.tre,s_tree.trees} for n = 1..replicates:
# all_gt.tre holds one gene tree per line, s_tree.trees the species tree. SimPhy's own .command, .params and
# .db files stay in the dataset folder as the record of how it was made.
#
# Engines (--engine):
#   simphy-x-compat  (default) SimPhy-X in compat mode: the same trees as SimPhy 1.0.2, multithreaded.
#   simphy-x-fast    SimPhy-X fast mode: the same model with independent per-locus random streams, so the
#                    trees are not those of SimPhy 1.0.2.
#   simphy-old       SimPhy 1.0.2 itself (single-threaded).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON:-python3}"

sb="0.000001"
spmin="500000"
spmax="1500000"
seed="42"
out_dir=""
data_base_dir="data"
replicates="10"
engine="simphy-x-compat"
threads="all"
simphy_bin=""
taxa_num=""
gene_trees=""

print_usage() {
  cat <<USAGE
Usage: $0 -t <taxa> -g <gene_trees> [options]

Required:
  -t, --taxa_num    number of taxa (e.g. 1000)
  -g, --gene_trees  number of gene trees (e.g. 1000)

Optional:
  -o, --out_dir     output folder (default: built from the parameters under --data_dir)
  -d, --data_dir    base data folder (default: ${data_base_dir})
      --replicates  number of replicates (default: ${replicates})
      --sb          speciation rate (default: ${sb})
      --spmin       population size minimum (default: ${spmin})
      --spmax       population size maximum (default: ${spmax})
      --seed        SimPhy random seed (default: ${seed})
      --engine      simphy-x-compat | simphy-x-fast | simphy-old (default: ${engine})
      --threads     threads for the SimPhy-X engines, N or all (default: ${threads})
      --simphy-bin  the simulator binary (default: simphy-x on PATH, or simphy_lnx64 for simphy-old)
  -h, --help        show this help and exit

Example (1,000 taxa, 1,000 gene trees, ILS-L2, five replicates, as in the paper):
  $0 -t 1000 -g 1000 --replicates 5 --spmin 100000 --spmax 200000 --simphy-bin ./simphy-x
USAGE
}

OPTS=$(getopt -o t:g:o:d:h --long taxa_num:,taxa-num:,gene_trees:,gene-trees:,out_dir:,out-dir:,data_dir:,data-dir:,replicates:,sb:,spmin:,spmax:,seed:,engine:,threads:,simphy-bin:,help -n 'simulate.sh' -- "$@") || {
  echo "Failed parsing options." >&2
  exit 1
}
eval set -- "$OPTS"

while true; do
  case "$1" in
    -t|--taxa_num|--taxa-num) taxa_num="$2"; shift 2 ;;
    -g|--gene_trees|--gene-trees) gene_trees="$2"; shift 2 ;;
    -o|--out_dir|--out-dir) out_dir="$2"; shift 2 ;;
    -d|--data_dir|--data-dir) data_base_dir="$2"; shift 2 ;;
    --replicates) replicates="$2"; shift 2 ;;
    --sb) sb="$2"; shift 2 ;;
    --spmin) spmin="$2"; shift 2 ;;
    --spmax) spmax="$2"; shift 2 ;;
    --seed) seed="$2"; shift 2 ;;
    --engine) engine="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --simphy-bin) simphy_bin="$2"; shift 2 ;;
    -h|--help) print_usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Unexpected argument to simulate.sh: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "${taxa_num}" || -z "${gene_trees}" ]]; then
  echo "Error: --taxa_num (-t) and --gene_trees (-g) are required." >&2
  print_usage
  exit 1
fi
for item in "taxa_num:$taxa_num" "gene_trees:$gene_trees" "replicates:$replicates" "spmin:$spmin" "spmax:$spmax" "seed:$seed"; do
  if [[ ! "${item#*:}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: ${item%%:*} must be a positive integer (got '${item#*:}')." >&2
    exit 1
  fi
done
if (( spmin > spmax )); then
  echo "Error: --spmin cannot exceed --spmax." >&2
  exit 1
fi

engine_args=()
case "$engine" in
  simphy-x-compat) default_bin="simphy-x"; engine_args=(-XM compat) ;;
  simphy-x-fast)   default_bin="simphy-x"; engine_args=(-XM fast) ;;
  simphy-old|simphy) engine="simphy-old"; default_bin="simphy_lnx64"; engine_args=() ;;
  *)
    echo "Error: unknown --engine '$engine' (expected simphy-x-compat, simphy-x-fast or simphy-old)." >&2
    exit 1
    ;;
esac
if [[ "$engine" != simphy-old && "$threads" != all ]]; then
  if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --threads must be a positive integer or 'all' (got '$threads')." >&2
    exit 1
  fi
  engine_args+=(-XT "$threads")
fi
if [[ -z "$simphy_bin" ]]; then
  simphy_bin="$(command -v "$default_bin" || true)"
  if [[ -z "$simphy_bin" ]]; then
    echo "Error: $default_bin is not on PATH; give the binary with --simphy-bin (see README, Simulating the datasets)." >&2
    exit 1
  fi
fi
if [[ ! -x "$simphy_bin" ]]; then
  echo "Error: SimPhy binary not found or not executable: $simphy_bin" >&2
  exit 1
fi

if [[ -z "${out_dir}" ]]; then
  out_dir="${data_base_dir%/}/t_${taxa_num}_g_${gene_trees}_sb_${sb}_spmin_${spmin}_spmax_${spmax}"
else
  out_dir="${out_dir%/}"
fi
mkdir -p "${out_dir}"

echo "Running with parameters:"
echo "  taxa_num      = ${taxa_num}"
echo "  gene_trees    = ${gene_trees}"
echo "  replicates    = ${replicates}"
echo "  sb            = ${sb}"
echo "  spmin         = ${spmin}"
echo "  spmax         = ${spmax}"
echo "  seed          = ${seed}"
echo "  engine        = ${engine} (${simphy_bin}${engine_args[0]+ ${engine_args[*]}})"
echo "  out_dir       = ${out_dir}"
echo ""

# SimPhy parameters
#   -sb speciation rate            -ld/-lb/-lt loss/duplication/HGT rates (off)
#   -rs replicates                 -rl gene trees per species tree
#   -rg gene trees per locus tree  -sp population size ~ U(spmin, spmax)
#   -su substitution rate ~ LN     -sg generation time (fixed 1)
#   -sl taxa                       -st species tree height ~ LN(16.2, 1)
#   -om/-od/-op/-oc/-on output options (mapping, SQLite db, params, command)
#   -cs seed
"${simphy_bin}" \
  ${engine_args[@]+"${engine_args[@]}"} \
  -sb "f:${sb}" \
  -ld f:0 \
  -lb f:0 \
  -lt f:0 \
  -rs "${replicates}" \
  -rl "f:${gene_trees}" \
  -rg 1 \
  -o "${out_dir}" \
  -sp "u:${spmin},${spmax}" \
  -su ln:-17.27461,0.6931472 \
  -sg f:1 \
  -sl "f:${taxa_num}" \
  -st ln:16.2,1 \
  -om 1 \
  -v 2 \
  -od 1 \
  -op 1 \
  -oc 1 \
  -on 1 \
  -cs "${seed}"

# Replicate folders without leading zeros.
for dir in "${out_dir}"/*/; do
  [[ -d "$dir" ]] || continue
  dirname=$(basename "$dir")
  if [[ "$dirname" =~ ^0*([1-9][0-9]*|0)$ ]]; then
    canonical_name="${BASH_REMATCH[1]}"
    if [[ "$dirname" != "$canonical_name" ]]; then
      mv "${out_dir}/${dirname}" "${out_dir}/${canonical_name}"
    fi
  fi
done

# The per-locus gene trees of every replicate, one per line, in all_gt.tre.
for i in $(seq 1 "${replicates}"); do
  if [[ -d "${out_dir}/${i}" ]]; then
    "$PYTHON_BIN" "${SCRIPT_DIR}/concat_gene_trees.py" "${out_dir}/${i}"
  else
    echo "Warning: replicate folder ${out_dir}/${i} not found." >&2
  fi
done

# all_gt.tre and s_tree.trees into R<n>/; the raw numbered folders are removed.
"$PYTHON_BIN" "${SCRIPT_DIR}/reorganize_trees.py" "${out_dir}"

echo "Done. Output in ${out_dir}"
