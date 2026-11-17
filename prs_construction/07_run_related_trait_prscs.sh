#!/bin/bash

#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=10G
#SBATCH -t 2-
#SBATCH --output=logs/%x_%A_%a.log

# Locate the repository independently of the caller's working directory.
_mdd_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_mdd_repo="$_mdd_here"
while [[ ! -f "$_mdd_repo/config/paths.sh" && "$_mdd_repo" != / ]]; do
  _mdd_repo="$(dirname "$_mdd_repo")"
done
if [[ ! -f "$_mdd_repo/config/paths.sh" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  _mdd_repo="$SLURM_SUBMIT_DIR"
  while [[ ! -f "$_mdd_repo/config/paths.sh" && "$_mdd_repo" != / ]]; do
    _mdd_repo="$(dirname "$_mdd_repo")"
  done
fi
[[ -f "$_mdd_repo/config/paths.sh" ]] || { echo 'Cannot locate config/paths.sh' >&2; exit 1; }
source "$_mdd_repo/config/paths.sh"
cd "$MDD_CODE_ROOT/prs_construction"
set -euo pipefail
trait="${1:?Supply the model/group argument}" #ADHD, BIP, SCZ
group=$2 #hapmap_logbf_over4, hapmap_cr95set, hapmap

head -n 1 "${MDD_PROJECT}/05.PRS_CS_construction/data/clean_sumstat/${trait}_NGWAS.txt"
N_gwas=$(head -n 1 "${MDD_PROJECT}/05.PRS_CS_construction/data/clean_sumstat/${trait}_NGWAS.txt" )
echo "${N_gwas}"

home_dir="${PRSCS_DIR}"


"$PRSCS_PYTHON" -u "${home_dir}/PRScs.py" --ref_dir="${MDD_LD_REF}" \
--bim_prefix="${MDD_PROJECT}/05.PRS_CS_construction/data/bim/UKB_hapmap_bim" \
--sst_file="${MDD_PROJECT}/05.PRS_CS_construction/data/clean_sumstat/${trait}_${group}.txt" \
--n_gwas=${N_gwas} \
--out_dir="${MDD_PROJECT}/05.PRS_CS_construction/data/prs-cs_output_redo/${trait}_1kg_LD_ref_${group}_redo" \
--seed=0
