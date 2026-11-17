#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH -t 1-
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
"$R_BIN" 23_calculate_allele_frequencies.R
