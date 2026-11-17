#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=20GB
#SBATCH -t 1-
#SBATCH --array=1-10
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
cd "$MDD_CODE_ROOT/elastic_net/cross_trait_variant_sets"
set -euo pipefail
"$R_BIN" 01_train_models.R "${SLURM_ARRAY_TASK_ID:?Set a split/group index or submit as an array}" "${1:?Supply the model/group argument}"
