#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$repo/config/paths.sh"
if [[ $# -lt 1 ]]; then
  echo 'Usage: bash scripts/run_r.sh REPOSITORY_RELATIVE_SCRIPT.R [arguments...]' >&2
  exit 2
fi
target="$1"; shift
[[ "$target" != /* && "$target" != *..* && "$target" == *.R ]] || {
  echo 'Supply an R script path relative to this repository.' >&2; exit 2;
}
[[ -f "$repo/$target" ]] || { echo "Script not found: $target" >&2; exit 2; }
cd "$repo"
exec "$R_BIN" --vanilla "$repo/$target" "$@"
