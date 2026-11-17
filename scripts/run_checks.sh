#!/usr/bin/env bash
# Static core-release checks; no participant data or model fits.
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$repo/config/paths.sh"
exec python3 "$repo/scripts/check_release.py"
