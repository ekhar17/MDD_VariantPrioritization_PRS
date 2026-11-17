#!/usr/bin/env bash
# Shared defaults. Put machine-specific overrides in ignored paths.local.sh.
export MDD_CODE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -f "$MDD_CODE_ROOT/config/paths.local.sh" ]]; then
  source "$MDD_CODE_ROOT/config/paths.local.sh"
fi
export MDD_PROJECT="${MDD_PROJECT:-$MDD_CODE_ROOT/workspace}"
export MDD_SOURCE_DATA="${MDD_SOURCE_DATA:-$MDD_PROJECT/source_data}"
export MDD_REFERENCES="${MDD_REFERENCES:-$MDD_PROJECT/references}"
export MDD_TOOLS="${MDD_TOOLS:-$MDD_PROJECT/tools}"
export MDD_REPORTS="${MDD_REPORTS:-$MDD_PROJECT/reports}"
export MDD_HAPMAP="${MDD_HAPMAP:-$MDD_REFERENCES/snpinfo_1kg_hm3}"
export MDD_LD_REF="${MDD_LD_REF:-$MDD_REFERENCES/ldblk_1kg_eur}"
export MDD_BACKING="${MDD_BACKING:-$MDD_PROJECT/genotypes/imputed_bgen_hapmap3_allind_mult_1kg.rds}"
export MDD_SAMPLE="${MDD_SAMPLE:-$MDD_PROJECT/genotypes/ukb.sample}"
export R_BIN="${R_BIN:-Rscript}"
export PRSCS_DIR="${PRSCS_DIR:-$MDD_TOOLS/PRScs}"
export PRSCS_PYTHON="${PRSCS_PYTHON:-python3}"
