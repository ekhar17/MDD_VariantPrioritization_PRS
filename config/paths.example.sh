#!/usr/bin/env bash
# Copy to config/paths.local.sh and replace the example values.
# Keep the local file and all participant data outside the GitHub upload.
export MDD_PROJECT="/absolute/path/to/private/analysis_workspace"
export MDD_SOURCE_DATA="/absolute/path/to/source_data"
export MDD_HAPMAP="/absolute/path/to/snpinfo_1kg_hm3"
export MDD_LD_REF="/absolute/path/to/ldblk_1kg_eur"
# .rds and its companion .bk backing file must remain together.
export MDD_BACKING="/absolute/path/to/genotypes.rds"
export MDD_SAMPLE="/absolute/path/to/ukb.sample"
# Optional: use already-installed, compatible interpreters/tools.
# export PRSCS_DIR="/absolute/path/to/PRScs"
# export PRSCS_PYTHON="/absolute/path/to/prscs_environment/bin/python"
# export R_BIN="/absolute/path/to/Rscript"
