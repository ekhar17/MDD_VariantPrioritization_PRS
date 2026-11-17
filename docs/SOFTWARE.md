# Software

The manuscript reports R 4.1.0 for the original analysis. No complete original
lockfile or session record was supplied, so historical package versions cannot
be fully recovered from the imports.

Run `bash scripts/run_r.sh scripts/01_check_dependencies.R` to inspect your R
environment. The command reports installed versions; it does not install them.

| Purpose | Packages/tools |
| --- | --- |
| Data preparation and scoring | `data.table`, `tidyverse`, `purrr`, `bigsnpr`, `bigstatsr` |
| Regression and evaluation | `fmsb`, `pscl`, `caret`, `glmnet`, `glmnetUtils`; some environments also use `pROC` |
| Figures | `ggplot2`, `ggpattern`, `colorspace`, `UpSetR`, `ComplexUpset`; retained scripts import the full `tidyverse` |
| PRS-CS fitting | PRS-CS, a compatible Python environment, and the specified reference LD data |
| Launching/checking | Bash; Python 3 for `scripts/check_release.py`; SLURM only for the batch examples |

Install PRS-CS separately from <https://github.com/getian107/PRScs>. Set
`PRSCS_DIR` to its directory and `PRSCS_PYTHON` to a compatible interpreter in
`config/paths.local.sh`. Defaults are `$MDD_TOOLS/PRScs` and `python3`.
The original scripts did not record the exact PRS-CS commit. Retain the version
and environment information from your actual run.

The original LDSC enrichment/prioritization steps require their own software
and environment. Their scripts are outside this selected upload; supply their
original outputs as described in [INPUTS.md](INPUTS.md).

Full genotype scoring and repeated model fitting can require substantial RAM,
storage and cluster time. Scheduler partition/module names are not hardcoded;
load the required tools for your own environment before submitting a batch job.
