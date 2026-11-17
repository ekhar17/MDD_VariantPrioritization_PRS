# Scope and reproducibility

This upload is limited to the core PRS construction, original elastic-net models
and original figure scripts selected for sharing. The broader benchmark,
sensitivity and layout-assembly modules remain outside this upload.

## What changed during preparation

- Original research sources remain unchanged. `SCRIPT_PROVENANCE.json` records
  their source paths and hashes and the cleaned release hashes.
- Duplicated files, reduced-CV QUICK variants, scratch work and obsolete
  output-version comparisons were excluded.
- Original machine paths were replaced with shared configuration; site-specific
  module loads and partition settings were removed from launchers.
- Undefined error sentinels were replaced with explicit `stop()` calls and
  plotting imports were made explicit.
- The old stage-04 formatting destination was normalized to the stage-05 data
  paths consumed downstream, and `_hapamap.txt` was corrected to `_hapmap.txt`.
- Incomplete interactive plotting tails were removed while retaining complete
  saved plots.
- Code names were made descriptive, then execution-order prefixes were added.
  Launchers and references were updated together. `RENAMING_MAP.json` records
  the code-name mapping; `RUN_ORDER.json` records the numbered steps.

Numbering changes do not alter model formulas, hyperparameter grids, seeds,
CV settings, sample memberships, scoring conventions or data/output names.
Matched R/shell files represent alternative ways of running the same step.

## Required precomputed inputs

The original cohort QC/ancestry/split generation, initial fGWAS/TWAS
prioritization, MetaSubtract and stratified-LDSC scripts were not found in the
source collection. Their saved outputs remain required. Some plotting
intermediate generators are also missing; see [INPUTS.md](INPUTS.md).

Exact historical package versions and the original PRS-CS commit are not fully
recorded. No original lockfile has been invented.

## Checks and limits

`bash scripts/run_checks.sh` checks syntax, source-only contents, private-path
patterns, common credential patterns, provenance hashes and code-map targets.
It requires R, Bash and Python 3, and does not fit models or read participant
data. Packaging checks also confirm launcher argument forwarding and stable
external data/output paths.

The full UK Biobank analyses were not rerun during packaging. These checks do
not establish numerical reproduction of all published results or final layouts.
The separate synthetic suites previously used for additional analyses are not
part of this selected core upload; they are not claimed as validation of the
original core models.
