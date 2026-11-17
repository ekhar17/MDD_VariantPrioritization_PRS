# Run order

Numbers indicate execution order **within each folder**. Run PRS construction
before elastic-net models, then generate plots from completed results. A matching
`.R` and `.sh` pair is one step: run R directly or use its shell/batch launcher.
Wait for all jobs at a step to finish before starting a step that reads their
outputs. Numbering is a guide, not an automatic job dependency system.

## Setup

From the repository root:

```bash
cp config/paths.example.sh config/paths.local.sh
# Edit paths.local.sh to point to your authorized data and installed tools.
bash scripts/run_r.sh scripts/01_check_dependencies.R
bash scripts/run_r.sh scripts/02_init_workspace.R
```

The utility `scripts/run_r.sh` launches an R script with the shared path settings.
Configuration files and utility launchers are not analysis steps, so they do not
have execution-order prefixes. See [INPUTS.md](INPUTS.md) for the original QC
cohorts, split memberships, genotype files and prioritized SNP lists required
before starting. Existing external data/output paths keep their original names.

## PRS construction

Run the scripts in `prs_construction/` in the following order. Trait/group steps
must be repeated for each required trait and variant set. The launchers retain
the original PRS-CS settings and seed; inspect their trait/group arguments.

| Step | Script | Purpose |
| --- | --- | --- |
| 01 | `01_format_mdd_gwas.R` | Format MDD GWAS and target BIM inputs. |
| 02 | `02_format_related_trait_gwas.R` | Format related-trait GWAS inputs. |
| 03 | `03_format_related_trait_variant_sets.R` | Prepare related-trait SNP sets with MDD effects. |
| 04 | `04_format_related_trait_effect_sizes.R` | Prepare MDD SNP sets with related-trait effects. |
| 05 | `05_run_mdd_prscs.sh` | Fit the full MDD score. |
| 06 | `06_run_mdd_variant_sets.sh` | Fit MDD scores for each prioritized set. |
| 07 | `07_run_related_trait_prscs.sh` | Fit related-trait scores. |
| 08 | `08_run_related_trait_variant_sets.sh` | Fit related-trait SNP sets with MDD effects. |
| 09 | `09_run_related_trait_effect_sizes.sh` | Fit MDD SNP sets with related-trait effects. |
| 10 | `10_combine_chromosome_weights.R` | Combine chromosome weights after fitting completes. |
| 11 | `11_score_mdd_prs.R` | Calculate the full MDD score. |
| 12 | `12_score_mdd_variant_sets.R` | Calculate prioritized MDD scores. |
| 13 | `13_score_related_trait_prs.R` | Calculate related-trait scores. |
| 14 | `14_score_related_trait_variant_sets.R` | Score related-trait SNP sets with MDD effects. |
| 15 | `15_score_related_trait_effect_sizes.R` | Score MDD SNP sets with related-trait effects. |
| 16 | `16_collect_mdd_scores.R` | Collect MDD predictions. |
| 17 | `17_collect_related_trait_scores.R` | Collect related-trait predictions. |
| 18 | `18_prepare_mdd_evaluation_data.R` | Prepare the MDD evaluation data. |
| 19 | `19_prepare_ancestry_evaluation_data.R` | Prepare ancestry-specific MDD data. |
| 20 | `20_prepare_related_trait_ancestry_data.R` | Prepare related-trait ancestry data. |
| 21 | `21_prepare_related_trait_eur_data.R` | Prepare related-trait EUR data. |
| 22 | `22_evaluate_ancestry_prs.R` | Evaluate standalone PRS in other ancestries. |
| 23 | `23_calculate_allele_frequencies.R` | Optional allele-frequency calculation. |
| 24 | `24_summarize_allele_frequencies.R` | Optional frequency summaries; follows step 23. |

Steps 11-15 and 23 also have matching `.sh` launchers. Example:

```bash
bash scripts/run_r.sh prs_construction/01_format_mdd_gwas.R
bash prs_construction/06_run_mdd_variant_sets.sh cr95set
```

PRS-CS must be installed at `PRSCS_DIR`; set `PRSCS_PYTHON` to its compatible
interpreter. The original stage names under `MDD_PROJECT`, including
`05.PRS_CS_construction/`, identify data locations, not code folders.

Examples of the positional arguments required by the shell launchers:

| Step | Arguments after the script filename |
| --- | --- |
| 05 | None |
| 06 | `cr95set` |
| 07 | `ADHD hapmap_cr95set` |
| 08–09 | `ADHD cr95set` |
| 11 | `1kg` |
| 12 | `1kg cr95set` |
| 13 | `1kg ADHD hapmap_cr95set` |
| 14–15 | `1kg ADHD cr95set` |
| 23 | None |

Substitute the required trait and set names for each run. The R scoring scripts
take the same arguments as their matching shell launchers.

## Elastic-net models

Choose a model folder:

| Folder | Models |
| --- | --- |
| `elastic_net/mdd_only/` | MDD-only combined scores. |
| `elastic_net/multi_trait/` | Scores from MDD and related traits. |
| `elastic_net/cross_trait_variant_sets/` | Related-trait SNP sets with MDD effects. |
| `elastic_net/cross_trait_effect_sizes/` | MDD SNP sets with related-trait effects. |

Each folder follows the same sequence:

1. `01_train_models.R` — train each required model for each split.
2. `02_predict_prs.R` — apply the fitted models.
3. `03_collect_predictions.R` — collect saved predictions.
4. `04_evaluate_models.R` — calculate EUR performance.
5. `05_summarize_coefficients.R` — optional coefficient summaries.
6. `06_predict_other_ancestries.R` — apply fitted models to AFR/SAS.
7. `07_collect_ancestry_predictions.R` — organize AFR/SAS predictions.
8. `08_evaluate_other_ancestries.R` — calculate AFR/SAS performance.

Training and prediction also have matching `.sh` array launchers. Repeat the
required operations for splits 1-10 and the model settings used by the study.
The production training uses 10-fold CV repeated ten times.

```bash
# MDD-only, split 1, covariates plus cr95set/logbf4 scores.
bash scripts/run_r.sh elastic_net/mdd_only/01_train_models.R 1 ALLCovar_NoTWAS
# Related-trait model, split 1, cr95set-only scores with covariates.
bash scripts/run_r.sh elastic_net/multi_trait/01_train_models.R 1 cr95setonly_ALLCovar
```

For multi/cross-trait models, the other main setting is
`cr95setlogbf4_ALLCovar`. Historical `NoCovar` and TWAS alternatives are retained
where used by the original scripts. Some readers load all historical score
columns even when a particular model uses only a subset; retain the input schema.

For SLURM, load your environment, create `logs/` in the relevant script directory
and submit from that directory. The batch launchers do not automatically wait
for earlier steps. Adapt resource requests to your cluster.

## Figures

Plots have their own result inputs and can be run separately after those inputs
exist. Their numbers organize the main comparisons followed by supplementary
summaries; they do not create missing upstream data.

| Step | Script | Content |
| --- | --- | --- |
| 01 | `figures/01_plot_variant_set_overlap.R` | SNP-set overlap plots. |
| 02 | `figures/02_plot_partial_prs_auc.R` | Standalone partial-PRS AUC comparisons. |
| 03 | `figures/03_plot_combined_prs_performance.R` | Combined-model EUR performance. |
| 04 | `figures/04_plot_variant_set_and_effect_comparisons.R` | SNP-set/effect-source comparisons. |
| 05 | `figures/05_plot_ancestry_transfer.R` | AFR/SAS transfer results. |
| 06 | `figures/06_plot_auc_per_variant.R` | Incremental AUC per variant. |
| 07 | `figures/07_plot_partial_prs_risk_metrics.R` | R2 and tail-odds summaries. |
| 08 | `figures/08_plot_logistic_coefficients.R` | Optional historical logistic-coefficient comparisons; requires the original OLD/P-version inputs. |

The plotting inputs `auc_eur_res.csv` and `nagr2_diff_eur_res.csv` must be supplied
from the original analysis. Their generators were not found in the source
collection. See [MANUSCRIPT_MAP.md](MANUSCRIPT_MAP.md) for figure coverage.

## Checks

```bash
bash scripts/run_checks.sh
```

This checks syntax, release-file integrity and path patterns without fitting
models or accessing participant data. It does not rerun the study.
