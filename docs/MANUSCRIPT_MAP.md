# Core figure and result map

This selected upload contains the original PRS construction, elastic-net and
figure code. Its scope does not include additional benchmark/sensitivity modules
or final image-layout assembly. The plotting scripts include historical layouts;
their presence does not establish pixel-identical regeneration of every panel.

| Result | Included code |
| --- | --- |
| Figure 2A and Figure S5 SNP-set overlaps | `figures/01_plot_variant_set_overlap.R` |
| Figure 3 partial-PRS AUC comparisons | `figures/02_plot_partial_prs_auc.R` |
| Figure 4 combined EUR models | `elastic_net/*/04_evaluate_models.R`, `figures/03_plot_combined_prs_performance.R` |
| Figure 5 SNP-set/effect-source comparisons | `figures/04_plot_variant_set_and_effect_comparisons.R` |
| Figure 6 ancestry transfer | `elastic_net/*/08_evaluate_other_ancestries.R`, `figures/05_plot_ancestry_transfer.R` |
| Figure S2 AUC per variant | `figures/06_plot_auc_per_variant.R` |
| Figure S4 risk summaries | `figures/07_plot_partial_prs_risk_metrics.R` |
| Historical logistic-coefficient comparisons | `figures/08_plot_logistic_coefficients.R`; extra OLD/P-version inputs required |

The diagram source, initial statistical panel-generation scripts for Figure 2B/C,
original stratified-LDSC enrichment code and some intermediate-table generators
were not present. No substitute statistical results are generated for them.

AFR/SAS model evaluations reuse fixed ancestry cohorts across EUR-trained model
fits. Ten fits do not represent ten independent population samples.
