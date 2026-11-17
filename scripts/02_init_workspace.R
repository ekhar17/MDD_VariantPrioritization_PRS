# Create output directories; no input files, participants or results are created.
repo <- Sys.getenv("MDD_CODE_ROOT", getwd())
source(file.path(repo, "config", "paths.R"))
project <- Sys.getenv("MDD_PROJECT")
directories <- c("01.TrainTestSplit/data/splits", "05.PRS_CS_construction/data/bim",
  "05.PRS_CS_construction/data/clean_sumstat", "05.PRS_CS_construction/data/prs-cs_output_redo",
  "05.PRS_CS_construction/data/prs-cs_prediction_redo", "05.PRS_CS_construction/data/clean_output_redo/ldref_1kg",
  "08.Generate_Figures/new_fig", paste0("08.Generate_Figures/fig_v", 1:4),
  "reports/auc")
for (group in c("single_trait", "multi_trait", "muti_trait_MDDcoef_alttraitset", "muti_trait_MDDset_alttraitcoef")) {
  directories <- c(directories, file.path("06.Elastic_Net", group,
    c("results", "results_altancestry", "slurm_output")))
}
for (d in directories) dir.create(file.path(project, d), recursive=TRUE, showWarnings=FALSE)
dir.create(Sys.getenv("MDD_REPORTS"), recursive=TRUE, showWarnings=FALSE)
cat("Created workspace directories under:", project, "\n")
