# R entry points also work when MDD_CODE_ROOT is exported without the shell loader.
.mdd_set_default <- function(key, value) {
  if (!nzchar(Sys.getenv(key))) do.call(Sys.setenv, setNames(list(value), key))
  invisible(Sys.getenv(key))
}
.mdd_set_default("MDD_CODE_ROOT", getwd())
.mdd_set_default("MDD_PROJECT", file.path(Sys.getenv("MDD_CODE_ROOT"), "workspace"))
.mdd_set_default("MDD_SOURCE_DATA", file.path(Sys.getenv("MDD_PROJECT"), "source_data"))
.mdd_set_default("MDD_REFERENCES", file.path(Sys.getenv("MDD_PROJECT"), "references"))
.mdd_set_default("MDD_TOOLS", file.path(Sys.getenv("MDD_PROJECT"), "tools"))
.mdd_set_default("MDD_REPORTS", file.path(Sys.getenv("MDD_PROJECT"), "reports"))
.mdd_set_default("MDD_HAPMAP", file.path(Sys.getenv("MDD_REFERENCES"), "snpinfo_1kg_hm3"))
.mdd_set_default("MDD_LD_REF", file.path(Sys.getenv("MDD_REFERENCES"), "ldblk_1kg_eur"))
.mdd_set_default("MDD_BACKING", file.path(Sys.getenv("MDD_PROJECT"), "genotypes", "imputed_bgen_hapmap3_allind_mult_1kg.rds"))
.mdd_set_default("MDD_SAMPLE", file.path(Sys.getenv("MDD_PROJECT"), "genotypes", "ukb.sample"))
.mdd_set_default("PRSCS_DIR", file.path(Sys.getenv("MDD_TOOLS"), "PRScs"))
.mdd_set_default("PRSCS_PYTHON", "python3")
mdd_path <- function(key, suffix = "") {
  value <- Sys.getenv(key)
  if (!nzchar(value)) stop("Set the path variable ", key, " before running this script.")
  # Keep trailing separators: historical paste() calls concatenate file prefixes.
  paste0(sub("/+$", "", value), suffix)
}
