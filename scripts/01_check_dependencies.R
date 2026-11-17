# Report the current environment without installing or changing packages.
packages <- c("data.table", "tidyverse", "ggplot2", "dplyr", "tidyr", "purrr",
              "fmsb", "pscl", "bigsnpr", "bigstatsr", "caret", "glmnet",
              "glmnetUtils", "ggpattern", "colorspace", "UpSetR", "ComplexUpset",
              "pROC")
rows <- lapply(packages, function(p) {
  present <- requireNamespace(p, quietly=TRUE)
  data.frame(package=p, installed=present,
             version=if (present) as.character(packageVersion(p)) else NA_character_)
})
cat(R.version.string, "\n\n")
print(do.call(rbind, rows), row.names=FALSE)
cat("\nOnly install dependencies needed for the analysis module you run.\n")
