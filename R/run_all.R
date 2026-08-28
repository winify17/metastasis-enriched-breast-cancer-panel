# Package: survival. Prerequisite: cleaned TCGA survival input.

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", arguments, value = TRUE)
if (length(file_argument) == 1) {
  script_directory <- dirname(normalizePath(sub("^--file=", "", file_argument)))
  setwd(script_directory)
}

for (script in c(
  "00_setup.R",
  "01_cox_models.R",
  "02_individual_genes.R",
  "03_km_curves.R",
  "04_sensitivity_analysis.R",
  "05_validate_expected_results.R"
)) {
  cat("Running", script, "\n")
  source(script, echo = FALSE)
}

cat("R survival pipeline completed.\n")
