# Package: survival. Prerequisite: data/tcga_processed/tcga_brca_panel_survival_input.csv.

required_packages <- c("survival")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required R packages before running: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(survival))

repo_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
input_file <- file.path(repo_root, "data", "tcga_processed", "tcga_brca_panel_survival_input.csv")
tables_dir <- file.path(repo_root, "results", "tables")
figures_dir <- file.path(repo_root, "results", "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

analysis_data <- read.csv(input_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
analysis_data$stage <- factor(analysis_data$stage, levels = c("I", "II", "III", "IV"))
analysis_data$mol_subtype <- factor(
  analysis_data$mol_subtype,
  levels = c("HR+/HER2-", "HR+/HER2+", "HER2+", "Triple Negative")
)
analysis_data$mol_subtype <- addNA(analysis_data$mol_subtype)
levels(analysis_data$mol_subtype)[is.na(levels(analysis_data$mol_subtype))] <- "Unknown"
analysis_data$mol_subtype <- relevel(analysis_data$mol_subtype, ref = "HR+/HER2-")
analysis_data$log_tmb <- log1p(analysis_data$tmb)
panel_genes <- c("TP53", "FLG", "PIK3CA", "CACNA1A", "GOLGB1", "COL14A1", "MUC5B")

if (nrow(analysis_data) != 1083 || sum(analysis_data$event) != 151) {
  stop("The TCGA OS input does not match the dissertation cohort (expected N=1083, events=151).")
}
