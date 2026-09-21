# Package: base R. Prerequisite: source scripts 01-04.

check_round <- function(observed, expected, digits, label) {
  if (length(observed) != 1 || is.na(observed) || round(observed, digits) != round(expected, digits)) {
    stop(sprintf("%s: observed %s, expected %s at %d decimal places", label, observed, expected, digits))
  }
}

check_integer <- function(observed, expected, label) check_round(observed, expected, 0, label)

check_p <- function(observed, expected, label) {
  if (is.na(expected)) {
    if (is.na(observed) || observed >= 0.001) stop(label, ": expected p<0.001, observed ", observed)
  } else {
    check_round(observed, expected, 3, label)
  }
}

validate_cox <- function(observed, expected, label) {
  if (!identical(as.character(observed$variable), as.character(expected$variable))) {
    stop(label, " coefficient order or labels differ from the expected table")
  }
  for (index in seq_len(nrow(expected))) {
    variable <- expected$variable[index]
    check_round(observed$HR[index], expected$HR[index], 3, paste(label, variable, "HR"))
    check_round(observed$HR_lower[index], expected$lower[index], 3, paste(label, variable, "CI lower"))
    check_round(observed$HR_upper[index], expected$upper[index], 3, paste(label, variable, "CI upper"))
    check_p(observed$p[index], expected$p[index], paste(label, variable, "p"))
  }
}

expected_m1 <- data.frame(
  variable = c("panel_burden", "age", "stageII", "stageIII", "stageIV", "mol_subtypeHR+/HER2+", "mol_subtypeHER2+", "mol_subtypeTriple Negative", "mol_subtypeUnknown"),
  HR = c(1.164, 1.037, 1.705, 3.658, 12.246, 1.332, 2.207, 1.966, 1.023),
  lower = c(0.917, 1.023, 0.979, 2.039, 5.933, 0.733, 0.976, 1.117, 0.670),
  upper = c(1.479, 1.051, 2.969, 6.563, 25.277, 2.421, 4.987, 3.462, 1.564),
  p = c(0.212, NA, 0.059, NA, NA, 0.347, 0.057, 0.019, 0.916)
)
expected_m2 <- data.frame(
  variable = c("log_tmb", "age", "stageII", "stageIII", "stageIV", "mol_subtypeHR+/HER2+", "mol_subtypeHER2+", "mol_subtypeTriple Negative", "mol_subtypeUnknown"),
  HR = c(1.133, 1.037, 1.714, 3.937, 10.518, 1.370, 1.938, 2.050, 1.014),
  lower = c(0.777, 1.023, 0.967, 2.160, 4.830, 0.735, 0.807, 1.136, 0.649),
  upper = c(1.654, 1.051, 3.038, 7.176, 22.905, 2.557, 4.653, 3.699, 1.585),
  p = c(0.516, NA, 0.065, NA, NA, 0.322, 0.139, 0.017, 0.951)
)
expected_m3 <- data.frame(
  variable = c("panel_burden", "log_tmb", "age", "stageII", "stageIII", "stageIV", "mol_subtypeHR+/HER2+", "mol_subtypeHER2+", "mol_subtypeTriple Negative", "mol_subtypeUnknown"),
  HR = c(1.200, 1.033, 1.038, 1.753, 4.107, 10.842, 1.368, 1.878, 2.059, 1.021),
  lower = c(0.922, 0.689, 1.024, 0.988, 2.243, 4.964, 0.733, 0.781, 1.141, 0.653),
  upper = c(1.563, 1.547, 1.053, 3.111, 7.520, 23.677, 2.556, 4.517, 3.715, 1.596),
  p = c(0.175, 0.876, NA, 0.055, NA, NA, 0.325, 0.159, 0.016, 0.927)
)
expected_m4 <- data.frame(
  variable = c("panel_binary", "log_tmb", "age", "stageII", "stageIII", "stageIV", "mol_subtypeHR+/HER2+", "mol_subtypeHER2+", "mol_subtypeTriple Negative", "mol_subtypeUnknown"),
  HR = c(1.141, 1.094, 1.038, 1.736, 4.012, 10.488, 1.368, 1.913, 2.020, 1.017),
  lower = c(0.783, 0.737, 1.023, 0.978, 2.195, 4.813, 0.732, 0.796, 1.117, 0.651),
  upper = c(1.662, 1.624, 1.052, 3.080, 7.333, 22.855, 2.554, 4.598, 3.651, 1.590),
  p = c(0.493, 0.656, NA, 0.060, NA, NA, 0.326, 0.147, 0.020, 0.940)
)

validate_cox(res1, expected_m1, "Supplementary Table 4")
validate_cox(res2, expected_m2, "Supplementary Table 5")
validate_cox(res3, expected_m3, "Table 5")
validate_cox(res4, expected_m4, "Supplementary Table 6")

expected_comparison <- data.frame(
  N = c(1057, 941, 941, 941), Events = c(140, 129, 129, 129), EPV = c(15.6, 14.3, 12.9, 12.9),
  C_index = c(0.769, 0.767, 0.766, 0.768),
  Panel_HR = c(1.16, NA, 1.20, 1.14), Panel_p = c(0.212, NA, 0.175, 0.493),
  TMB_HR = c(NA, 1.13, 1.03, 1.09), TMB_p = c(NA, 0.516, 0.876, 0.656),
  PH_global_p = c(0.027, 0.011, 0.007, 0.006)
)
for (index in seq_len(nrow(expected_comparison))) {
  check_integer(model_comparison$N[index], expected_comparison$N[index], paste("Supplementary Table 3 model", index, "N"))
  check_integer(model_comparison$Events[index], expected_comparison$Events[index], paste("Supplementary Table 3 model", index, "events"))
  check_round(model_comparison$EPV[index], expected_comparison$EPV[index], 1, paste("Supplementary Table 3 model", index, "EPV"))
  check_round(model_comparison$C_index[index], expected_comparison$C_index[index], 3, paste("Supplementary Table 3 model", index, "C-index"))
  for (column in c("Panel_HR", "Panel_p", "TMB_HR", "TMB_p")) {
    expected_value <- expected_comparison[[column]][index]
    observed_value <- model_comparison[[column]][index]
    if (is.na(expected_value)) {
      if (!is.na(observed_value)) stop("Supplementary Table 3 model ", index, " ", column, " should be NA")
    } else {
      digits <- if (grepl("_HR$", column)) 2 else 3
      check_round(observed_value, expected_value, digits, paste("Supplementary Table 3 model", index, column))
    }
  }
  check_round(model_comparison$PH_global_p[index], expected_comparison$PH_global_p[index], 3, paste("Supplementary Table 3 model", index, "PH p"))
}

expected_sensitivity <- data.frame(
  HR = c(1.200, 1.010, 1.374, 1.081), lower = c(0.922, 0.735, 0.988, 0.648),
  upper = c(1.563, 1.386, 1.910, 1.805), p = c(0.175, 0.953, 0.059, 0.765),
  C_index = c(0.766, 0.767, 0.768, 0.766), PH = c(0.007, 0.010, 0.008, 0.011)
)
for (index in seq_len(nrow(expected_sensitivity))) {
  check_round(sensitivity_results$HR[index], expected_sensitivity$HR[index], 3, paste("Table 6 row", index, "HR"))
  check_round(sensitivity_results$HR_lower[index], expected_sensitivity$lower[index], 3, paste("Table 6 row", index, "CI lower"))
  check_round(sensitivity_results$HR_upper[index], expected_sensitivity$upper[index], 3, paste("Table 6 row", index, "CI upper"))
  check_round(sensitivity_results$p[index], expected_sensitivity$p[index], 3, paste("Table 6 row", index, "p"))
  check_round(sensitivity_results$C_index[index], expected_sensitivity$C_index[index], 3, paste("Table 6 row", index, "C-index"))
  check_integer(sensitivity_results$N[index], 941, paste("Table 6 row", index, "N"))
  check_integer(sensitivity_results$Events[index], 129, paste("Table 6 row", index, "events"))
  check_round(sensitivity_results$PH_global_p[index], expected_sensitivity$PH[index], 3, paste("Table 6 row", index, "PH p"))
}

expected_genes <- data.frame(
  gene = c("TP53", "FLG", "PIK3CA", "CACNA1A", "GOLGB1", "COL14A1", "MUC5B"),
  n_mut = c(296, 43, 313, 17, 13, 16, 35), n_wt = c(787, 1040, 770, 1066, 1070, 1067, 1048),
  uni_HR = c(1.26, 0.92, 0.97, 1.23, 1.41, 0.82, 0.65), uni_p = c(0.187, 0.820, 0.853, 0.728, 0.627, 0.775, 0.465),
  multi_HR = c(1.49, 1.18, 0.94, 2.24, 2.57, 0.85, 0.52), multi_lower = c(1.01, 0.55, 0.64, 0.69, 0.62, 0.21, 0.13),
  multi_upper = c(2.18, 2.54, 1.39, 7.22, 10.59, 3.44, 2.11), multi_p = c(0.043, 0.674, 0.774, 0.178, 0.192, 0.816, 0.358),
  tmb_HR = c(1.54, 1.15, 0.97, 2.21, 2.41, 0.75, 0.51), tmb_lower = c(1.02, 0.52, 0.65, 0.68, 0.57, 0.18, 0.12),
  tmb_upper = c(2.32, 2.50, 1.44, 7.21, 10.26, 3.20, 2.08), tmb_p = c(0.038, 0.733, 0.881, 0.189, 0.233, 0.699, 0.348),
  multi_cindex = c(0.775, 0.768, 0.769, 0.768, 0.767, 0.768, 0.770)
)
if (!identical(as.character(gene_results$gene), as.character(expected_genes$gene))) stop("Supplementary Table 7 gene order differs")
for (index in seq_len(nrow(expected_genes))) {
  gene <- expected_genes$gene[index]
  check_integer(gene_results$n_mut[index], expected_genes$n_mut[index], paste("Supplementary Table 7", gene, "n_mut"))
  check_integer(gene_results$n_wt[index], expected_genes$n_wt[index], paste("Supplementary Table 7", gene, "n_wt"))
  for (column in c("uni_HR", "multi_HR", "multi_lower", "multi_upper", "tmb_HR", "tmb_lower", "tmb_upper")) {
    check_round(gene_results[[column]][index], expected_genes[[column]][index], 2, paste("Supplementary Table 7", gene, column))
  }
  for (column in c("uni_p", "multi_p", "tmb_p")) {
    check_round(gene_results[[column]][index], expected_genes[[column]][index], 3, paste("Supplementary Table 7", gene, column))
  }
  check_round(gene_results$multi_cindex[index], expected_genes$multi_cindex[index], 3, paste("Supplementary Table 7", gene, "C-index"))
  check_integer(gene_results$multi_n[index], 1057, paste("Supplementary Table 7", gene, "multivariable N"))
  check_integer(gene_results$multi_events[index], 140, paste("Supplementary Table 7", gene, "multivariable events"))
  check_integer(gene_results$tmb_n[index], 941, paste("Supplementary Table 7", gene, "TMB-adjusted N"))
  check_integer(gene_results$tmb_events[index], 129, paste("Supplementary Table 7", gene, "TMB-adjusted events"))
}

expected_ph <- data.frame(
  model = c("Model 1", "Model 1", "Model 1", "Model 1", "Model 1", "Model 3", "Model 3", "Model 3", "Model 3", "Model 3", "Model 3"),
  variable = c("panel_burden", "age", "stage", "mol_subtype", "GLOBAL", "panel_burden", "log_tmb", "age", "stage", "mol_subtype", "GLOBAL"),
  chisq = c(2.750, 0.039, 10.958, 3.340, 18.830, 4.112, 6.242, 0.105, 9.913, 3.302, 24.221),
  p = c(0.097, 0.843, 0.012, 0.503, 0.027, 0.043, 0.012, 0.746, 0.019, 0.509, 0.007)
)
for (index in seq_len(nrow(expected_ph))) {
  row <- ph_results[ph_results$model == expected_ph$model[index] & ph_results$variable == expected_ph$variable[index], ]
  if (nrow(row) != 1) stop("Supplementary Table 8A missing or duplicate row: ", expected_ph$model[index], " ", expected_ph$variable[index])
  check_round(row$chisq, expected_ph$chisq[index], 3, paste("Supplementary Table 8A", expected_ph$model[index], expected_ph$variable[index], "chisq"))
  check_round(row$p, expected_ph$p[index], 3, paste("Supplementary Table 8A", expected_ph$model[index], expected_ph$variable[index], "p"))
}
if (!identical(
  names(ph_table8a),
  c("Variable", "Model 1 chisq", "Model 1 p", "Model 3 chisq", "Model 3 p")
) || !identical(
  as.character(ph_table8a$Variable),
  c("Panel burden", "log1p(TMB)", "Age", "Stage", "Molecular subtype", "GLOBAL")
)) {
  stop("Supplementary Table 8A layout or row labels differ from the dissertation")
}

expected_time_varying <- data.frame(
  term = rep(c("Panel burden term", "Tumour mutational burden term"), each = 4),
  configuration = rep(c("7-gene full panel", "5-gene reduced", "TP53-only", "PIK3CA-only"), 2),
  HR_at_1_year = c(1.266, 0.911, 2.015, 0.954, 1.467, 1.792, 1.702, 1.736),
  HR_at_1_year_lower = c(0.836, 0.466, 1.052, 0.483, 0.912, 1.067, 1.125, 1.164),
  HR_at_1_year_upper = c(1.915, 1.779, 3.859, 1.885, 2.361, 3.007, 2.576, 2.589),
  main_p = c(0.265, 0.784, 0.035, 0.892, 0.114, 0.027, 0.012, 0.007),
  time_interaction_ratio = c(0.905, 0.962, 0.818, 0.988, 0.606, 0.558, 0.526, 0.560),
  time_interaction_lower = c(0.675, 0.584, 0.504, 0.616, 0.412, 0.369, 0.361, 0.395),
  time_interaction_upper = c(1.213, 1.587, 1.326, 1.585, 0.890, 0.842, 0.766, 0.795),
  time_interaction_p = c(0.505, 0.880, 0.414, 0.961, 0.011, 0.005, 0.001, 0.001),
  joint_wald_chisq = c(1.34, 0.26, 6.02, 0.090, 6.56, 8.09, 12.64, 12.36),
  joint_p = c(0.512, 0.879, 0.049, 0.955, 0.038, 0.018, 0.002, 0.002)
)
if (!identical(as.character(time_varying_results$term), expected_time_varying$term) ||
    !identical(as.character(time_varying_results$configuration), expected_time_varying$configuration)) {
  stop("Supplementary Table 8B term or configuration order differs")
}
for (index in seq_len(nrow(expected_time_varying))) {
  label <- paste("Supplementary Table 8B", expected_time_varying$term[index], expected_time_varying$configuration[index])
  for (column in c(
    "HR_at_1_year", "HR_at_1_year_lower", "HR_at_1_year_upper", "main_p",
    "time_interaction_ratio", "time_interaction_lower", "time_interaction_upper",
    "time_interaction_p", "joint_p"
  )) {
    check_round(time_varying_results[[column]][index], expected_time_varying[[column]][index], 3, paste(label, column))
  }
  check_round(time_varying_results$joint_wald_chisq[index], expected_time_varying$joint_wald_chisq[index], 2, paste(label, "joint Wald chisq"))
  check_integer(time_varying_results$N[index], 941, paste(label, "N"))
  check_integer(time_varying_results$Events[index], 129, paste(label, "events"))
}

reported_time_specific <- time_specific_results[
  time_specific_results$term == "Panel burden term" &
    time_specific_results$configuration == "7-gene full panel",
]
expected_time_specific <- data.frame(
  time_years = c(1, 2.4, 5),
  HR = c(1.27, 1.16, 1.08),
  HR_lower = c(0.84, 0.89, 0.80),
  HR_upper = c(1.92, 1.52, 1.46)
)
if (nrow(reported_time_specific) != nrow(expected_time_specific) ||
    any(reported_time_specific$time_years != expected_time_specific$time_years)) {
  stop("Time-specific full-panel estimates do not contain the reported 1-, 2.4-, and 5-year rows")
}
for (index in seq_len(nrow(expected_time_specific))) {
  for (column in c("HR", "HR_lower", "HR_upper")) {
    check_round(reported_time_specific[[column]][index], expected_time_specific[[column]][index], 2, paste("Section 3.8 time-specific panel", expected_time_specific$time_years[index], "years", column))
  }
}

expected_km_groups <- c(
  "0 (500), 1 (455), >=2 (128)",
  "Wild-type (500), Mutated (583)",
  "Wild-type (787), Mutated (296)",
  "Low (333), Mid (313), High (318)"
)
expected_km_p <- c(0.959, 0.826, 0.186, 0.281)
if (!identical(as.character(km_results$groups), expected_km_groups)) stop("Supplementary Table 9 group counts or labels differ")
for (index in seq_along(expected_km_p)) check_round(km_results$log_rank_p[index], expected_km_p[index], 3, paste("Supplementary Table 9 row", index, "p"))

validate_risk_table <- function(filename, groups, expected_counts, label) {
  observed <- read.csv(file.path(tables_dir, filename), stringsAsFactors = FALSE)
  times <- seq(0, 25, 5)
  expected_groups <- rep(groups, each = length(times))
  expected_times <- rep(times, times = length(groups))
  expected_n <- unlist(expected_counts, use.names = FALSE)
  if (!identical(as.character(observed$group), expected_groups) ||
      any(observed$time_years != expected_times) ||
      any(observed$n_at_risk != expected_n)) {
    stop(label, " numbers at risk differ from the dissertation")
  }
}

validate_risk_table(
  "table7_numbers_at_risk_panel_burden.csv",
  c("0", "1", ">=2"),
  list(c(500, 113, 19, 4, 1, 0), c(455, 103, 21, 6, 4, 0), c(128, 37, 2, 2, 1, 0)),
  "Table 7"
)
validate_risk_table(
  "table8_numbers_at_risk_panel_status.csv",
  c("Wild-type", "Mutated"),
  list(c(500, 113, 19, 4, 1, 0), c(583, 140, 23, 8, 5, 0)),
  "Table 8"
)
validate_risk_table(
  "table9_numbers_at_risk_tp53.csv",
  c("Wild-type", "Mutated"),
  list(c(787, 177, 31, 7, 3, 0), c(296, 76, 11, 5, 3, 0)),
  "Table 9"
)
validate_risk_table(
  "table10_numbers_at_risk_tmb_tertiles.csv",
  c("Low", "Mid", "High"),
  list(c(333, 96, 11, 5, 1, 0), c(313, 71, 15, 4, 3, 0), c(318, 67, 12, 3, 2, 0)),
  "Table 10"
)

required_outputs <- c(
  "table5_cox_model3.csv", "table6_sensitivity_configurations.csv",
  "table7_numbers_at_risk_panel_burden.csv", "table8_numbers_at_risk_panel_status.csv",
  "table9_numbers_at_risk_tp53.csv", "table10_numbers_at_risk_tmb_tertiles.csv",
  "supplementary_table3_model_comparison.csv", "supplementary_table4_cox_model1.csv",
  "supplementary_table5_cox_model2.csv", "supplementary_table6_cox_model4.csv",
  "supplementary_table7_individual_genes.csv", "supplementary_table8a_ph_tests.csv",
  "supplementary_table8b_time_varying_cox.csv", "supplementary_table9_km_summary.csv",
  "time_specific_extended_cox_estimates.csv"
)
missing_tables <- required_outputs[!file.exists(file.path(tables_dir, required_outputs))]
if (length(missing_tables) > 0) stop("Missing R table outputs: ", paste(missing_tables, collapse = ", "))

figure_stems <- c(
  "fig4_forest_model3", "fig5_individual_gene_forest", "fig6_cindex_comparison",
  "fig7_km_panel_burden_groups", "fig8_km_panel_binary", "fig9_km_tp53", "fig10_km_tmb_tertiles"
)
required_figures <- unlist(lapply(figure_stems, function(stem) paste0(stem, c(".png", ".svg"))))
figure_info <- file.info(file.path(figures_dir, required_figures))
missing_figures <- required_figures[is.na(figure_info$size) | figure_info$size == 0]
if (length(missing_figures) > 0) stop("Missing R figure outputs: ", paste(missing_figures, collapse = ", "))

cat("All displayed R manuscript values passed validation.\n")
