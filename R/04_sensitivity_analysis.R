# Package: survival. Prerequisite: source 00_setup.R.

if (!exists("analysis_data")) source("00_setup.R")

panel_configurations <- list(
  "Full seven-gene panel" = panel_genes,
  "TP53-excluded six-gene panel" = setdiff(panel_genes, "TP53"),
  "PIK3CA-excluded six-gene panel" = setdiff(panel_genes, "PIK3CA"),
  "TP53/PIK3CA-excluded five-gene panel" = setdiff(panel_genes, c("TP53", "PIK3CA"))
)

fit_configuration <- function(label, genes) {
  model_data <- analysis_data
  model_data$test_burden <- rowSums(model_data[paste0(genes, "_mutated")])
  model <- coxph(
    Surv(time_years, event) ~ test_burden + log_tmb + age + stage + mol_subtype,
    data = model_data, na.action = na.exclude
  )
  model_summary <- summary(model)
  ph <- cox.zph(model)
  data.frame(
    configuration = label,
    included_genes = paste(genes, collapse = "; "),
    n_genes = length(genes),
    HR = model_summary$coefficients["test_burden", "exp(coef)"],
    HR_lower = model_summary$conf.int["test_burden", "lower .95"],
    HR_upper = model_summary$conf.int["test_burden", "upper .95"],
    p = model_summary$coefficients["test_burden", "Pr(>|z|)"],
    C_index = unname(model_summary$concordance[1]),
    N = model$n,
    Events = model$nevent,
    EPV = model$nevent / length(coef(model)),
    PH_global_p = unname(ph$table[nrow(ph$table), "p"]),
    row.names = NULL
  )
}

sensitivity_results <- do.call(
  rbind,
  lapply(names(panel_configurations), function(label) fit_configuration(label, panel_configurations[[label]]))
)
write.csv(sensitivity_results, file.path(tables_dir, "table6_sensitivity_configurations.csv"), row.names = FALSE)
cat("Generated Table 6.\n")
