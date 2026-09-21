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

time_varying_configurations <- list(
  "7-gene full panel" = panel_genes,
  "5-gene reduced" = setdiff(panel_genes, c("TP53", "PIK3CA")),
  "TP53-only" = "TP53",
  "PIK3CA-only" = "PIK3CA"
)

extract_time_varying_term <- function(model, configuration, term_label, main_term, time_term) {
  coefficients <- coef(model)
  variance <- vcov(model)
  model_summary <- summary(model)$coefficients
  main_estimate <- unname(coefficients[main_term])
  time_estimate <- unname(coefficients[time_term])
  main_se <- sqrt(variance[main_term, main_term])
  time_se <- sqrt(variance[time_term, time_term])
  joint_terms <- c(main_term, time_term)
  joint_coefficients <- coefficients[joint_terms]
  joint_variance <- variance[joint_terms, joint_terms, drop = FALSE]
  joint_chisq <- as.numeric(t(joint_coefficients) %*% solve(joint_variance, joint_coefficients))

  data.frame(
    term = term_label,
    configuration = configuration,
    HR_at_1_year = exp(main_estimate),
    HR_at_1_year_lower = exp(main_estimate - qnorm(0.975) * main_se),
    HR_at_1_year_upper = exp(main_estimate + qnorm(0.975) * main_se),
    main_p = unname(model_summary[main_term, "Pr(>|z|)"]),
    time_interaction_ratio = exp(time_estimate),
    time_interaction_lower = exp(time_estimate - qnorm(0.975) * time_se),
    time_interaction_upper = exp(time_estimate + qnorm(0.975) * time_se),
    time_interaction_p = unname(model_summary[time_term, "Pr(>|z|)"]),
    joint_wald_chisq = joint_chisq,
    joint_p = pchisq(joint_chisq, df = 2, lower.tail = FALSE),
    N = model$n,
    Events = model$nevent,
    row.names = NULL,
    check.names = FALSE
  )
}

extract_time_specific_term <- function(model, configuration, term_label, main_term, time_term, times) {
  coefficients <- coef(model)[c(main_term, time_term)]
  variance <- vcov(model)[c(main_term, time_term), c(main_term, time_term), drop = FALSE]
  do.call(rbind, lapply(times, function(time_years) {
    contrast <- c(1, log(time_years))
    log_hr <- sum(contrast * coefficients)
    log_hr_se <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))
    data.frame(
      term = term_label,
      configuration = configuration,
      time_years = time_years,
      HR = exp(log_hr),
      HR_lower = exp(log_hr - qnorm(0.975) * log_hr_se),
      HR_upper = exp(log_hr + qnorm(0.975) * log_hr_se),
      row.names = NULL,
      check.names = FALSE
    )
  }))
}

fit_time_varying_configuration <- function(label, genes) {
  model_data <- analysis_data
  model_data$test_burden <- rowSums(model_data[paste0(genes, "_mutated")])
  model <- coxph(
    Surv(time_years, event) ~
      test_burden + log_tmb + age + mol_subtype + strata(stage) +
      tt(test_burden) + tt(log_tmb),
    data = model_data,
    na.action = na.exclude,
    tt = function(x, t, ...) x * log(t)
  )
  list(
    model = model,
    summary = rbind(
      extract_time_varying_term(model, label, "Panel burden term", "test_burden", "tt(test_burden)"),
      extract_time_varying_term(model, label, "Tumour mutational burden term", "log_tmb", "tt(log_tmb)")
    ),
    time_specific = rbind(
      extract_time_specific_term(model, label, "Panel burden term", "test_burden", "tt(test_burden)", c(1, 2.4, 5)),
      extract_time_specific_term(model, label, "Tumour mutational burden term", "log_tmb", "tt(log_tmb)", c(1, 2.4, 5))
    )
  )
}

time_varying_fits <- lapply(
  names(time_varying_configurations),
  function(label) fit_time_varying_configuration(label, time_varying_configurations[[label]])
)

time_varying_results <- do.call(rbind, lapply(time_varying_fits, function(result) result$summary))
time_specific_results <- do.call(rbind, lapply(time_varying_fits, function(result) result$time_specific))

term_order <- c("Panel burden term", "Tumour mutational burden term")
configuration_order <- names(time_varying_configurations)
time_varying_results <- time_varying_results[order(
  match(time_varying_results$term, term_order),
  match(time_varying_results$configuration, configuration_order)
), ]
time_specific_results <- time_specific_results[order(
  match(time_specific_results$term, term_order),
  match(time_specific_results$configuration, configuration_order),
  time_specific_results$time_years
), ]
rownames(time_varying_results) <- NULL
rownames(time_specific_results) <- NULL

write.csv(
  time_varying_results,
  file.path(tables_dir, "supplementary_table8b_time_varying_cox.csv"),
  row.names = FALSE
)
write.csv(
  time_specific_results,
  file.path(tables_dir, "time_specific_extended_cox_estimates.csv"),
  row.names = FALSE
)

cat("Generated Table 6, Supplementary Table 8B, and time-specific extended Cox estimates.\n")
