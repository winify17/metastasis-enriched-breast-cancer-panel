# Package: survival. Prerequisite: source 00_setup.R.

if (!exists("analysis_data")) source("00_setup.R")

extract_cox <- function(model, model_name) {
  model_summary <- summary(model)
  coefficients <- model_summary$coefficients
  confidence <- model_summary$conf.int
  ph <- cox.zph(model)
  data.frame(
    model = model_name,
    variable = rownames(coefficients),
    HR = coefficients[, "exp(coef)"],
    HR_lower = confidence[, "lower .95"],
    HR_upper = confidence[, "upper .95"],
    p = coefficients[, "Pr(>|z|)"],
    cindex = unname(model_summary$concordance[1]),
    n = model$n,
    events = model$nevent,
    epv = model$nevent / length(coef(model)),
    ph_global_p = unname(ph$table[nrow(ph$table), "p"]),
    row.names = NULL,
    check.names = FALSE
  )
}

m1 <- coxph(Surv(time_years, event) ~ panel_burden + age + stage + mol_subtype,
            data = analysis_data, na.action = na.exclude)
m2 <- coxph(Surv(time_years, event) ~ log_tmb + age + stage + mol_subtype,
            data = analysis_data, na.action = na.exclude)
m3 <- coxph(Surv(time_years, event) ~ panel_burden + log_tmb + age + stage + mol_subtype,
            data = analysis_data, na.action = na.exclude)
m4 <- coxph(Surv(time_years, event) ~ panel_binary + log_tmb + age + stage + mol_subtype,
            data = analysis_data, na.action = na.exclude)

res1 <- extract_cox(m1, "M1: Panel + Clinical")
res2 <- extract_cox(m2, "M2: TMB + Clinical")
res3 <- extract_cox(m3, "M3: Panel + TMB + Clinical")
res4 <- extract_cox(m4, "M4: Binary Panel + TMB + Clinical")

model_comparison <- do.call(rbind, lapply(list(res1, res2, res3, res4), function(result) {
  data.frame(
    Model = result$model[1], N = result$n[1], Events = result$events[1],
    EPV = result$epv[1], C_index = result$cindex[1],
    Panel_HR = if (any(result$variable == "panel_burden")) result$HR[result$variable == "panel_burden"] else if (any(result$variable == "panel_binary")) result$HR[result$variable == "panel_binary"] else NA,
    Panel_p = if (any(result$variable == "panel_burden")) result$p[result$variable == "panel_burden"] else if (any(result$variable == "panel_binary")) result$p[result$variable == "panel_binary"] else NA,
    TMB_HR = if (any(result$variable == "log_tmb")) result$HR[result$variable == "log_tmb"] else NA,
    TMB_p = if (any(result$variable == "log_tmb")) result$p[result$variable == "log_tmb"] else NA,
    PH_global_p = result$ph_global_p[1]
  )
}))

write.csv(res3, file.path(tables_dir, "table5_cox_model3.csv"), row.names = FALSE)
write.csv(model_comparison, file.path(tables_dir, "supplementary_table3_model_comparison.csv"), row.names = FALSE)
write.csv(res1, file.path(tables_dir, "supplementary_table4_cox_model1.csv"), row.names = FALSE)
write.csv(res2, file.path(tables_dir, "supplementary_table5_cox_model2.csv"), row.names = FALSE)
write.csv(res4, file.path(tables_dir, "supplementary_table6_cox_model4.csv"), row.names = FALSE)

ph_rows <- function(model, model_name) {
  values <- as.data.frame(cox.zph(model)$table)
  values$variable <- rownames(values)
  values$model <- model_name
  rownames(values) <- NULL
  values[, c("model", "variable", "chisq", "p")]
}
ph_results <- rbind(ph_rows(m1, "Model 1"), ph_rows(m3, "Model 3"))
write.csv(ph_results, file.path(tables_dir, "supplementary_table8_ph_tests.csv"), row.names = FALSE)

coefficient_labels <- c(
  panel_burden = "Panel burden (per gene)", log_tmb = "ln(TMB + 1)", age = "Age (per year)",
  stageII = "Stage II vs I", stageIII = "Stage III vs I", stageIV = "Stage IV vs I",
  `mol_subtypeHR+/HER2+` = "HR+/HER2+ vs HR+/HER2-", `mol_subtypeHER2+` = "HER2+ vs HR+/HER2-",
  `mol_subtypeTriple Negative` = "Triple negative vs HR+/HER2-", `mol_subtypeUnknown` = "Unknown vs HR+/HER2-"
)
forest_data <- res3
forest_data$label <- unname(coefficient_labels[forest_data$variable])
draw_forest <- function() {
  y <- seq_len(nrow(forest_data))
  colours <- ifelse(forest_data$p < 0.05, "#FF9000", "#1683E6")
  p_text <- ifelse(forest_data$p < 0.001, "p<0.001", sprintf("p=%.3f", forest_data$p))
  stars <- ifelse(forest_data$p < 0.001, "***", ifelse(forest_data$p < 0.01, "**", ifelse(forest_data$p < 0.05, "*", "")))
  annotations <- sprintf("%.2f  (%.2f-%.2f)    %s %s", forest_data$HR, forest_data$HR_lower, forest_data$HR_upper, p_text, stars)
  par(mar = c(5, 14, 4, 4))
  plot(forest_data$HR, y, xlim = c(0, 46), pch = 19, col = colours,
       yaxt = "n", ylab = "", xlab = "Hazard ratio (95% CI)",
       main = "Model 3: panel, TMB, and clinical covariates")
  axis(2, at = y, labels = forest_data$label, las = 1)
  segments(forest_data$HR_lower, y, forest_data$HR_upper, y, col = colours, lwd = 2)
  abline(v = 1, lty = 2, col = "grey50")
  text(27, y, annotations, adj = 0, family = "mono", cex = 0.78)
  mtext("TCGA-BRCA overall survival (N=941, 129 events)", side = 3, line = 0.4)
}
png(file.path(figures_dir, "fig4_forest_model3.png"), width = 3300, height = 1800, res = 300)
draw_forest(); dev.off()
svg(file.path(figures_dir, "fig4_forest_model3.svg"), width = 11, height = 6)
draw_forest(); dev.off()

draw_cindex <- function() {
  labels <- c("M1: Panel\n+ Clinical", "M2: TMB\n+ Clinical", "M3: Panel + TMB\n+ Clinical", "M4: Binary panel\n+ TMB + Clinical")
  colours <- c("#1683E6", "#8CAF43", "#FF9F1C", "#EF8BDB")
  par(mar = c(7, 5, 4, 2))
  positions <- barplot(model_comparison$C_index, names.arg = labels, col = colours,
                       ylim = c(0.70, 0.785), ylab = "C-index", las = 1,
                       main = "Model discrimination comparison")
  text(positions, model_comparison$C_index + 0.003, sprintf("%.3f", model_comparison$C_index), font = 2)
  mtext("In-sample Harrell concordance", side = 3, line = 0.4)
}
png(file.path(figures_dir, "fig6_cindex_comparison.png"), width = 2700, height = 1800, res = 300)
draw_cindex(); dev.off()
svg(file.path(figures_dir, "fig6_cindex_comparison.svg"), width = 9, height = 6)
draw_cindex(); dev.off()

cat("Generated Table 5, Supplementary Tables 3-6 and 8, and Figures 4 and 6.\n")
