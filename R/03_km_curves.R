# Package: survival. Prerequisite: source 00_setup.R.

if (!exists("analysis_data")) source("00_setup.R")

log_rank_p <- function(formula, data) {
  test <- survdiff(formula, data = data)
  1 - pchisq(test$chisq, length(test$n) - 1)
}

analysis_data$burden_group <- factor(
  ifelse(analysis_data$panel_burden == 0, "0", ifelse(analysis_data$panel_burden == 1, "1", ">=2")),
  levels = c("0", "1", ">=2")
)
analysis_data$panel_status <- factor(analysis_data$panel_binary, levels = c(0, 1), labels = c("Wild-type", "Mutated"))
analysis_data$tp53_status <- factor(analysis_data$TP53_mutated, levels = c(0, 1), labels = c("Wild-type", "Mutated"))
analysis_data$tmb_tertile <- cut(
  analysis_data$tmb,
  breaks = quantile(analysis_data$tmb, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE),
  labels = c("Low", "Mid", "High"), include.lowest = TRUE
)

km_specs <- list(
  list(variable = "burden_group", title = "TCGA-BRCA OS by panel burden", legend = "Panel burden", colours = c("#1683E6", "#FF9000", "#EF8BDB"), stem = "fig7_km_panel_burden_groups"),
  list(variable = "panel_status", title = "TCGA-BRCA OS by panel mutation status", legend = "Panel status", colours = c("#1683E6", "#FF9000"), stem = "fig8_km_panel_binary"),
  list(variable = "tp53_status", title = "TCGA-BRCA OS by TP53 mutation status", legend = "TP53", colours = c("#1683E6", "#FF9000"), stem = "fig9_km_tp53"),
  list(variable = "tmb_tertile", title = "TCGA-BRCA OS by TMB tertile", legend = "TMB tertile", colours = c("#1683E6", "#8CAF43", "#FF9000"), stem = "fig10_km_tmb_tertiles")
)

draw_km <- function(spec, fit, p_value) {
  groups <- levels(droplevels(analysis_data[[spec$variable]]))
  times <- seq(0, 25, 5)
  layout(matrix(c(1, 2), ncol = 1), heights = c(3.2, 1.2))
  par(mar = c(1, 5, 4, 2))
  plot(fit, col = spec$colours, lwd = 2, mark.time = TRUE, xlim = c(0, 25), ylim = c(0, 1),
       xlab = "", ylab = "Overall survival", main = spec$title, xaxt = "n")
  legend("top", legend = groups, col = spec$colours, lty = 1, lwd = 2, horiz = TRUE, bty = "n", title = spec$legend)
  text(1, 0.18, sprintf("Log-rank p = %.3f", p_value), adj = 0, cex = 1.1)

  par(mar = c(4, 7, 1, 2))
  plot.new()
  plot.window(xlim = c(-7, 25), ylim = c(0.5, length(groups) + 0.8))
  axis(1, at = times)
  title(xlab = "Time (years)")
  mtext("Number at risk", side = 3, line = -0.2, adj = 0, font = 2)
  for (index in seq_along(groups)) {
    group_data <- analysis_data[analysis_data[[spec$variable]] == groups[index] & !is.na(analysis_data[[spec$variable]]), ]
    group_fit <- survfit(Surv(time_years, event) ~ 1, data = group_data)
    risks <- summary(group_fit, times = times, extend = TRUE)$n.risk
    y <- length(groups) - index + 1
    text(-1.4, y, groups[index], adj = 1, col = spec$colours[index], font = 2)
    text(times, y, risks)
  }
  box()
}

km_rows <- lapply(km_specs, function(spec) {
  formula <- as.formula(paste("Surv(time_years, event) ~", spec$variable))
  fit <- survfit(formula, data = analysis_data)
  p_value <- log_rank_p(formula, analysis_data)
  png(file.path(figures_dir, paste0(spec$stem, ".png")), width = 2400, height = 2000, res = 300)
  draw_km(spec, fit, p_value); dev.off()
  svg(file.path(figures_dir, paste0(spec$stem, ".svg")), width = 8, height = 7)
  draw_km(spec, fit, p_value); dev.off()
  counts <- table(analysis_data[[spec$variable]], useNA = "no")
  group_text <- paste0(names(counts), " (", as.integer(counts), ")", collapse = ", ")
  data.frame(variable = spec$variable, groups = group_text, log_rank_p = p_value)
})
km_results <- do.call(rbind, km_rows)
write.csv(km_results, file.path(tables_dir, "supplementary_table9_km_summary.csv"), row.names = FALSE)

cat("Generated Supplementary Table 9 and Figures 7-10.\n")
