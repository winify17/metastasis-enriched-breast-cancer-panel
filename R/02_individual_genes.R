# Package: survival. Prerequisite: source 00_setup.R.

if (!exists("analysis_data")) source("00_setup.R")

extract_gene_term <- function(model, term) {
  model_summary <- summary(model)
  c(
    HR = model_summary$coefficients[term, "exp(coef)"],
    lower = model_summary$conf.int[term, "lower .95"],
    upper = model_summary$conf.int[term, "upper .95"],
    p = model_summary$coefficients[term, "Pr(>|z|)"],
    cindex = unname(model_summary$concordance[1])
  )
}

gene_results <- lapply(panel_genes, function(gene) {
  term <- paste0(gene, "_mutated")
  univariable <- coxph(as.formula(paste("Surv(time_years, event) ~", term)), data = analysis_data, na.action = na.exclude)
  clinical <- coxph(as.formula(paste("Surv(time_years, event) ~", term, "+ age + stage + mol_subtype")), data = analysis_data, na.action = na.exclude)
  tmb_adjusted <- coxph(as.formula(paste("Surv(time_years, event) ~", term, "+ log_tmb + age + stage + mol_subtype")), data = analysis_data, na.action = na.exclude)
  u <- extract_gene_term(univariable, term)
  c <- extract_gene_term(clinical, term)
  t <- extract_gene_term(tmb_adjusted, term)
  data.frame(
    gene = gene,
    n_mut = sum(analysis_data[[term]] == 1),
    n_wt = sum(analysis_data[[term]] == 0),
    uni_HR = u["HR"], uni_p = u["p"],
    multi_HR = c["HR"], multi_lower = c["lower"], multi_upper = c["upper"], multi_p = c["p"], multi_cindex = c["cindex"],
    tmb_HR = t["HR"], tmb_lower = t["lower"], tmb_upper = t["upper"], tmb_p = t["p"], tmb_cindex = t["cindex"],
    multi_n = clinical$n, multi_events = clinical$nevent,
    tmb_n = tmb_adjusted$n, tmb_events = tmb_adjusted$nevent,
    row.names = NULL
  )
})
gene_results <- do.call(rbind, gene_results)
write.csv(gene_results, file.path(tables_dir, "supplementary_table7_individual_genes.csv"), row.names = FALSE)

draw_gene_forest <- function() {
  y <- seq_len(nrow(gene_results))
  clinical_y <- y + 0.12
  tmb_y <- y - 0.12
  labels <- paste0(gene_results$gene, " (n=", gene_results$n_mut, ")")
  clinical_p <- ifelse(gene_results$multi_p < 0.001, "p<0.001", sprintf("p=%.3f", gene_results$multi_p))
  tmb_p <- ifelse(gene_results$tmb_p < 0.001, "p<0.001", sprintf("p=%.3f", gene_results$tmb_p))
  clinical_p <- paste(clinical_p, ifelse(gene_results$multi_p < 0.05, "*", ""))
  tmb_p <- paste(tmb_p, ifelse(gene_results$tmb_p < 0.05, "*", ""))
  par(mar = c(8, 9, 4, 2))
  plot(gene_results$multi_HR, clinical_y, xlim = c(0, 15.5), ylim = c(0.5, 7.5),
       pch = 15, col = "#1683E6", yaxt = "n", ylab = "", xlab = "Hazard ratio (95% CI)",
       main = "Individual gene Cox analysis")
  axis(2, at = y, labels = labels, las = 1)
  segments(gene_results$multi_lower, clinical_y, gene_results$multi_upper, clinical_y, col = "#1683E6", lwd = 2)
  points(gene_results$tmb_HR, tmb_y, pch = 18, col = "#FF9000")
  segments(gene_results$tmb_lower, tmb_y, gene_results$tmb_upper, tmb_y, col = "#FF9000", lwd = 2)
  text(12.2, clinical_y, clinical_p, adj = 0, col = "#1683E6", family = "mono", cex = 0.85)
  text(12.2, tmb_y, tmb_p, adj = 0, col = "#FF9000", family = "mono", cex = 0.85)
  abline(v = 1, lty = 2, col = "grey50")
  legend("bottom", inset = c(0, -0.32), xpd = NA, horiz = TRUE,
         legend = c("Clinical adjustment", "Clinical and TMB adjustment"),
         col = c("#1683E6", "#FF9000"), pch = c(15, 18), lty = 1, bty = "n")
  mtext("TCGA-BRCA overall survival", side = 3, line = 0.4)
}
png(file.path(figures_dir, "fig5_individual_gene_forest.png"), width = 3000, height = 1800, res = 300)
draw_gene_forest(); dev.off()
svg(file.path(figures_dir, "fig5_individual_gene_forest.svg"), width = 10, height = 6)
draw_gene_forest(); dev.off()

cat("Generated Supplementary Table 7 and Figure 5.\n")
