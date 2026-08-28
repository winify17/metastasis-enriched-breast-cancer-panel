# Manuscript output manifest

All generated files are written to `results/` and are excluded from the public repository.

| Manuscript item | Output stem | Generator | Required input |
|---|---|---|---|
| Figure 1 | `fig1_panel_construction` | `python/05_manuscript_descriptives.py` | AURORA McNemar results, final panel |
| Figure 2 | `fig2_burden_boxplot` | `python/05_manuscript_descriptives.py` | Aggregated paired AURORA burdens |
| Figure 3 | `fig3_tcga_distribution_tmb` | `python/05_manuscript_descriptives.py` | TCGA survival input |
| Figure 4 | `fig4_forest_model3` | `R/01_cox_models.R` | TCGA survival input |
| Figure 5 | `fig5_individual_gene_forest` | `R/02_individual_genes.R` | TCGA survival input |
| Figure 6 | `fig6_cindex_comparison` | `R/01_cox_models.R` | TCGA survival input |
| Figure 7 | `fig7_km_panel_burden_groups` | `R/03_km_curves.R` | TCGA survival input |
| Figure 8 | `fig8_km_panel_binary` | `R/03_km_curves.R` | TCGA survival input |
| Figure 9 | `fig9_km_tp53` | `R/03_km_curves.R` | TCGA survival input |
| Figure 10 | `fig10_km_tmb_tertiles` | `R/03_km_curves.R` | TCGA survival input |
| Table 1 | `table1_aurora_cohort_summary.csv` | `python/05_manuscript_descriptives.py` | AURORA cohort summary |
| Table 2 | `table2_mcnemar_summary.csv` | `python/05_manuscript_descriptives.py` | AURORA McNemar results |
| Table 3 | `table3_panel_composition.csv` | `python/05_manuscript_descriptives.py` | AURORA McNemar results, panel selection, TCGA survival input |
| Table 4 | `table4_aurora_burden_summary.csv` | `python/05_manuscript_descriptives.py` | Aggregated paired AURORA burdens |
| Table 5 | `table5_cox_model3.csv` | `R/01_cox_models.R` | TCGA survival input |
| Table 6 | `table6_sensitivity_configurations.csv` | `R/04_sensitivity_analysis.R` | TCGA survival input |
| Supplementary Table 1 | `supplementary_table1_tcga_summary.csv` | `python/05_manuscript_descriptives.py` | TCGA survival input |
| Supplementary Table 2 | `supplementary_table2_panel_distribution.csv` | `python/05_manuscript_descriptives.py` | TCGA survival input |
| Supplementary Table 3 | `supplementary_table3_model_comparison.csv` | `R/01_cox_models.R` | TCGA survival input |
| Supplementary Table 4 | `supplementary_table4_cox_model1.csv` | `R/01_cox_models.R` | TCGA survival input |
| Supplementary Table 5 | `supplementary_table5_cox_model2.csv` | `R/01_cox_models.R` | TCGA survival input |
| Supplementary Table 6 | `supplementary_table6_cox_model4.csv` | `R/01_cox_models.R` | TCGA survival input |
| Supplementary Table 7 | `supplementary_table7_individual_genes.csv` | `R/02_individual_genes.R` | TCGA survival input |
| Supplementary Table 8 | `supplementary_table8_ph_tests.csv` | `R/01_cox_models.R` | TCGA survival input |
| Supplementary Table 9 | `supplementary_table9_km_summary.csv` | `R/03_km_curves.R` | TCGA survival input |

`python/run_public_pipeline.py` runs the public Python stages. `R/run_all.R` runs the complete survival workflow and checks the displayed manuscript estimates.
