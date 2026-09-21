# R analysis

Required package: `survival`. All figures use the R graphics devices supplied with R.

Required input: `data/tcga_processed/tcga_brca_panel_survival_input.csv`.

From the repository root, run:

```bash
Rscript R/run_all.R
```

The scripts create manuscript Tables 5-10, Supplementary Tables 3-7, 8A, 8B, and 9, and Figures 4-10 in `results/`. They also write the time-specific extended Cox estimates reported in Section 3.8.

The numerical workflow was checked with R 4.5.3 and `survival` 3.8-6. Exact core package versions are recorded in `renv.lock`.
