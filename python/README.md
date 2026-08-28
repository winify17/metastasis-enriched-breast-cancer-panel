# Python workflow

Required packages: `pandas`, `numpy`, `scipy`, and `matplotlib`.

Required public inputs: the CSV files under `data/aurora/` and `data/tcga_raw/`.

Run the archived public workflow from the repository root:

```bash
python python/run_public_pipeline.py
```

The optional restricted-data workflow begins with `01_aurora_mcnemar.py`. `02_tcga_download.py` requires internet access and replaces the archived TCGA snapshot. Script `03_panel_construction.py` applies the TCGA patient-count threshold, `04_tcga_cleaning.py` rebuilds the survival input, and `05_manuscript_descriptives.py` creates the descriptive tables and figures.
