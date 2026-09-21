# Metastasis-enriched mutation panel

This repository reproduces the statistical workflow for a seven-gene mutation panel derived from paired primary and metastatic breast tumours in AURORA and assessed for overall survival in TCGA-BRCA.

Panel genes: `TP53`, `FLG`, `PIK3CA`, `CACNA1A`, `GOLGB1`, `COL14A1`, and `MUC5B`.

## Scope

The workflow generates all 10 statistical figures, all 10 in-text tables, and all 10 supplementary tables described in the dissertation. Generated files are written to `results/`, which is excluded from version control.

Before publication, AURORA-derived files are retained in the local analysis workspace only to document the methods and demonstrate the reproducibility workflow. They are not included as an independent data release, and access and reuse remain subject to the terms of the original AURORA source.

## Requirements

- Python 3.11 or later with the packages in `requirements.txt`
- R 4.4 or later with `survival`; the validated core environment is recorded in `renv.lock`

## Reproduction with cleared AURORA inputs

From the repository root:

```bash
python -m pip install -r requirements.txt
python python/run_public_pipeline.py
Rscript R/run_all.R
```

The Python pipeline reconstructs the seven-gene panel, rebuilds the TCGA analysis input, and generates Tables 1-4, Supplementary Tables 1-2, and Figures 1-3. The R pipeline fits the survival models and generates Tables 5-10, Supplementary Tables 3-7, 8A, 8B, and 9, and Figures 4-10. Tables 7-10 contain the numbers at risk displayed within Figures 7-10. Both pipelines stop when manuscript values differ from their expected values.

The default GitHub release is code-first and does not track the AURORA aggregate CSV files. The AURORA Supplementary Tables 1-8 can be downloaded from the [Garcia-Recio et al. Nature Cancer article](https://www.nature.com/articles/s43018-022-00491-x) or as a [direct supplementary archive](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs43018-022-00491-x/MediaObjects/43018_2022_491_MOESM2_ESM.rar). Supplementary Table 2 contains the sample-level molecular information used to identify the WES samples and classify them as primary or metastatic.

## Restricted AURORA reproduction

Place the authorised protein-altering mutation file at `data/aurora/aurora_protein_altering_mutations.tsv`, then run:

```bash
python python/01_aurora_mcnemar.py
python python/02_tcga_download.py
python python/03_panel_construction.py --aurora-raw data/aurora/aurora_protein_altering_mutations.tsv
python python/04_tcga_cleaning.py
python python/05_manuscript_descriptives.py
Rscript R/run_all.R
```

`02_tcga_download.py` queries the public cBioPortal API and may return a later data snapshot. The archived inputs in this repository are retained to reproduce the dissertation values exactly.

## Analysis definition

Ten AURORA candidates were selected by requiring metastatic enrichment and at least five mutation-positive metastatic patients. Seven candidates remained after requiring at least ten mutation-positive patients in the TCGA-BRCA overall-survival cohort. The paired AURORA burden comparison is descriptive because panel selection and burden comparison use the same 39 pairs.

See `docs/output_manifest.md` for the complete table and figure map. See `docs/method_alignment.md` for the dissertation corrections needed to align the written account with the reproducible analysis.

## Data use

TCGA data are public through cBioPortal. AURORA raw data, patient-level derivatives, and small-cell aggregate files are not part of the default public package. Users remain responsible for the terms attached to each source dataset.

## License

Code is released under the MIT License. Source datasets retain their original access conditions and licences.
