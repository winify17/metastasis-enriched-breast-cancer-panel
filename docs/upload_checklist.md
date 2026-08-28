# GitHub release checklist

Complete each item before making the repository public.

- Apply the dissertation corrections listed in `method_alignment.md`.
- Obtain written AURORA disclosure approval. Confirm whether genome-wide small counts and paired burden cells may be released.
- Complete `data/aurora/source_provenance_template.yml` from the authorised source documentation.
- Keep all restricted or uncleared AURORA CSV and TSV files excluded from version control.
- Install the pinned Python packages and run `python python/run_public_pipeline.py` on a host with Matplotlib.
- Run `Rscript R/run_all.R`. Confirm that both validators report that all displayed manuscript values passed.
- Inspect all ten regenerated figures, including censoring marks, risk tables, confidence intervals, p-values, and axis labels.
- Confirm that the working tree contains no `results/`, `__pycache__/`, `.Rhistory`, patient-level AURORA file, credentials, tokens, or local paths.
- Verify included public inputs against `data/SHA256SUMS`; verify controlled AURORA inputs separately in the authorised environment.
- Review `docs/output_manifest.md` against the final dissertation captions and numbering.

The default `.gitignore` deliberately excludes generated results and all AURORA CSV files except the published final gene list. Do not force-add an AURORA file until its exact release status has been approved.
