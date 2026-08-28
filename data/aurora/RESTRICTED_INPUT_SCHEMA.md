# Restricted AURORA input and provenance

The authorised source export must be normalised to `aurora_protein_altering_mutations.tsv` before running `python/01_aurora_mcnemar.py`. The normalised file remains restricted and is excluded from version control.

## Required analysis columns

| Column | Definition |
|---|---|
| `patient_id` | Stable AURORA patient identifier |
| `sample_id` | Stable tumour sample identifier |
| `sample_type` | Exactly `Primary` or `Metastasis` after metadata mapping |
| `gene` | Approved gene symbol used for gene-level aggregation |

The restricted audit copy should also retain `consequence`, `chromosome`, `position`, `reference_allele`, `alternate_allele`, and the source variant identifier when those fields are available. They are not required by the paired gene-level script, but they are needed to audit filtering and deduplication.

## Required restricted preprocessing record

1. Record the source study accession or release, export date, source file checksum, access approval, and data custodian.
2. Map each sample to one patient and one disease compartment using the authorised sample metadata. Normalise compartments to `Primary` or `Metastasis` and record any excluded or ambiguous samples.
3. Retain somatic variants with an eligible protein-altering consequence: missense, nonsense or stop-gained, frameshift, splice-site, in-frame insertion or deletion, translation-start or start-lost, and nonstop or stop-lost variants.
4. Apply no additional variant allele fraction threshold.
5. Define and record the variant-level deduplication key before gene-level aggregation. Transcript duplicates must not create additional patient-gene mutation states. If a single variant has several transcript consequences, document how the retained consequence is chosen.
6. Confirm the normalised export totals before analysis: 49 patients, 134 samples, 47,776 eligible mutation records, and 13,648 mutated genes. The paired subset must contain 39 patients with at least one primary and one metastatic sample.

The source-specific column mapping and deduplication rule cannot be inferred from the derived repository files. Complete `source_provenance_template.yml` from the authorised source documentation before claiming full raw-data reproducibility.
