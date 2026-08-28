"""Packages: pandas, numpy, scipy. Prerequisite: restricted AURORA protein-altering mutation TSV."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2

from common import AURORA_DIR, ensure_output_dirs, require_columns


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def build_matrix(data):
    rows = []
    for patient_id, patient in data.groupby("patient_id", sort=True):
        primary = set(patient.loc[patient["sample_type"] == "Primary", "gene"].dropna())
        metastasis = set(patient.loc[patient["sample_type"] == "Metastasis", "gene"].dropna())
        for gene in sorted(primary | metastasis):
            rows.append(
                {
                    "patient_id": patient_id,
                    "gene": gene,
                    "primary_mutated": int(gene in primary),
                    "met_mutated": int(gene in metastasis),
                }
            )
    return pd.DataFrame(rows)


def run_mcnemar(matrix, paired_n):
    records = []
    for gene, gene_data in matrix.groupby("gene", sort=True):
        n_primary = int(gene_data["primary_mutated"].sum())
        n_metastasis = int(gene_data["met_mutated"].sum())
        if max(n_primary, n_metastasis) < 2:
            continue

        both = int(((gene_data.primary_mutated == 1) & (gene_data.met_mutated == 1)).sum())
        primary_only = int(((gene_data.primary_mutated == 1) & (gene_data.met_mutated == 0)).sum())
        metastasis_only = int(((gene_data.primary_mutated == 0) & (gene_data.met_mutated == 1)).sum())
        discordant = primary_only + metastasis_only
        both_wild_type = paired_n - both - primary_only - metastasis_only

        if discordant == 0:
            statistic = np.nan
            p_value = np.nan
        else:
            statistic = (abs(primary_only - metastasis_only) - 1) ** 2 / discordant
            if discordant < 25:
                p_value = binomtest(metastasis_only, discordant, 0.5).pvalue
            else:
                p_value = float(chi2.sf(statistic, 1))

        if metastasis_only > primary_only:
            direction = "met_enriched"
        elif primary_only > metastasis_only:
            direction = "met_depleted"
        else:
            direction = "no_change"

        records.append(
            {
                "gene": gene,
                "a_both_mut": both,
                "b_prim_only": primary_only,
                "c_met_only": metastasis_only,
                "d_both_wt": both_wild_type,
                "n_prim": n_primary,
                "n_met": n_metastasis,
                "prim_freq": n_primary / paired_n,
                "met_freq": n_metastasis / paired_n,
                "freq_diff": (n_metastasis - n_primary) / paired_n,
                "discordant": discordant,
                "mcnemar_stat": statistic,
                "mcnemar_p": p_value,
                "direction": direction,
            }
        )

    results = pd.DataFrame(records)
    valid = results["mcnemar_p"].notna()
    results.loc[valid, "mcnemar_fdr"] = bh_adjust(results.loc[valid, "mcnemar_p"])
    return results.sort_values(["mcnemar_p", "gene"], na_position="last")


def cohort_summary(data, paired_ids):
    sample_counts = data.groupby("patient_id")["sample_id"].nunique()
    patient_types = data.groupby("patient_id")["sample_type"].agg(set)
    primary_only = sum(types == {"Primary"} for types in patient_types)
    metastasis_only = sum(types == {"Metastasis"} for types in patient_types)
    return pd.DataFrame(
        [
            ("Total patients", data.patient_id.nunique()),
            ("Total samples", data.sample_id.nunique()),
            ("Total protein-altering mutations", len(data)),
            ("Unique genes with mutations", data.gene.nunique()),
            ("Paired patients (primary + metastasis)", len(paired_ids)),
            ("Primary-only patients", primary_only),
            ("Metastasis-only patients", metastasis_only),
            ("Samples per patient median", sample_counts.median()),
            ("Samples per patient minimum", sample_counts.min()),
            ("Samples per patient maximum", sample_counts.max()),
        ],
        columns=["characteristic", "value"],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=AURORA_DIR / "aurora_protein_altering_mutations.tsv",
    )
    args = parser.parse_args()
    ensure_output_dirs()

    data = pd.read_csv(args.input, sep="\t")
    require_columns(data, ["patient_id", "sample_id", "sample_type", "gene"], args.input)
    patient_types = data.groupby("patient_id")["sample_type"].agg(set)
    paired_ids = sorted(
        patient_id
        for patient_id, types in patient_types.items()
        if {"Primary", "Metastasis"}.issubset(types)
    )
    if len(paired_ids) != 39:
        raise ValueError(f"Expected 39 paired AURORA patients, found {len(paired_ids)}")

    paired = data[data.patient_id.isin(paired_ids)].copy()
    matrix = build_matrix(paired)
    results = run_mcnemar(matrix, len(paired_ids))
    candidates = results[(results.direction == "met_enriched") & (results.n_met >= 5)]

    if len(results) != 7503 or results.mcnemar_p.notna().sum() != 7491:
        raise ValueError("AURORA gene-test totals do not match the dissertation cohort")
    if len(candidates) != 10:
        raise ValueError(f"Expected 10 AURORA candidates, found {len(candidates)}")

    summary = cohort_summary(data, paired_ids)
    summary_values = summary.set_index("characteristic").value
    expected_totals = {
        "Total patients": 49,
        "Total samples": 134,
        "Total protein-altering mutations": 47776,
        "Unique genes with mutations": 13648,
        "Paired patients (primary + metastasis)": 39,
        "Primary-only patients": 1,
        "Metastasis-only patients": 9,
    }
    mismatches = {
        label: (summary_values[label], expected)
        for label, expected in expected_totals.items()
        if int(summary_values[label]) != expected
    }
    if mismatches:
        raise ValueError(f"AURORA normalised input totals do not match the dissertation: {mismatches}")

    matrix.to_csv(AURORA_DIR / "matrix_union.csv", index=False)
    results.to_csv(AURORA_DIR / "aurora_mcnemar_results.csv", index=False)
    candidates.to_csv(AURORA_DIR / "aurora_candidate_genes.csv", index=False)
    summary.to_csv(AURORA_DIR / "aurora_cohort_summary.csv", index=False)
    print(f"Paired patients: {len(paired_ids)}")
    print(f"Genes tested: {len(results)}; valid p values: {results.mcnemar_p.notna().sum()}")
    print("AURORA candidates:", ", ".join(candidates.gene))


if __name__ == "__main__":
    main()
