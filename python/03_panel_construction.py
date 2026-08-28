"""Packages: pandas, scipy. Prerequisites: AURORA candidates, TCGA raw files, and an AURORA burden input."""

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon

from common import AURORA_DIR, TCGA_RAW_DIR, ensure_output_dirs, require_columns, tcga_sample_type


def valid_os_patients(clinical):
    time = pd.to_numeric(clinical["OS_MONTHS"], errors="coerce")
    status = clinical["OS_STATUS"].fillna("").str.upper()
    valid_status = status.str.contains("DECEASED") | status.str.contains("LIVING")
    return set(clinical.loc[time.gt(0) & valid_status, "patientId"])


def burden_from_restricted_mutations(path, panel_genes):
    data = pd.read_csv(path, sep="\t")
    require_columns(data, ["patient_id", "sample_type", "gene"], path)
    patient_types = data.groupby("patient_id")["sample_type"].agg(set)
    paired = sorted(
        patient_id
        for patient_id, types in patient_types.items()
        if {"Primary", "Metastasis"}.issubset(types)
    )
    rows = []
    for patient_id in paired:
        patient = data[data.patient_id.eq(patient_id)]
        primary = set(patient.loc[patient.sample_type.eq("Primary"), "gene"])
        metastasis = set(patient.loc[patient.sample_type.eq("Metastasis"), "gene"])
        rows.append(
            {
                "primary_burden": sum(gene in primary for gene in panel_genes),
                "metastatic_burden": sum(gene in metastasis for gene in panel_genes),
            }
        )
    pairs = pd.DataFrame(rows)
    return pairs.value_counts(sort=False).rename("patient_count").reset_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--aurora-raw", type=Path)
    args = parser.parse_args()
    ensure_output_dirs()

    candidates = pd.read_csv(AURORA_DIR / "aurora_candidate_genes.csv")
    clinical = pd.read_csv(TCGA_RAW_DIR / "tcga_brca_clinical_raw.csv")
    mutations = pd.read_csv(TCGA_RAW_DIR / "tcga_brca_mutations_raw.csv")
    require_columns(candidates, ["gene", "n_met", "direction"], "aurora_candidate_genes.csv")
    require_columns(clinical, ["patientId", "OS_MONTHS", "OS_STATUS"], "tcga_brca_clinical_raw.csv")
    require_columns(mutations, ["sampleId", "patientId", "gene"], "tcga_brca_mutations_raw.csv")

    candidates = candidates[(candidates.n_met >= 5) & candidates.direction.eq("met_enriched")].copy()
    os_ids = valid_os_patients(clinical)
    primary_mutations = mutations[tcga_sample_type(mutations.sampleId).eq("01")]
    counts = (
        primary_mutations[primary_mutations.patientId.isin(os_ids) & primary_mutations.gene.isin(candidates.gene)]
        .groupby("gene")["patientId"]
        .nunique()
    )
    candidates["tcga_mutated_n"] = candidates.gene.map(counts).fillna(0).astype(int)
    candidates["included_in_panel"] = candidates.tcga_mutated_n.ge(10)
    candidates.to_csv(AURORA_DIR / "panel_selection.csv", index=False)

    panel_genes = candidates.loc[candidates.included_in_panel, "gene"].tolist()
    manuscript_order = ["TP53", "FLG", "PIK3CA", "CACNA1A", "GOLGB1", "COL14A1", "MUC5B"]
    if set(panel_genes) != set(manuscript_order):
        raise ValueError(f"Derived panel does not match the dissertation: {panel_genes}")
    pd.DataFrame({"gene": manuscript_order}).to_csv(AURORA_DIR / "final_panel_genes.csv", index=False)

    burden_file = AURORA_DIR / "aurora_paired_burden_counts.csv"
    if args.aurora_raw:
        burden_counts = burden_from_restricted_mutations(args.aurora_raw, manuscript_order)
        burden_counts.to_csv(burden_file, index=False)
    else:
        burden_counts = pd.read_csv(burden_file)
    require_columns(burden_counts, ["primary_burden", "metastatic_burden", "patient_count"], burden_file)
    expanded = burden_counts.loc[burden_counts.index.repeat(burden_counts.patient_count)].reset_index(drop=True)
    if len(expanded) != 39:
        raise ValueError(f"Expected 39 AURORA burden pairs, found {len(expanded)}")
    _, p_value = wilcoxon(expanded.metastatic_burden, expanded.primary_burden)
    if abs(p_value - 0.001232) > 0.000001:
        raise ValueError(f"Unexpected AURORA Wilcoxon p value: {p_value}")

    print("Final panel:", ", ".join(manuscript_order))
    print(f"AURORA paired burden Wilcoxon p: {p_value:.6f}")


if __name__ == "__main__":
    main()
