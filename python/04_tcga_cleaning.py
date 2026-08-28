"""Packages: pandas, numpy, scipy. Prerequisites: TCGA clinical, mutation, and TMB CSV files plus final panel genes."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import AURORA_DIR, TCGA_PROCESSED_DIR, TCGA_RAW_DIR, ensure_output_dirs, require_columns, tcga_sample_type


def parse_stage(value):
    stage_groups = {
        "I": {"Stage I", "Stage IA", "Stage IB"},
        "II": {"Stage II", "Stage IIA", "Stage IIB"},
        "III": {"Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"},
        "IV": {"Stage IV"},
    }
    text = str(value).strip() if pd.notna(value) else ""
    for group, labels in stage_groups.items():
        if text in labels:
            return group
    return np.nan


def clean_receptor(value):
    text = str(value).strip().lower() if pd.notna(value) else ""
    if "positive" in text:
        return "Positive"
    if "negative" in text:
        return "Negative"
    return np.nan


def derive_subtype(row):
    er, pr, her2 = row.er_status, row.pr_status, row.her2_status
    if pd.isna(er) and pd.isna(pr) and pd.isna(her2):
        return np.nan
    hr_positive = er == "Positive" or pr == "Positive"
    if hr_positive and her2 == "Negative":
        return "HR+/HER2-"
    if hr_positive and her2 == "Positive":
        return "HR+/HER2+"
    if not hr_positive and her2 == "Positive":
        return "HER2+"
    if not hr_positive and her2 == "Negative":
        return "Triple Negative"
    return np.nan


def main():
    ensure_output_dirs()
    clinical_path = TCGA_RAW_DIR / "tcga_brca_clinical_raw.csv"
    mutation_path = TCGA_RAW_DIR / "tcga_brca_mutations_raw.csv"
    tmb_path = TCGA_RAW_DIR / "tcga_brca_tmb.csv"
    panel_path = AURORA_DIR / "final_panel_genes.csv"
    clinical = pd.read_csv(clinical_path)
    mutations = pd.read_csv(mutation_path)
    tmb = pd.read_csv(tmb_path)
    panel_genes = pd.read_csv(panel_path).gene.tolist()

    require_columns(
        clinical,
        ["patientId", "OS_MONTHS", "OS_STATUS", "AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "ER_STATUS_BY_IHC", "PR_STATUS_BY_IHC", "IHC_HER2"],
        clinical_path,
    )
    require_columns(mutations, ["sampleId", "patientId", "gene"], mutation_path)
    require_columns(tmb, ["sampleId", "tmb"], tmb_path)
    if len(panel_genes) != 7:
        raise ValueError(f"Expected seven panel genes, found {len(panel_genes)}")

    clinical["time_years"] = pd.to_numeric(clinical.OS_MONTHS, errors="coerce") / 12
    status = clinical.OS_STATUS.fillna("").str.upper()
    clinical["event"] = np.select(
        [status.str.contains("DECEASED"), status.str.contains("LIVING")],
        [1, 0],
        default=np.nan,
    )
    clinical["age"] = pd.to_numeric(clinical.AGE, errors="coerce")
    clinical["stage"] = clinical.AJCC_PATHOLOGIC_TUMOR_STAGE.map(parse_stage)
    clinical["er_status"] = clinical.ER_STATUS_BY_IHC.map(clean_receptor)
    clinical["pr_status"] = clinical.PR_STATUS_BY_IHC.map(clean_receptor)
    clinical["her2_status"] = clinical.IHC_HER2.map(clean_receptor)
    clinical["mol_subtype"] = clinical.apply(derive_subtype, axis=1)
    clinical = clinical.rename(columns={"patientId": "sample_id"})
    clinical = clinical[
        clinical.time_years.notna() & clinical.event.notna() & clinical.time_years.gt(0)
    ].copy()

    primary_mutations = mutations[tcga_sample_type(mutations.sampleId).eq("01")].copy()
    for gene in panel_genes:
        mutated_ids = set(primary_mutations.loc[primary_mutations.gene.eq(gene), "patientId"])
        clinical[f"{gene}_mutated"] = clinical.sample_id.isin(mutated_ids).astype(int)
    mutation_columns = [f"{gene}_mutated" for gene in panel_genes]
    clinical["panel_burden"] = clinical[mutation_columns].sum(axis=1)
    clinical["panel_binary"] = clinical.panel_burden.gt(0).astype(int)

    tmb["patientId"] = tmb.sampleId.str[:12]
    tmb["sample_type"] = tcga_sample_type(tmb.sampleId)
    primary_tmb = tmb[tmb.sample_type.eq("01")].sort_values(["patientId", "sampleId"]).copy()
    duplicate_primary = primary_tmb.groupby("patientId").tmb.nunique()
    if duplicate_primary.gt(1).any():
        patients = ", ".join(duplicate_primary[duplicate_primary.gt(1)].index)
        raise ValueError(f"Multiple primary samples have different TMB values for: {patients}")
    patient_tmb = primary_tmb.drop_duplicates("patientId", keep="first")[["patientId", "tmb"]]
    patient_tmb = patient_tmb.rename(columns={"patientId": "sample_id"})
    clinical = clinical.merge(patient_tmb, on="sample_id", how="left")

    output_columns = [
        "sample_id", "time_years", "event", "age", "stage", "mol_subtype",
        "panel_burden", "panel_binary", "tmb", *mutation_columns,
    ]
    survival_input = clinical[output_columns].copy()
    survival_input.to_csv(TCGA_PROCESSED_DIR / "tcga_brca_panel_survival_input.csv", index=False)

    complete = survival_input[survival_input.tmb.notna()]
    rho, p_value = spearmanr(complete.panel_burden, complete.tmb)
    expected_counts = {"TP53": 296, "FLG": 43, "PIK3CA": 313, "CACNA1A": 17, "GOLGB1": 13, "COL14A1": 16, "MUC5B": 35}
    observed_counts = {gene: int(survival_input[f"{gene}_mutated"].sum()) for gene in panel_genes}
    if len(survival_input) != 1083 or int(survival_input.event.sum()) != 151:
        raise ValueError("TCGA OS cohort does not match the dissertation (expected N=1083, events=151)")
    if observed_counts != expected_counts:
        raise ValueError(f"Panel mutation counts do not match the dissertation: {observed_counts}")
    if abs(rho - 0.321) > 0.001:
        raise ValueError(f"Panel-TMB correlation does not match the dissertation: rho={rho}")

    print(f"TCGA OS cohort: {len(survival_input)} patients, {int(survival_input.event.sum())} events")
    print(f"TMB available: {survival_input.tmb.notna().sum()}; Spearman rho={rho:.3f}, p={p_value:.3g}")


if __name__ == "__main__":
    main()
