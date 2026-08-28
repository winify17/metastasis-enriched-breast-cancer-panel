"""Package: pandas. Prerequisite: generated Python manuscript tables and, optionally, Figures 1-3."""

import argparse

import pandas as pd

from common import FIGURES_DIR, TABLES_DIR


def assert_equal(observed, expected, label):
    if observed != expected:
        raise AssertionError(f"{label}: observed {observed!r}, expected {expected!r}")


def assert_rounded(observed, expected, digits, label):
    if pd.isna(observed) or round(float(observed), digits) != round(float(expected), digits):
        raise AssertionError(
            f"{label}: observed {observed!r}, expected {expected!r} at {digits} decimal places"
        )


def validate_table1():
    observed = pd.read_csv(TABLES_DIR / "table1_aurora_cohort_summary.csv").set_index("characteristic").value
    expected = {
        "Total patients": 49,
        "Total samples": 134,
        "Total protein-altering mutations": 47776,
        "Unique genes with mutations": 13648,
        "Paired patients (primary + metastasis)": 39,
        "Primary-only patients": 1,
        "Metastasis-only patients": 9,
        "Samples per patient median": 2,
        "Samples per patient minimum": 1,
        "Samples per patient maximum": 9,
    }
    assert_equal(set(observed.index), set(expected), "Table 1 row labels")
    for label, value in expected.items():
        assert_rounded(observed[label], value, 0, f"Table 1 {label}")


def validate_table2():
    observed = pd.read_csv(TABLES_DIR / "table2_mcnemar_summary.csv").set_index("parameter").value
    expected_counts = {
        "Genes tested": 7503,
        "Genes with valid p-value": 7491,
        "Met-enriched genes (FDR < 0.05)": 0,
        "Met-depleted genes (FDR < 0.05)": 0,
    }
    for label, value in expected_counts.items():
        assert_rounded(observed[label], value, 0, f"Table 2 {label}")
    assert_equal(observed["Most significant gene"], "TP53", "Table 2 most significant gene")
    assert_rounded(observed["Most significant McNemar p"], 0.000244, 6, "Table 2 TP53 McNemar p")
    assert_rounded(observed["Most significant FDR"], 0.366, 3, "Table 2 TP53 FDR")


def validate_table3():
    observed = pd.read_csv(TABLES_DIR / "table3_panel_composition.csv").set_index("gene")
    expected = {
        "TP53": (13, 26, 33.3, 66.7, 33.3, 0.0002, 4, 0.366, 296, 27.3),
        "FLG": (5, 11, 12.8, 28.2, 15.4, 0.180, 3, 0.710, 43, 4.0),
        "PIK3CA": (5, 9, 12.8, 23.1, 10.3, 0.125, 3, 0.710, 313, 28.9),
        "CACNA1A": (4, 8, 10.3, 20.5, 10.3, 0.219, 3, 0.710, 17, 1.6),
        "GOLGB1": (2, 7, 5.1, 17.9, 12.8, 0.180, 3, 0.710, 13, 1.2),
        "COL14A1": (4, 7, 10.3, 17.9, 7.7, 0.508, 3, 0.720, 16, 1.5),
        "MUC5B": (4, 5, 10.3, 12.8, 2.6, 1.000, 3, 1.000, 35, 3.2),
    }
    assert_equal(set(observed.index), set(expected), "Table 3 genes")
    for gene, values in expected.items():
        n_primary, n_metastasis, primary_percent, metastatic_percent, difference_percent, p_value, p_digits, fdr, tcga_count, tcga_percent = values
        row = observed.loc[gene]
        assert_rounded(row.n_prim, n_primary, 0, f"Table 3 {gene} primary count")
        assert_rounded(row.n_met, n_metastasis, 0, f"Table 3 {gene} metastatic count")
        assert_rounded(100 * row.prim_freq, primary_percent, 1, f"Table 3 {gene} primary percent")
        assert_rounded(100 * row.met_freq, metastatic_percent, 1, f"Table 3 {gene} metastatic percent")
        assert_rounded(100 * row.freq_diff, difference_percent, 1, f"Table 3 {gene} frequency difference")
        assert_rounded(row.mcnemar_p, p_value, p_digits, f"Table 3 {gene} McNemar p")
        assert_rounded(row.mcnemar_fdr, fdr, 3, f"Table 3 {gene} FDR")
        assert_rounded(row.tcga_mutated_n, tcga_count, 0, f"Table 3 {gene} TCGA count")
        assert_rounded(row.tcga_total_n, 1083, 0, f"Table 3 {gene} TCGA denominator")
        assert_rounded(row.tcga_mutated_percent, tcga_percent, 1, f"Table 3 {gene} TCGA percent")


def validate_table4():
    observed = pd.read_csv(TABLES_DIR / "table4_aurora_burden_summary.csv").set_index("statistic")
    continuous = {
        "Mean": (0.95, 1.87, 2),
        "Median": (1, 2, 0),
        "Minimum": (0, 0, 0),
        "Maximum": (3, 5, 0),
    }
    for label, (primary, metastasis, digits) in continuous.items():
        assert_rounded(observed.loc[label, "primary"], primary, digits, f"Table 4 {label} primary")
        assert_rounded(observed.loc[label, "metastasis"], metastasis, digits, f"Table 4 {label} metastasis")
    categories = {
        "Burden = 0": (13, 33.3, 5, 12.8),
        "Burden = 1": (17, 43.6, 12, 30.8),
        "Burden = 2": (7, 17.9, 11, 28.2),
        "Burden >= 3": (2, 5.1, 11, 28.2),
    }
    for label, (primary_n, primary_percent, metastasis_n, metastasis_percent) in categories.items():
        assert_rounded(observed.loc[label, "primary"], primary_n, 0, f"Table 4 {label} primary count")
        assert_rounded(observed.loc[label, "primary_percent"], primary_percent, 1, f"Table 4 {label} primary percent")
        assert_rounded(observed.loc[label, "metastasis"], metastasis_n, 0, f"Table 4 {label} metastatic count")
        assert_rounded(observed.loc[label, "metastasis_percent"], metastasis_percent, 1, f"Table 4 {label} metastatic percent")
    assert_rounded(observed.loc["Wilcoxon signed-rank p", "metastasis"], 0.001, 3, "Table 4 Wilcoxon p")


def validate_supplementary_table1():
    observed = pd.read_csv(TABLES_DIR / "supplementary_table1_tcga_summary.csv").set_index("characteristic")
    count_rows = {
        "Total patients": (1083, None),
        "Overall survival events": (151, 13.9),
        "Stage I": (183, 16.9), "Stage II": (611, 56.4), "Stage III": (245, 22.6),
        "Stage IV": (20, 1.8), "Stage missing": (24, 2.2),
        "Subtype HR+/HER2-": (439, 40.5), "Subtype HR+/HER2+": (122, 11.3),
        "Subtype HER2+": (38, 3.5), "Subtype Triple Negative": (119, 11.0),
        "Subtype Unknown": (365, 33.7), "Samples with TMB": (964, 89.0),
    }
    for label, (count, percent) in count_rows.items():
        assert_rounded(observed.loc[label, "value"], count, 0, f"Supplementary Table 1 {label} count")
        if percent is not None:
            assert_rounded(observed.loc[label, "percent"], percent, 1, f"Supplementary Table 1 {label} percent")
    continuous = {
        "Median follow-up (years)": (2.4, 1), "Age median": (58, 0), "Age minimum": (26, 0), "Age maximum": (90, 0),
        "TMB median": (1.07, 2), "TMB mean": (2.18, 2), "TMB minimum": (0.03, 2), "TMB maximum": (147.10, 2),
        "TMB Q1": (0.7, 1), "TMB Q3": (1.9, 1),
    }
    for label, (value, digits) in continuous.items():
        assert_rounded(observed.loc[label, "value"], value, digits, f"Supplementary Table 1 {label}")


def validate_supplementary_table2():
    observed = pd.read_csv(TABLES_DIR / "supplementary_table2_panel_distribution.csv").set_index("measure")
    expected = {
        "Panel burden 0": (500, 46.2), "Panel burden 1": (455, 42.0), "Panel burden 2": (112, 10.3),
        "Panel burden 3": (12, 1.1), "Panel burden 4": (3, 0.3), "Panel burden 6": (1, 0.1),
        "Panel wild-type": (500, 46.2), "Panel mutated": (583, 53.8),
    }
    for label, (count, percent) in expected.items():
        assert_rounded(observed.loc[label, "value"], count, 0, f"Supplementary Table 2 {label} count")
        assert_rounded(observed.loc[label, "percent"], percent, 1, f"Supplementary Table 2 {label} percent")
    assert_rounded(observed.loc["Spearman rho", "value"], 0.321, 3, "Supplementary Table 2 Spearman rho")
    if float(observed.loc["Spearman p", "value"]) >= 0.0001:
        raise AssertionError("Supplementary Table 2 Spearman p is not <0.0001")


def validate_figures():
    for stem in ("fig1_panel_construction", "fig2_burden_boxplot", "fig3_tcga_distribution_tmb"):
        for suffix in ("png", "svg"):
            path = FIGURES_DIR / f"{stem}.{suffix}"
            if not path.exists() or path.stat().st_size == 0:
                raise AssertionError(f"Missing figure output: {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-figures", action="store_true")
    args = parser.parse_args()
    validate_table1()
    validate_table2()
    validate_table3()
    validate_table4()
    validate_supplementary_table1()
    validate_supplementary_table2()
    if args.require_figures:
        validate_figures()
    print("All displayed Python manuscript values passed validation.")


if __name__ == "__main__":
    main()
