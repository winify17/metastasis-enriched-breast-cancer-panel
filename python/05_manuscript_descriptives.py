"""Packages: pandas, numpy, scipy, matplotlib. Prerequisites: archived AURORA aggregates and cleaned TCGA input."""

import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from common import AURORA_DIR, FIGURES_DIR, TABLES_DIR, TCGA_PROCESSED_DIR, ensure_output_dirs


BLUE = "#1683E6"
ORANGE = "#FF9000"
PINK = "#F58ADC"
GREEN = "#67A357"


def save_figure(figure, stem):
    figure.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    figure.savefig(FIGURES_DIR / f"{stem}.svg", bbox_inches="tight")
    figure.clear()


def get_pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot
    return pyplot


def expand_burden_pairs(counts):
    return counts.loc[counts.index.repeat(counts.patient_count)].reset_index(drop=True)


def write_aurora_tables(cohort, mcnemar, panel, selection, burden, tcga_n):
    cohort.to_csv(TABLES_DIR / "table1_aurora_cohort_summary.csv", index=False)

    valid = mcnemar.mcnemar_p.notna()
    most_significant = mcnemar.loc[valid & mcnemar.direction.eq("met_enriched")].sort_values(["mcnemar_p", "gene"]).iloc[0]
    summary = pd.DataFrame(
        [
            ("Genes tested", len(mcnemar)),
            ("Genes with valid p-value", int(valid.sum())),
            ("Met-enriched genes (FDR < 0.05)", int(((mcnemar.mcnemar_fdr < 0.05) & mcnemar.direction.eq("met_enriched")).sum())),
            ("Met-depleted genes (FDR < 0.05)", int(((mcnemar.mcnemar_fdr < 0.05) & mcnemar.direction.eq("met_depleted")).sum())),
            ("Most significant gene", most_significant.gene),
            ("Most significant McNemar p", most_significant.mcnemar_p),
            ("Most significant FDR", most_significant.mcnemar_fdr),
        ],
        columns=["parameter", "value"],
    )
    summary.to_csv(TABLES_DIR / "table2_mcnemar_summary.csv", index=False)

    panel_rows = mcnemar.set_index("gene").loc[panel].reset_index()
    tcga_counts = selection.set_index("gene").tcga_mutated_n
    panel_rows["tcga_mutated_n"] = panel_rows.gene.map(tcga_counts).astype(int)
    panel_rows["tcga_total_n"] = tcga_n
    panel_rows["tcga_mutated_percent"] = 100 * panel_rows.tcga_mutated_n / tcga_n
    panel_rows[["gene", "n_prim", "n_met", "prim_freq", "met_freq", "freq_diff", "mcnemar_p", "mcnemar_fdr", "tcga_mutated_n", "tcga_total_n", "tcga_mutated_percent"]].to_csv(
        TABLES_DIR / "table3_panel_composition.csv", index=False
    )

    primary = burden.primary_burden
    metastasis = burden.metastatic_burden
    _, p_value = wilcoxon(metastasis, primary)
    burden_summary = pd.DataFrame(
        [
            ("Mean", primary.mean(), metastasis.mean()),
            ("Median", primary.median(), metastasis.median()),
            ("Minimum", primary.min(), metastasis.min()),
            ("Maximum", primary.max(), metastasis.max()),
            ("Burden = 0", int(primary.eq(0).sum()), int(metastasis.eq(0).sum())),
            ("Burden = 1", int(primary.eq(1).sum()), int(metastasis.eq(1).sum())),
            ("Burden = 2", int(primary.eq(2).sum()), int(metastasis.eq(2).sum())),
            ("Burden >= 3", int(primary.ge(3).sum()), int(metastasis.ge(3).sum())),
            ("Wilcoxon signed-rank p", np.nan, p_value),
        ],
        columns=["statistic", "primary", "metastasis"],
    )
    burden_summary["primary_percent"] = np.nan
    burden_summary["metastasis_percent"] = np.nan
    category_rows = burden_summary.statistic.isin(["Burden = 0", "Burden = 1", "Burden = 2", "Burden >= 3"])
    burden_summary.loc[category_rows, "primary_percent"] = 100 * burden_summary.loc[category_rows, "primary"] / len(burden)
    burden_summary.loc[category_rows, "metastasis_percent"] = 100 * burden_summary.loc[category_rows, "metastasis"] / len(burden)
    burden_summary.to_csv(TABLES_DIR / "table4_aurora_burden_summary.csv", index=False)


def write_tcga_tables(tcga):
    stage_counts = tcga.stage.value_counts()
    subtype_counts = tcga.mol_subtype.fillna("Unknown").value_counts()
    tmb = tcga.tmb.dropna()
    supplementary_1 = pd.DataFrame(
        [
            ("Total patients", len(tcga)),
            ("Overall survival events", int(tcga.event.sum())),
            ("Median follow-up (years)", tcga.time_years.median()),
            ("Age median", tcga.age.median()),
            ("Age minimum", tcga.age.min()),
            ("Age maximum", tcga.age.max()),
            ("Stage I", stage_counts.get("I", 0)),
            ("Stage II", stage_counts.get("II", 0)),
            ("Stage III", stage_counts.get("III", 0)),
            ("Stage IV", stage_counts.get("IV", 0)),
            ("Stage missing", tcga.stage.isna().sum()),
            ("Subtype HR+/HER2-", subtype_counts.get("HR+/HER2-", 0)),
            ("Subtype HR+/HER2+", subtype_counts.get("HR+/HER2+", 0)),
            ("Subtype HER2+", subtype_counts.get("HER2+", 0)),
            ("Subtype Triple Negative", subtype_counts.get("Triple Negative", 0)),
            ("Subtype Unknown", subtype_counts.get("Unknown", 0)),
            ("Samples with TMB", len(tmb)),
            ("TMB median", tmb.median()),
            ("TMB mean", tmb.mean()),
            ("TMB minimum", tmb.min()),
            ("TMB maximum", tmb.max()),
            ("TMB Q1", tmb.quantile(0.25)),
            ("TMB Q3", tmb.quantile(0.75)),
        ],
        columns=["characteristic", "value"],
    )
    supplementary_1["percent"] = np.nan
    percentage_labels = [
        "Overall survival events", "Stage I", "Stage II", "Stage III", "Stage IV", "Stage missing",
        "Subtype HR+/HER2-", "Subtype HR+/HER2+", "Subtype HER2+", "Subtype Triple Negative", "Subtype Unknown",
        "Samples with TMB",
    ]
    supplementary_1.loc[supplementary_1.characteristic.isin(percentage_labels), "percent"] = (
        100 * supplementary_1.loc[supplementary_1.characteristic.isin(percentage_labels), "value"] / len(tcga)
    )
    supplementary_1.to_csv(TABLES_DIR / "supplementary_table1_tcga_summary.csv", index=False)

    burden_counts = tcga.panel_burden.value_counts().sort_index()
    rho, p_value = spearmanr(tcga.loc[tcga.tmb.notna(), "panel_burden"], tcga.loc[tcga.tmb.notna(), "tmb"])
    supplementary_2 = pd.DataFrame(
        [(f"Panel burden {int(score)}", int(count)) for score, count in burden_counts.items()]
        + [
            ("Panel wild-type", int(tcga.panel_binary.eq(0).sum())),
            ("Panel mutated", int(tcga.panel_binary.eq(1).sum())),
            ("Spearman rho", rho),
            ("Spearman p", p_value),
        ],
        columns=["measure", "value"],
    )
    supplementary_2["percent"] = np.nan
    count_rows = supplementary_2.measure.str.startswith("Panel burden") | supplementary_2.measure.isin(["Panel wild-type", "Panel mutated"])
    supplementary_2.loc[count_rows, "percent"] = 100 * supplementary_2.loc[count_rows, "value"] / len(tcga)
    supplementary_2.to_csv(TABLES_DIR / "supplementary_table2_panel_distribution.csv", index=False)


def figure_1(mcnemar, panel):
    plt = get_pyplot()
    rows = mcnemar.set_index("gene").loc[panel].reset_index()
    x = np.arange(len(rows))
    width = 0.36
    figure, axis = plt.subplots(figsize=(10, 6))
    axis.bar(x - width / 2, 100 * rows.prim_freq, width, label="Primary", color=BLUE)
    axis.bar(x + width / 2, 100 * rows.met_freq, width, label="Metastasis", color=ORANGE)
    for index, row in rows.iterrows():
        height = 100 * max(row.prim_freq, row.met_freq)
        axis.text(index, height + 2, f"p={row.mcnemar_p:.3g}", ha="center", fontsize=8)
    axis.set_xticks(x, rows.gene, rotation=30, ha="right")
    axis.set_ylabel("Mutation-positive patients (%)")
    axis.set_title("AURORA paired primary and metastatic mutation frequencies")
    axis.legend(frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    save_figure(figure, "fig1_panel_construction")


def figure_2(burden):
    plt = get_pyplot()
    figure, axis = plt.subplots(figsize=(7, 6))
    for _, row in burden.iterrows():
        difference = row.metastatic_burden - row.primary_burden
        color = GREEN if difference > 0 else (PINK if difference < 0 else "#A0A0A0")
        axis.plot([0, 1], [row.primary_burden, row.metastatic_burden], color=color, alpha=0.45, linewidth=1)
    axis.boxplot(
        [burden.primary_burden, burden.metastatic_burden],
        positions=[0, 1], widths=0.35, patch_artist=True,
        boxprops={"facecolor": "white", "edgecolor": "black"},
        medianprops={"color": "black", "linewidth": 2},
    )
    _, p_value = wilcoxon(burden.metastatic_burden, burden.primary_burden)
    axis.set_xticks([0, 1], ["Primary", "Metastasis"])
    axis.set_ylabel("Seven-gene panel burden")
    axis.set_title("AURORA paired panel burden")
    axis.text(0.5, 0.96, f"Wilcoxon signed-rank p={p_value:.4f}", transform=axis.transAxes, ha="center", va="top")
    axis.spines[["top", "right"]].set_visible(False)
    save_figure(figure, "fig2_burden_boxplot")


def figure_3(tcga):
    plt = get_pyplot()
    available = tcga[tcga.tmb.notna()].copy()
    rho, p_value = spearmanr(available.panel_burden, available.tmb)
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    counts = tcga.panel_burden.value_counts().sort_index()
    axes[0].bar(counts.index.astype(str), counts.values, color=BLUE)
    axes[0].set_xlabel("Panel burden")
    axes[0].set_ylabel("Patients")
    axes[0].set_title("Panel burden distribution")
    rng = np.random.default_rng(2024)
    jitter = rng.normal(0, 0.055, len(available))
    axes[1].scatter(available.panel_burden + jitter, available.tmb, s=12, alpha=0.45, color=ORANGE, edgecolors="none")
    medians = available.groupby("panel_burden").tmb.median()
    axes[1].plot(medians.index, medians.values, color=BLUE, marker="o", linewidth=2, label="Median TMB")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Panel burden")
    axes[1].set_ylabel("TMB (mutations/Mb, log scale)")
    axes[1].set_title("Panel burden and tumour mutational burden")
    axes[1].text(0.04, 0.96, f"Spearman rho={rho:.3f}\np={p_value:.2g}", transform=axes[1].transAxes, va="top")
    axes[1].legend(frameon=False)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    save_figure(figure, "fig3_tcga_distribution_tmb")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-only", action="store_true")
    args = parser.parse_args()
    ensure_output_dirs()
    cohort = pd.read_csv(AURORA_DIR / "aurora_cohort_summary.csv")
    mcnemar = pd.read_csv(AURORA_DIR / "aurora_mcnemar_results.csv")
    panel = pd.read_csv(AURORA_DIR / "final_panel_genes.csv").gene.tolist()
    selection = pd.read_csv(AURORA_DIR / "panel_selection.csv")
    burden_counts = pd.read_csv(AURORA_DIR / "aurora_paired_burden_counts.csv")
    burden = expand_burden_pairs(burden_counts)
    tcga = pd.read_csv(TCGA_PROCESSED_DIR / "tcga_brca_panel_survival_input.csv")

    write_aurora_tables(cohort, mcnemar, panel, selection, burden, len(tcga))
    write_tcga_tables(tcga)
    if not args.tables_only:
        figure_1(mcnemar, panel)
        figure_2(burden)
        figure_3(tcga)
        print("Generated Tables 1-4, Supplementary Tables 1-2, and Figures 1-3.")
    else:
        print("Generated Tables 1-4 and Supplementary Tables 1-2.")


if __name__ == "__main__":
    main()
