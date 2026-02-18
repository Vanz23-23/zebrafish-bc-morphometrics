"""
Morphometric Analysis of Zebrafish Retinal Bipolar Cells
=========================================================
Standalone analysis script for reproducibility. Loads pre-processed
morphometric data and generates statistical comparisons and figures.

Author: Vanz Labitad, University of Sussex
Dissertation: Regional Morphological Specialisation of Monostratifying
              Bipolar Cells in the Larval Zebrafish Retina (2025)
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATA_PATH = Path("data/morphometrics_clean.csv")
OUTPUT_DIR = Path("figures")

METRICS = ["lateral_spread_um", "ipl_depth_pct", "volume_um3", "surface_area_um2"]

METRIC_LABELS = {
    "lateral_spread_um": "Lateral Spread (um)",
    "ipl_depth_pct": "IPL Depth (%)",
    "volume_um3": "Volume (um3)",
    "surface_area_um2": "Surface Area (um2)",
}

# Colorblind-friendly palette (Tol bright scheme)
COLORS = {
    "S2": "#EE6677",   # rose — OFF-pathway
    "S4": "#4477AA",   # blue — ON-pathway
    "Dorsal": "#228833",
    "Ventral": "#CCBB44",
}

# Publication-quality defaults
plt.rcParams.update({
    "font.size": 11,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def load_data(path):
    """Load the cleaned morphometrics dataset."""
    df = pd.read_csv(path)
    print(f"Loaded {len(df)} cells from {path}")
    return df


def print_dataset_summary(df):
    """Print cell counts per annotator x region x cell_type."""
    print("\n=== DATASET SUMMARY ===")
    summary = (
        df.groupby(["annotator", "region", "cell_type"])
        .size()
        .reset_index(name="n")
    )
    print(summary.to_string(index=False))
    print(f"\nTotal: {len(df)} cells")
    print(f"  S2: {len(df[df['cell_type'] == 'S2'])}")
    print(f"  S4: {len(df[df['cell_type'] == 'S4'])}")
    print(f"  Dorsal: {len(df[df['region'] == 'Dorsal'])}")
    print(f"  Ventral: {len(df[df['region'] == 'Ventral'])}")


def cohens_d(a, b):
    """
    Compute Cohen's d using pooled standard deviation.

    d = (mean_a - mean_b) / SD_pooled
    where SD_pooled = sqrt(((n_a-1)*s_a^2 + (n_b-1)*s_b^2) / (n_a+n_b-2))
    """
    na, nb = len(a), len(b)
    sa, sb = np.std(a, ddof=1), np.std(b, ddof=1)
    pooled_sd = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    if pooled_sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled_sd


def run_comparisons(df):
    """
    Run Mann-Whitney U tests and compute Cohen's d for three comparison sets:
      1. S2 vs S4 (all cells) — validates cell type classification
      2. Dorsal vs Ventral (all cells) — broad regional differences
      3. Dorsal S4 vs Ventral S4 — key regional comparison for ON-pathway
    """
    comparisons = [
        ("S2 vs S4", df[df["cell_type"] == "S2"], df[df["cell_type"] == "S4"]),
        ("Dorsal vs Ventral", df[df["region"] == "Dorsal"], df[df["region"] == "Ventral"]),
        (
            "Dorsal S4 vs Ventral S4",
            df[(df["region"] == "Dorsal") & (df["cell_type"] == "S4")],
            df[(df["region"] == "Ventral") & (df["cell_type"] == "S4")],
        ),
    ]

    rows = []
    for label, group_a, group_b in comparisons:
        for metric in METRICS:
            a = group_a[metric].dropna().values
            b = group_b[metric].dropna().values

            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            d = cohens_d(a, b)

            rows.append({
                "comparison": label,
                "metric": metric,
                "n_a": len(a),
                "n_b": len(b),
                "mean_a": np.mean(a),
                "mean_b": np.mean(b),
                "U_statistic": stat,
                "p_value": p,
                "cohens_d": d,
                "significant": p < 0.05,
            })

    results = pd.DataFrame(rows)
    return results


def print_results(results):
    """Print statistical results in a clean tabular format."""
    print("\n=== STATISTICAL RESULTS ===")
    print(f"{'Comparison':<28} {'Metric':<22} {'n_a':>4} {'n_b':>4} "
          f"{'U':>10} {'p':>12} {'d':>8} {'Sig':>5}")
    print("-" * 100)

    for _, row in results.iterrows():
        sig_marker = "***" if row["p_value"] < 0.001 else (
            "**" if row["p_value"] < 0.01 else (
            "*" if row["p_value"] < 0.05 else "ns"))
        print(
            f"{row['comparison']:<28} {row['metric']:<22} "
            f"{row['n_a']:>4} {row['n_b']:>4} "
            f"{row['U_statistic']:>10.1f} {row['p_value']:>12.2e} "
            f"{row['cohens_d']:>8.3f} {sig_marker:>5}"
        )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_ipl_depth_by_type(df):
    """
    Figure: IPL depth distribution by cell type (S2 vs S4).

    S2 cells (OFF-pathway) should stratify in the upper IPL (~10-45%),
    while S4 cells (ON-pathway) stratify deeper (~45-80%).
    """
    fig, ax = plt.subplots(figsize=(6, 5))

    # Violin plot with embedded box
    parts = ax.violinplot(
        [df[df["cell_type"] == "S2"]["ipl_depth_pct"].dropna(),
         df[df["cell_type"] == "S4"]["ipl_depth_pct"].dropna()],
        positions=[0, 1],
        showmeans=False, showmedians=False, showextrema=False,
    )

    for i, (body, ct) in enumerate(zip(parts["bodies"], ["S2", "S4"])):
        body.set_facecolor(COLORS[ct])
        body.set_alpha(0.4)

    # Overlay box plots
    bp = ax.boxplot(
        [df[df["cell_type"] == "S2"]["ipl_depth_pct"].dropna(),
         df[df["cell_type"] == "S4"]["ipl_depth_pct"].dropna()],
        positions=[0, 1],
        widths=0.15, patch_artist=True,
        showfliers=False, zorder=3,
    )

    for i, (box, ct) in enumerate(zip(bp["boxes"], ["S2", "S4"])):
        box.set_facecolor(COLORS[ct])
        box.set_alpha(0.7)

    for element in ["whiskers", "caps", "medians"]:
        for line in bp[element]:
            line.set_color("black")
            line.set_linewidth(1.2)

    # Overlay individual points
    for i, ct in enumerate(["S2", "S4"]):
        vals = df[df["cell_type"] == ct]["ipl_depth_pct"].dropna()
        jitter = np.random.default_rng(42).uniform(-0.08, 0.08, len(vals))
        ax.scatter(i + jitter, vals, color=COLORS[ct], alpha=0.5, s=15,
                   edgecolors="white", linewidth=0.3, zorder=4)

    n_s2 = len(df[df["cell_type"] == "S2"])
    n_s4 = len(df[df["cell_type"] == "S4"])
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"S2 (OFF)\nn={n_s2}", f"S4 (ON)\nn={n_s4}"],
                       fontweight="bold")
    ax.set_ylabel("IPL Depth (%)")
    ax.set_title("IPL Stratification by Cell Type")

    # Reference zones: expected IPL strata
    ax.axhspan(10, 45, color=COLORS["S2"], alpha=0.06, label="Expected S2 range")
    ax.axhspan(45, 80, color=COLORS["S4"], alpha=0.06, label="Expected S4 range")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    plt.tight_layout()
    out = OUTPUT_DIR / "fig_ipl_depth_by_type.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_lateral_spread_regional(df, results):
    """
    Figure: Lateral spread of Dorsal S4 vs Ventral S4 bipolar cells.

    This is the key regional comparison. Ventral S4 cells are expected to
    have wider terminal arbors, reflecting adaptation to the bright
    dorsal visual field (sky).
    """
    dorsal_s4 = df[(df["region"] == "Dorsal") & (df["cell_type"] == "S4")]
    ventral_s4 = df[(df["region"] == "Ventral") & (df["cell_type"] == "S4")]

    # Get stats for annotation
    row = results[
        (results["comparison"] == "Dorsal S4 vs Ventral S4")
        & (results["metric"] == "lateral_spread_um")
    ].iloc[0]

    fig, ax = plt.subplots(figsize=(5.5, 5))

    # Box plot
    bp = ax.boxplot(
        [dorsal_s4["lateral_spread_um"].dropna(),
         ventral_s4["lateral_spread_um"].dropna()],
        positions=[0, 1],
        widths=0.4, patch_artist=True,
        showfliers=False, zorder=2,
    )

    for box, color in zip(bp["boxes"], [COLORS["Dorsal"], COLORS["Ventral"]]):
        box.set_facecolor(color)
        box.set_alpha(0.5)
    for element in ["whiskers", "caps", "medians"]:
        for line in bp[element]:
            line.set_color("black")
            line.set_linewidth(1.2)

    # Strip plot overlay
    for i, (subset, color) in enumerate([
        (dorsal_s4, COLORS["Dorsal"]),
        (ventral_s4, COLORS["Ventral"]),
    ]):
        vals = subset["lateral_spread_um"].dropna()
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(i + jitter, vals, color=color, alpha=0.6, s=20,
                   edgecolors="white", linewidth=0.3, zorder=3)

    n_d = len(dorsal_s4)
    n_v = len(ventral_s4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Dorsal S4\nn={n_d}", f"Ventral S4\nn={n_v}"],
                       fontweight="bold")
    ax.set_ylabel("Lateral Spread (um)")
    ax.set_title("Terminal Arbor Width: Dorsal vs Ventral S4")

    # Significance annotation
    y_max = max(dorsal_s4["lateral_spread_um"].max(),
                ventral_s4["lateral_spread_um"].max())
    bar_y = y_max * 1.05
    ax.plot([0, 0, 1, 1], [bar_y, bar_y * 1.02, bar_y * 1.02, bar_y],
            color="black", linewidth=1.2)

    p_str = "p < 0.001" if row["p_value"] < 0.001 else f"p = {row['p_value']:.3f}"
    ax.text(0.5, bar_y * 1.04, f"{p_str}\nCohen's d = {row['cohens_d']:.2f}",
            ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_ylim(top=bar_y * 1.15)
    plt.tight_layout()
    out = OUTPUT_DIR / "fig_lateral_spread_regional.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_effect_heatmap(results):
    """
    Figure: Heatmap of Cohen's d effect sizes across all comparisons and metrics.

    Colour scale: blue = group B larger, red = group A larger.
    Significance markers overlaid as asterisks.
    """
    pivot = results.pivot(index="metric", columns="comparison", values="cohens_d")

    # Reorder for readability
    col_order = ["S2 vs S4", "Dorsal vs Ventral", "Dorsal S4 vs Ventral S4"]
    row_order = ["lateral_spread_um", "ipl_depth_pct", "volume_um3", "surface_area_um2"]
    pivot = pivot.reindex(index=row_order, columns=col_order)

    # Readable row labels
    pivot.index = [METRIC_LABELS.get(m, m) for m in pivot.index]

    fig, ax = plt.subplots(figsize=(7, 4))
    sns.heatmap(
        pivot,
        annot=True, fmt=".2f",
        cmap="RdBu_r", center=0,
        vmin=-5, vmax=5,
        linewidths=1, linecolor="white",
        cbar_kws={"label": "Cohen's d", "shrink": 0.8},
        ax=ax,
    )

    # Overlay significance asterisks
    for i, metric in enumerate(row_order):
        for j, comp in enumerate(col_order):
            row = results[
                (results["metric"] == metric) & (results["comparison"] == comp)
            ]
            if len(row) > 0 and row.iloc[0]["significant"]:
                p = row.iloc[0]["p_value"]
                stars = "***" if p < 0.001 else ("**" if p < 0.01 else "*")
                ax.text(j + 0.5, i + 0.78, stars,
                        ha="center", va="center", fontsize=10,
                        fontweight="bold", color="black")

    ax.set_title("Effect Sizes Across Comparisons")
    ax.set_ylabel("")
    ax.set_xlabel("")
    plt.tight_layout()
    out = OUTPUT_DIR / "fig_effect_heatmap.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Run the full analysis: load data, compute stats, generate figures."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_data(DATA_PATH)
    print_dataset_summary(df)

    results = run_comparisons(df)
    print_results(results)

    plot_ipl_depth_by_type(df)
    plot_lateral_spread_regional(df, results)
    plot_effect_heatmap(results)

    print("\nDone. Figures saved to figures/")


if __name__ == "__main__":
    main()
