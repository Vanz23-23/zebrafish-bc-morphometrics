"""
Morphometric Analysis of Zebrafish Retinal Bipolar Cells
=========================================================
Standalone analysis script for reproducibility. Loads pre-processed
morphometric data and generates statistical comparisons and figures.

Author: Vanz Labitad, University of Sussex
Dissertation: Regional Morphological Specialisation of Monostratifying
              Bipolar Cells in the Larval Zebrafish Retina (2025)
"""

import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_HERE = Path(__file__).parent.resolve()
DATA_PATH = _HERE.parent / "data" / "morphometrics_clean.csv"
OUTPUT_DIR = _HERE.parent / "figures"

METRICS = ["lateral_spread_um", "ipl_depth_pct", "volume_um3", "surface_area_um2"]

METRIC_LABELS = {
    "lateral_spread_um": "Lateral Spread (um)",
    "ipl_depth_pct": "IPL Depth (%)",
    "volume_um3": "Volume (um3)",
    "surface_area_um2": "Surface Area (um2)",
}


COLORS = {
    "S2": "#EE6677",   # rose — OFF-pathway
    "S4": "#4477AA",   # blue — ON-pathway
    "Dorsal": "#228833",
    "Ventral": "#CCBB44",
}

# Significance threshold — uncorrected; see run_comparisons docstring for justification.
ALPHA = 0.05

# Minimum group size for statistical tests
MIN_N = 3

# Required columns in the input dataset
REQUIRED_COLUMNS = {"cell_type", "region", "annotator"} | set(METRICS)


def configure_plot_style() -> None:
    """Apply publication-quality matplotlib defaults. Call from main()."""
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


def log_package_versions() -> None:
    """Print Python and package versions for reproducibility."""
    import scipy
    print(f"\nPython: {sys.version}")
    print(f"numpy: {np.__version__}, pandas: {pd.__version__}, "
          f"scipy: {scipy.__version__}, matplotlib: {plt.matplotlib.__version__}, "
          f"seaborn: {sns.__version__}")


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def load_data(path: Path) -> pd.DataFrame:
    """Load the cleaned morphometrics dataset."""
    if not path.exists():
        raise FileNotFoundError(
            f"Data file not found: {path.resolve()}\n"
            "Run this script from the project root, or check DATA_PATH."
        )
    df = pd.read_csv(path)
    print(f"Loaded {len(df)} cells from {path.resolve()}")
    return df


def validate_data(df: pd.DataFrame) -> None:
    """Check that all required columns are present in the dataset."""
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def print_dataset_summary(df: pd.DataFrame) -> None:
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


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """
    Compute Cohen's d using pooled standard deviation.

    d = (mean_a - mean_b) / SD_pooled
    where SD_pooled = sqrt(((n_a-1)*s_a^2 + (n_b-1)*s_b^2) / (n_a+n_b-2))

    NOTE: Cohen's d is a parametric effect size and is technically mismatched with
    Mann-Whitney U (a non-parametric test). We report it here as a *descriptive*
    standardised mean difference because: (1) it is widely understood and directly
    comparable across metrics, (2) our sample sizes (~20-60 per group) are large
    enough that the mean/SD are meaningful summaries even when normality is violated,
    and (3) rank-biserial r is reported alongside for non-parametric completeness.
    Positive values indicate group a > group b.
    """
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan  # Cohen's d undefined for n < 2
    sa, sb = np.std(a, ddof=1), np.std(b, ddof=1)
    pooled_sd = np.sqrt(((na - 1) * sa**2 + (nb - 1) * sb**2) / (na + nb - 2))
    if pooled_sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / pooled_sd


def rank_biserial_r(u_stat: float, n_a: int, n_b: int) -> float:
    """
    Rank-biserial correlation: non-parametric effect size for Mann-Whitney U.
    r = 1 - (2U / (n_a * n_b)), bounded [-1, 1].
    """
    if n_a * n_b == 0:
        return np.nan
    return 1 - (2 * u_stat) / (n_a * n_b)


def run_comparisons(df: pd.DataFrame) -> pd.DataFrame:
    """
    Run Mann-Whitney U tests and compute effect sizes for three comparison sets:
      1. S2 vs S4 (all cells) — validates cell type classification
      2. Dorsal vs Ventral (all cells) — broad regional differences
      3. Dorsal S4 vs Ventral S4 — key regional comparison for ON-pathway

    NOTE on multiple comparisons: 3 comparisons × 4 metrics = 12 tests
    at α = 0.05. No family-wise correction (Bonferroni, Holm) or false discovery
    rate control (Benjamini-Hochberg) is applied. Justification:

      1. All comparisons are pre-specified and hypothesis-driven, not exploratory.
      2. The three comparison families test biologically motivated contrasts
         (cell type, region, and their interaction for S4 cells).
      3. With 12 tests at α = 0.05, ~0.6 false positives are expected by chance.
         Readers should weight results by effect size (Cohen's d, rank-biserial r)
         rather than relying on p-value thresholds alone.
      4. Applying Bonferroni (α/12 = 0.004) would be overly conservative given
         the correlated nature of morphometric variables.

    The 'significant' column uses uncorrected p < ALPHA. This decision and its
    rationale must be stated explicitly in the dissertation methods section.
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

            # Guard against empty or too-small groups
            if len(a) < MIN_N or len(b) < MIN_N:
                print(f"WARNING: Skipping {label} / {metric} — "
                      f"insufficient n ({len(a)}, {len(b)})")
                rows.append({
                    "comparison": label, "metric": metric,
                    "n_a": len(a), "n_b": len(b),
                    "mean_a": np.nan, "mean_b": np.nan,
                    "U_statistic": np.nan, "p_value": np.nan,
                    "cohens_d": np.nan, "rank_biserial_r": np.nan,
                    "significant": False,
                })
                continue

            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            d = cohens_d(a, b)
            r = rank_biserial_r(stat, len(a), len(b))

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
                "rank_biserial_r": r,
                "significant": p < ALPHA,
            })

    return pd.DataFrame(rows)


def print_results(results: pd.DataFrame) -> None:
    """Print statistical results in a clean tabular format."""
    print("\n=== STATISTICAL RESULTS ===")
    print(f"{'Comparison':<28} {'Metric':<22} {'n_a':>4} {'n_b':>4} "
          f"{'U':>10} {'p':>12} {'d':>8} {'r':>8} {'Sig':>5}")
    print("-" * 108)

    for _, row in results.iterrows():
        sig_marker = "***" if row["p_value"] < 0.001 else (
            "**" if row["p_value"] < 0.01 else (
            "*" if row["p_value"] < ALPHA else "ns"))
        print(
            f"{row['comparison']:<28} {row['metric']:<22} "
            f"{row['n_a']:>4} {row['n_b']:>4} "
            f"{row['U_statistic']:>10.1f} {row['p_value']:>12.2e} "
            f"{row['cohens_d']:>8.3f} {row['rank_biserial_r']:>8.3f} {sig_marker:>5}"
        )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def plot_ipl_depth_by_type(df: pd.DataFrame) -> None:
    """
    Figure: IPL depth distribution by cell type (S2 vs S4).

    S2 cells (OFF-pathway) should stratify in the upper IPL (~10-45%),
    while S4 cells (ON-pathway) stratify deeper (~45-80%).
    """
    fig, ax = plt.subplots(figsize=(6, 5))

    vals_s2 = df[df["cell_type"] == "S2"]["ipl_depth_pct"].dropna()
    vals_s4 = df[df["cell_type"] == "S4"]["ipl_depth_pct"].dropna()

    # Violin plot with embedded box
    parts = ax.violinplot(
        [vals_s2, vals_s4],
        positions=[0, 1],
        showmeans=False, showmedians=False, showextrema=False,
    )

    for i, (body, ct) in enumerate(zip(parts["bodies"], ["S2", "S4"])):
        body.set_facecolor(COLORS[ct])
        body.set_alpha(0.4)

    # Overlay box plots
    bp = ax.boxplot(
        [vals_s2, vals_s4],
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

    # Overlay individual points (single RNG instance, advances across groups)
    rng = np.random.default_rng(42)
    for i, (ct, vals) in enumerate([("S2", vals_s2), ("S4", vals_s4)]):
        jitter = rng.uniform(-0.08, 0.08, len(vals))
        ax.scatter(i + jitter, vals, color=COLORS[ct], alpha=0.5, s=15,
                   edgecolors="white", linewidth=0.3, zorder=4)

    # n computed from dropna'd series (matches plotted data)
    n_s2, n_s4 = len(vals_s2), len(vals_s4)
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


def plot_lateral_spread_regional(df: pd.DataFrame, results: pd.DataFrame) -> None:
    """
    Figure: Lateral spread of Dorsal S4 vs Ventral S4 bipolar cells.

    This is the key regional comparison. Ventral S4 cells are expected to
    have wider terminal arbors, reflecting adaptation to the bright
    dorsal visual field (sky).
    """
    dorsal_s4 = df[(df["region"] == "Dorsal") & (df["cell_type"] == "S4")]
    ventral_s4 = df[(df["region"] == "Ventral") & (df["cell_type"] == "S4")]

    # Get stats for annotation (guard against missing comparison)
    stat_rows = results[
        (results["comparison"] == "Dorsal S4 vs Ventral S4")
        & (results["metric"] == "lateral_spread_um")
    ]
    if stat_rows.empty:
        print("WARNING: No stats found for Dorsal S4 vs Ventral S4 lateral_spread_um — skipping figure")
        return
    row = stat_rows.iloc[0]

    fig, ax = plt.subplots(figsize=(5.5, 5))

    d_vals = dorsal_s4["lateral_spread_um"].dropna()
    v_vals = ventral_s4["lateral_spread_um"].dropna()

    # Box plot
    bp = ax.boxplot(
        [d_vals, v_vals],
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

    # Strip plot overlay (single RNG, advances across groups)
    rng = np.random.default_rng(42)
    for i, (vals, color) in enumerate([
        (d_vals, COLORS["Dorsal"]),
        (v_vals, COLORS["Ventral"]),
    ]):
        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(i + jitter, vals, color=color, alpha=0.6, s=20,
                   edgecolors="white", linewidth=0.3, zorder=3)

    # n computed from dropna'd series
    n_d, n_v = len(d_vals), len(v_vals)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Dorsal S4\nn={n_d}", f"Ventral S4\nn={n_v}"],
                       fontweight="bold")
    ax.set_ylabel("Lateral Spread (um)")
    ax.set_title("Terminal Arbor Width: Dorsal vs Ventral S4")

    # Significance annotation
    y_max = max(d_vals.max(), v_vals.max())
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


def plot_effect_heatmap(results: pd.DataFrame) -> None:
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

def main() -> None:
    """Run the full analysis: load data, compute stats, generate figures."""
    configure_plot_style()
    log_package_versions()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = load_data(DATA_PATH)
    validate_data(df)
    print_dataset_summary(df)

    results = run_comparisons(df)
    print_results(results)

    # Save statistical results to CSV for reproducibility
    results_path = OUTPUT_DIR / "stats_results.csv"
    results.to_csv(results_path, index=False)
    print(f"Saved: {results_path}")

    plot_ipl_depth_by_type(df)
    plot_lateral_spread_regional(df, results)
    plot_effect_heatmap(results)

    print("\nDone. Figures saved to figures/")


if __name__ == "__main__":
    main()
