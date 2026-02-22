"""
Sensitivity Analysis: Clean vs Raw Dataset
===========================================
Tests whether the key finding (ventral S4 cells have wider lateral spread
than dorsal S4 cells) holds when using the full raw dataset (pre-QC)
versus the QC-filtered clean dataset.

This matters because QC disproportionately excluded Ventral S4 cells
(23 of 25 S4 exclusions were Ventral, mostly IPL depth outliers).

Author: Vanz Labitad, University of Sussex
"""

import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths & Constants
# ---------------------------------------------------------------------------

RESULTS_DIR = Path(__file__).parent.parent / "morphometrics_results"
OUTPUT_DIR = Path(__file__).parent

COLORS = {"Dorsal": "#228833", "Ventral": "#CCBB44"}

# Significance threshold — uncorrected (two independent dataset comparisons).
ALPHA = 0.05

# Effect size retention threshold: finding considered "held" if raw d >= this
# fraction of clean d. Threshold is conservative but arbitrary; no established
# standard exists for this type of sensitivity comparison. Document in methods.
EFFECT_RETENTION_THRESHOLD = 0.8


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def load_and_filter(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load CSV and extract Dorsal S4 vs Ventral S4 lateral_spread_um arrays.
    Raises ValueError if either group is empty after filtering.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Data file not found: {path.resolve()}\n"
            "Check RESULTS_DIR path."
        )
    df = pd.read_csv(path)
    s4 = df[df["cell_type"] == "S4"]
    dorsal = s4[s4["region"] == "Dorsal"]["lateral_spread_um"].dropna().values
    ventral = s4[s4["region"] == "Ventral"]["lateral_spread_um"].dropna().values

    if len(dorsal) == 0 or len(ventral) == 0:
        raise ValueError(
            f"Empty group after filtering '{path.name}': "
            f"Dorsal n={len(dorsal)}, Ventral n={len(ventral)}. "
            f"Check 'cell_type' and 'region' column values."
        )

    return dorsal, ventral


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """
    Cohen's d with pooled standard deviation.
    Positive values indicate group a > group b (call as cohens_d(dorsal, ventral)).

    NOTE: Cohen's d is a parametric effect size and is technically mismatched with
    Mann-Whitney U (a non-parametric test). We report it here as a *descriptive*
    standardised mean difference because: (1) it is widely understood and directly
    comparable across metrics, (2) our sample sizes (~20-60 per group) are large
    enough that the mean/SD are meaningful summaries even when normality is violated,
    and (3) rank-biserial r is reported alongside for non-parametric completeness.
    This decision should be stated explicitly in the dissertation methods section.
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


def run_test(dorsal: np.ndarray, ventral: np.ndarray, label: str) -> dict:
    """Run Mann-Whitney U test and compute effect sizes. Returns result dict."""
    stat, p = stats.mannwhitneyu(dorsal, ventral, alternative="two-sided")
    d = cohens_d(dorsal, ventral)
    r = rank_biserial_r(stat, len(dorsal), len(ventral))
    return {
        "dataset": label,
        "n_dorsal": len(dorsal),
        "n_ventral": len(ventral),
        "median_dorsal": np.median(dorsal),
        "median_ventral": np.median(ventral),
        "U_stat": stat,
        "p_value": p,
        "cohens_d": d,
        "rank_biserial_r": r,
        "significant": p < ALPHA,
    }


def save_summary(results: list[dict], path: Path) -> None:
    """Write sensitivity_comparison.csv."""
    df = pd.DataFrame(results)
    df.to_csv(path, index=False)
    print(f"Saved: {path}")


def plot_comparison(clean_d: np.ndarray, clean_v: np.ndarray,
                    raw_d: np.ndarray, raw_v: np.ndarray,
                    clean_stats: dict, raw_stats: dict, path: Path) -> None:
    """Side-by-side boxplot: clean (left) vs raw (right), matched y-axis."""
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(9, 5))

    panels = [
        (ax1, clean_d, clean_v, clean_stats,
         f"Clean Dataset (n={clean_stats['n_dorsal'] + clean_stats['n_ventral']})"),
        (ax2, raw_d, raw_v, raw_stats,
         f"Raw Dataset (n={raw_stats['n_dorsal'] + raw_stats['n_ventral']})"),
    ]

    # Single RNG instance — advances continuously across both panels
    rng = np.random.default_rng(42)

    for ax, d_vals, v_vals, st, title in panels:
        # Boxplot
        bp = ax.boxplot(
            [d_vals, v_vals], positions=[0, 1], widths=0.4,
            patch_artist=True, showfliers=False, zorder=2,
        )
        for box, color in zip(bp["boxes"], [COLORS["Dorsal"], COLORS["Ventral"]]):
            box.set_facecolor(color)
            box.set_alpha(0.5)
        for element in ["whiskers", "caps", "medians"]:
            for line in bp[element]:
                line.set_color("black")
                line.set_linewidth(1.2)

        # Individual data points
        for i, (vals, color) in enumerate([
            (d_vals, COLORS["Dorsal"]),
            (v_vals, COLORS["Ventral"]),
        ]):
            jitter = rng.uniform(-0.12, 0.12, len(vals))
            ax.scatter(i + jitter, vals, color=color, alpha=0.6, s=18,
                       edgecolors="white", linewidth=0.3, zorder=3)

        ax.set_xticks([0, 1])
        ax.set_xticklabels([
            f"Dorsal S4\nn={st['n_dorsal']}",
            f"Ventral S4\nn={st['n_ventral']}",
        ], fontweight="bold", fontsize=10)
        ax.set_title(title, fontsize=12, fontweight="bold")

        # Annotation box with stats
        p_str = "p < 0.001" if st["p_value"] < 0.001 else f"p = {st['p_value']:.4f}"
        ax.text(
            0.97, 0.97,
            f"{p_str}\nCohen's d = {st['cohens_d']:.2f}",
            transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="gray", alpha=0.9),
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax1.set_ylabel("Lateral Spread (um)", fontsize=12)
    fig.suptitle(
        "Sensitivity Analysis: Dorsal vs Ventral S4 Lateral Spread",
        fontsize=13, fontweight="bold", y=1.01,
    )

    plt.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def interpret(clean_stats: dict, raw_stats: dict) -> None:
    """Print brief interpretation of sensitivity analysis."""
    print("\n=== INTERPRETATION ===")
    print(f"Clean dataset: d = {clean_stats['cohens_d']:.2f}, "
          f"p = {clean_stats['p_value']:.2e}, "
          f"n = {clean_stats['n_dorsal']}D + {clean_stats['n_ventral']}V")
    print(f"Raw dataset:   d = {raw_stats['cohens_d']:.2f}, "
          f"p = {raw_stats['p_value']:.2e}, "
          f"n = {raw_stats['n_dorsal']}D + {raw_stats['n_ventral']}V")

    # Compare effect sizes
    d_clean = abs(clean_stats["cohens_d"])
    d_raw = abs(raw_stats["cohens_d"])
    both_sig = clean_stats["significant"] and raw_stats["significant"]

    if both_sig and d_raw >= d_clean * EFFECT_RETENTION_THRESHOLD:
        verdict = "HELD"
        detail = ("The key finding is robust to QC exclusions. The effect "
                  "remains significant and large in both datasets.")
    elif both_sig and d_raw < d_clean * EFFECT_RETENTION_THRESHOLD:
        verdict = "WEAKENED"
        detail = ("The finding remains significant but the effect size "
                  "decreases substantially in the raw dataset, suggesting "
                  "QC exclusions may inflate the effect.")
    elif clean_stats["significant"] and not raw_stats["significant"]:
        verdict = "LOST SIGNIFICANCE"
        detail = ("The finding is not significant in the raw dataset, "
                  "indicating it may depend on QC exclusion choices.")
    elif not clean_stats["significant"] and raw_stats["significant"]:
        verdict = "RAW ONLY"
        detail = ("The finding is significant only in the raw dataset. "
                  "QC filtering may be removing genuine signal.")
    else:
        verdict = "INCONCLUSIVE"
        detail = "The finding is not significant in either dataset."

    print(f"\nVerdict: {verdict}")
    print(detail)


def log_package_versions() -> None:
    """Print Python and package versions for reproducibility."""
    import scipy
    print(f"Python: {sys.version}")
    print(f"numpy: {np.__version__}, pandas: {pd.__version__}, "
          f"scipy: {scipy.__version__}, matplotlib: {plt.matplotlib.__version__}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    log_package_versions()

    # Load both datasets
    clean_d, clean_v = load_and_filter(RESULTS_DIR / "morphometrics_clean.csv")
    raw_d, raw_v = load_and_filter(RESULTS_DIR / "morphometrics_raw.csv")

    print(f"Clean: {len(clean_d)} Dorsal S4, {len(clean_v)} Ventral S4")
    print(f"Raw:   {len(raw_d)} Dorsal S4, {len(raw_v)} Ventral S4")

    # Run tests
    clean_stats = run_test(clean_d, clean_v, "clean")
    raw_stats = run_test(raw_d, raw_v, "raw")

    # Save summary CSV
    save_summary([clean_stats, raw_stats],
                 OUTPUT_DIR / "sensitivity_comparison.csv")

    # Generate figure
    plot_comparison(clean_d, clean_v, raw_d, raw_v,
                    clean_stats, raw_stats,
                    OUTPUT_DIR / "sensitivity_boxplot.png")

    # Interpretation
    interpret(clean_stats, raw_stats)


if __name__ == "__main__":
    main()
