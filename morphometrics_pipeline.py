"""
Zebrafish Bipolar Cell 3D Morphometrics Pipeline
Author: Vanz Labitad - MSci Dissertation, University of Sussex
Dataset: webKnossos EM reconstruction - Retinal bipolar cells (S2 and S4)

This pipeline analyzes 3D mesh files (.stl) to compare morphological properties
of S2 (OFF-pathway) and S4 (ON-pathway) bipolar cells between Dorsal and Ventral regions.
"""

import struct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from scipy.spatial import ConvexHull
import os
import json
from pathlib import Path
import warnings
from datetime import datetime
import sys
import importlib.metadata

# color scheme
COLORS = {
    'S2': '#E74C3C',           # red — OFF pathway
    'S2_light': '#F1948A',     # light red
    'S2_dark': '#C0392B',      # dark red
    'S4': '#2980B9',           # blue — ON pathway
    'S4_light': '#85C1E2',     # light blue
    'S4_dark': '#1F618D',      # dark blue
    'Dorsal': '#27AE60',       # green
    'Ventral': '#F39C12',      # orange
    'VL': '#8E44AD',           # purple
    'KM': '#16A085',           # teal
    'gray': '#7F8C8D',         # neutral gray
    'background': '#F8F9F9',   # light background
    'grid': '#E5E7E9'          # grid color
}

# Terminal size thresholds (µm) for lateral_spread categorisation.
# Based on observed distribution in this dataset — no published standard exists.
# Treat as exploratory classification; document in methods if used.
TERMINAL_SIZE_SMALL_UM = 5.0
TERMINAL_SIZE_MEDIUM_UM = 15.0

# Significance threshold — uncorrected; see run_all_comparisons docstring for justification.
ALPHA = 0.05


def configure_plot_style() -> None:
    """Apply quality matplotlib and seaborn defaults. Call from main()."""
    plt.rcParams.update({
        'font.size': 11,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica'],
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 16,
        'axes.linewidth': 1.5,
        'grid.linewidth': 0.8,
        'lines.linewidth': 2.0,
        'patch.linewidth': 1.5,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'xtick.minor.width': 1.0,
        'ytick.minor.width': 1.0,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white'
    })
    sns.set_style("ticks", {
        'axes.grid': True,
        'grid.linestyle': ':',
        'grid.color': '#CCCCCC'
    })
    sns.set_context("paper", font_scale=1.2)


def log_package_versions(output_dir: Path) -> None:
    """Record Python and package versions alongside results for reproducibility."""
    packages = ['numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn']
    versions = {p: importlib.metadata.version(p) for p in packages}
    versions['python'] = sys.version
    with open(output_dir / 'package_versions.json', 'w') as f:
        json.dump(versions, f, indent=2)
    print("OK Package versions saved")

# ============================================================================
# STL BINARY PARSER
# ============================================================================

def read_binary_stl(filepath: str | Path) -> tuple[np.ndarray | None, list | None]:
    """
    Read binary STL file. Returns (normals, triangles) arrays or (None, None) if corrupt.
    Coordinates are in nanometres.
    """
    try:
        with open(filepath, 'rb') as f:
            f.read(80)  # skip header
            data = f.read(4)
            if len(data) < 4:
                return None, None
            num_tris = struct.unpack('<I', data)[0]
            if num_tris == 0 or num_tris > 10_000_000:
                return None, None
            normals, triangles = [], []
            for _ in range(num_tris):
                try:
                    n = struct.unpack('<3f', f.read(12))
                    v1 = struct.unpack('<3f', f.read(12))
                    v2 = struct.unpack('<3f', f.read(12))
                    v3 = struct.unpack('<3f', f.read(12))
                    f.read(2)  # attribute byte count
                    normals.append(n)
                    triangles.append((v1, v2, v3))
                except struct.error:
                    break
        if len(triangles) < 10:
            return None, None
        return np.array(normals), triangles
    except Exception as e:
        print(f"WARNING: Failed to read STL file {filepath}: {e}")
        return None, None

# ============================================================================
# FILENAME PARSER
# ============================================================================

def parse_filename(fname: str) -> tuple[str | None, str | None, str | None, str | None]:
    """
    Extract annotator, region, cell_type, cell_number from filename.
    Expected format: annotator_region_type_number (e.g. vl_dorsal_s2_01).
    Returns (None, None, None, None) and prints a warning if parsing fails.
    """
    stem = Path(fname).stem.lower()
    parts = stem.split('_')
    if len(parts) >= 4:
        annotator = parts[0].upper()   # VL or KM
        region = parts[1].capitalize()  # Dorsal or Ventral
        cell_type = parts[2].upper()    # S2 or S4
        cell_num = parts[3]
        return annotator, region, cell_type, cell_num
    print(f"WARNING: Could not parse filename '{fname}' (expected annotator_region_celltype_num)")
    return None, None, None, None

# ============================================================================
# MORPHOMETRIC COMPUTATION FUNCTIONS
# ============================================================================

def compute_volume_um3(triangles: list) -> float:
    """Calculate volume using signed tetrahedron method. Returns volume in µm³."""
    vol = 0.0
    for v1, v2, v3 in triangles:
        v1, v2, v3 = np.array(v1), np.array(v2), np.array(v3)
        vol += np.dot(v1, np.cross(v2, v3)) / 6.0
    return abs(vol) / 1e9  # nm³ → µm³

def compute_surface_area_um2(triangles: list) -> float:
    """Calculate surface area from triangle meshes. Returns area in µm²."""
    area = 0.0
    for v1, v2, v3 in triangles:
        v1, v2, v3 = np.array(v1), np.array(v2), np.array(v3)
        area += 0.5 * np.linalg.norm(np.cross(v2 - v1, v3 - v1))
    return area / 1e6  # nm² → µm²

def compute_convex_hull(triangles: list) -> tuple[float, float]:
    """
    Compute convex hull volume and surface area. Returns (volume_um3, area_um2).
    Vertices are converted from nm to µm before hull computation.
    If >50,000 vertices, a reproducible subsample is taken for performance.
    """
    verts_nm = np.array([v for tri in triangles for v in tri])
    # Subsample for speed if very large mesh (seeded for reproducibility)
    if len(verts_nm) > 50000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(verts_nm), 50000, replace=False)
        verts_nm = verts_nm[idx]
    verts_um = verts_nm / 1000  # nm → µm
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            hull = ConvexHull(verts_um)
        return hull.volume, hull.area  # µm³, µm²
    except Exception:
        return np.nan, np.nan

def compute_all_metrics(triangles: list) -> dict:
    """
    Compute all morphological metrics for a single cell.
    Returns dictionary of metrics.
    """
    # Basic metrics
    n_triangles = len(triangles)
    volume_um3 = compute_volume_um3(triangles)
    surface_area_um2 = compute_surface_area_um2(triangles)

    # Extract all vertices
    all_verts = np.array([v for tri in triangles for v in tri])  # shape: (N, 3)
    xs, ys, zs = all_verts[:,0] / 1000, all_verts[:,1] / 1000, all_verts[:,2] / 1000  # → µm

    # Bounding box metrics
    x_span = xs.max() - xs.min()   # µm
    y_span = ys.max() - ys.min()   # µm
    z_span = zs.max() - zs.min()   # µm
    x_centroid = xs.mean()         # µm
    y_centroid = ys.mean()         # µm
    z_centroid = zs.mean()         # µm
    lateral_spread = max(x_span, y_span)  # µm
    xy_footprint = x_span * y_span        # µm²

    # Convex hull metrics
    hull_volume, hull_area = compute_convex_hull(triangles)

    # Shape complexity metrics
    sa_v_ratio = surface_area_um2 / volume_um3 if volume_um3 > 0 else np.nan
    sphericity = (np.pi**(1/3) * (6*volume_um3)**(2/3)) / surface_area_um2 if surface_area_um2 > 0 else np.nan
    aspect_ratio = z_span / lateral_spread if lateral_spread > 0 else np.nan
    # Custom metric: surface complexity relative to XY footprint.
    # Higher values indicate more elaborated 3D structure per unit retinal area.
    # Not a published standard — treat as exploratory.
    branching_index = surface_area_um2 / xy_footprint if xy_footprint > 0 else np.nan

    # Terminal size category (see module-level constants for thresholds)
    if lateral_spread < TERMINAL_SIZE_SMALL_UM:
        terminal_size = 'Small'
    elif lateral_spread <= TERMINAL_SIZE_MEDIUM_UM:
        terminal_size = 'Medium'
    else:
        terminal_size = 'Large'

    return {
        'n_triangles': n_triangles,
        'volume_um3': volume_um3,
        'surface_area_um2': surface_area_um2,
        'x_span_um': x_span,
        'y_span_um': y_span,
        'z_span_um': z_span,
        'x_centroid_um': x_centroid,
        'y_centroid_um': y_centroid,
        'z_centroid_um': z_centroid,
        'lateral_spread_um': lateral_spread,
        'xy_footprint_um2': xy_footprint,
        'hull_volume_um3': hull_volume,
        'hull_area_um2': hull_area,
        'sa_v_ratio': sa_v_ratio,
        'sphericity': sphericity,
        'aspect_ratio': aspect_ratio,
        'branching_index': branching_index,
        'terminal_size': terminal_size
    }

# ============================================================================
# IPL STRATIFICATION NORMALIZATION
# ============================================================================

def normalize_ipl_depth(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Apply empirical IPL stratification normalization.
    Uses 2nd and 98th percentile of z_centroid to define IPL boundaries.
    Automatically inverts z-axis if S2 > S4 (incorrect orientation).
    Returns modified dataframe with ipl_depth_pct column and calibration info.
    """
    df = df.copy()
    all_z_centroids = df['z_centroid_um'].values
    ipl_z_min = np.percentile(all_z_centroids, 2)   # ~INL/IPL border
    ipl_z_max = np.percentile(all_z_centroids, 98)  # ~IPL/GCL border
    ipl_span = ipl_z_max - ipl_z_min

    df['ipl_depth_pct'] = ((df['z_centroid_um'] - ipl_z_min) / ipl_span * 100).clip(0, 100)

    print(f"\n=== IPL CALIBRATION ===")
    print(f"Empirical IPL boundaries: {ipl_z_min:.2f} to {ipl_z_max:.2f} um (span={ipl_span:.2f} um)")
    print(f"Expected: S2 cells ~17-34%, S4 cells ~50-67% (Connaughton 2004)")

    # Check if z-axis needs inversion
    s2_mean = df[df['cell_type']=='S2']['ipl_depth_pct'].mean()
    s4_mean = df[df['cell_type']=='S4']['ipl_depth_pct'].mean()

    inverted = False
    if s2_mean > s4_mean:
        df['ipl_depth_pct'] = 100 - df['ipl_depth_pct']
        print("WARNING  NOTE: Z-axis inverted - corrected automatically")
        inverted = True
        s2_mean = df[df['cell_type']=='S2']['ipl_depth_pct'].mean()
        s4_mean = df[df['cell_type']=='S4']['ipl_depth_pct'].mean()

    print(f"S2 mean depth: {s2_mean:.1f}%")
    print(f"S4 mean depth: {s4_mean:.1f}%")

    # Validation
    s2_valid = 17 <= s2_mean <= 34
    s4_valid = 50 <= s4_mean <= 67

    if s2_valid and s4_valid:
        status = "VALID"
    else:
        status = "WARNING  WARNING: Outside expected range"

    print(f"Calibration status: {status}")

    calibration_info = {
        'ipl_z_min_um': float(ipl_z_min),
        'ipl_z_max_um': float(ipl_z_max),
        'ipl_span_um': float(ipl_span),
        's2_mean_depth_pct': float(s2_mean),
        's4_mean_depth_pct': float(s4_mean),
        'z_axis_inverted': bool(inverted),
        'calibration_valid': bool(s2_valid and s4_valid),
        'timestamp': datetime.now().isoformat()
    }

    return df, calibration_info

# ============================================================================
# QUALITY CONTROL FILTERS
# ============================================================================

def apply_qc_filters(df: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    """
    Apply QC filters and return (clean_df, exclusions_list).
    """
    exclusions = []

    # 1. Too few triangles
    mask_tris = df['n_triangles'] < 100
    for idx, row in df[mask_tris].iterrows():
        exclusions.append({
            'filename': row['filename'],
            'cell_type': row['cell_type'],
            'region': row['region'],
            'annotator': row['annotator'],
            'reason': 'Too few triangles (<100)'
        })

    # 2. Volume outliers (>3 SD from group mean)
    for (cell_type, region), group in df.groupby(['cell_type', 'region']):
        mean_v = group['volume_um3'].mean()
        sd_v = group['volume_um3'].std()
        outliers = group[group['volume_um3'] > mean_v + 3*sd_v]
        for idx, row in outliers.iterrows():
            exclusions.append({
                'filename': row['filename'],
                'cell_type': row['cell_type'],
                'region': row['region'],
                'annotator': row['annotator'],
                'reason': f'Volume outlier >3SD (group mean={mean_v:.1f}, SD={sd_v:.1f})'
            })

    # 3. Lateral spread too large
    mask_spread = df['lateral_spread_um'] > 30
    for idx, row in df[mask_spread].iterrows():
        exclusions.append({
            'filename': row['filename'],
            'cell_type': row['cell_type'],
            'region': row['region'],
            'annotator': row['annotator'],
            'reason': 'Lateral spread >30µm — likely merge artefact'
        })

    # 4. IPL depth implausible
    mask_depth = (df['ipl_depth_pct'] < 5) | (df['ipl_depth_pct'] > 95)
    for idx, row in df[mask_depth].iterrows():
        exclusions.append({
            'filename': row['filename'],
            'cell_type': row['cell_type'],
            'region': row['region'],
            'annotator': row['annotator'],
            'reason': f'IPL depth {row["ipl_depth_pct"]:.1f}% outside expected range (5-95%)'
        })

    # Remove excluded cells (deduplicate — a cell can fail multiple filters)
    excluded_filenames = list({e['filename'] for e in exclusions})
    df_clean = df[~df['filename'].isin(excluded_filenames)].copy()

    print(f"\n=== QC FILTERING ===")
    print(f"Total cells: {len(df)}")
    print(f"Excluded: {len(excluded_filenames)} cells ({len(exclusions)} total filter hits)")
    print(f"Remaining: {len(df_clean)}")

    return df_clean, exclusions

# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

def run_comparison(group_a: pd.DataFrame, group_b: pd.DataFrame,
                   label_a: str, label_b: str, metric: str) -> dict | None:
    """
    Run statistical comparison between two groups for a given metric.
    Returns dict with test results, or None if insufficient data.
    """
    a = group_a[metric].dropna().values
    b = group_b[metric].dropna().values

    if len(a) < 3 or len(b) < 3:
        return None

    # Shapiro-Wilk
    norm_a = stats.shapiro(a).pvalue > 0.05 if len(a) >= 3 else False
    norm_b = stats.shapiro(b).pvalue > 0.05 if len(b) >= 3 else False

    # Choose test
    if norm_a and norm_b and len(a) >= 5 and len(b) >= 5:
        stat, p = stats.ttest_ind(a, b)
        test_used = 't-test'
    else:
        stat, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        test_used = 'Mann-Whitney U'

    # Effect size: Cohen's d (pooled SD)
    # NOTE: Cohen's d is a parametric effect size and is technically mismatched with
    # Mann-Whitney U (a non-parametric test). We report it here as a *descriptive*
    # standardised mean difference because: (1) it is widely understood and directly
    # comparable across metrics, (2) our sample sizes (~20-60 per group) are large
    # enough that the mean/SD are meaningful summaries even when normality is violated,
    # and (3) rank-biserial r is reported alongside for non-parametric completeness.
    # This decision should be stated explicitly in the dissertation methods section.
    n_a, n_b = len(a), len(b)
    pooled_sd = np.sqrt(((n_a-1)*np.std(a,ddof=1)**2 + (n_b-1)*np.std(b,ddof=1)**2) / (n_a+n_b-2)) if (n_a + n_b) > 2 else 0.0
    cohens_d = (np.mean(a) - np.mean(b)) / pooled_sd if pooled_sd > 0 else np.nan

    # Rank-biserial correlation: non-parametric effect size for Mann-Whitney U
    # r = 1 - (2U)/(n_a * n_b), bounded [-1, 1]
    rank_biserial_r = 1 - (2 * stat) / (n_a * n_b) if (n_a * n_b) > 0 else np.nan

    # Effect magnitude (based on Cohen's d thresholds)
    if abs(cohens_d) < 0.5:
        effect_mag = 'small'
    elif abs(cohens_d) < 0.8:
        effect_mag = 'medium'
    else:
        effect_mag = 'large'

    return {
        'comparison': f'{label_a} vs {label_b}',
        'metric': metric,
        'n_a': n_a, 'n_b': n_b,
        'mean_a': np.mean(a), 'sd_a': np.std(a, ddof=1),
        'mean_b': np.mean(b), 'sd_b': np.std(b, ddof=1),
        'median_a': np.median(a), 'median_b': np.median(b),
        'test': test_used,
        'statistic': stat, 'p_value': p,
        'cohens_d': cohens_d,
        'rank_biserial_r': rank_biserial_r,
        'significant': p < ALPHA,
        'effect_magnitude': effect_mag
    }

def run_all_comparisons(df: pd.DataFrame) -> pd.DataFrame:
    """
    Run all statistical comparisons.

    NOTE on multiple comparisons: 8 metrics × 4 comparison families = 32 tests.
    No FDR/Bonferroni correction is applied here because these are pre-specified,
    hypothesis-driven comparisons (not exploratory screening). The 'significant'
    column uses uncorrected p < 0.05. Note to self - this decision must be stated explicitly in
    the dissertation methods section.
    """
    metrics = ['volume_um3', 'surface_area_um2', 'lateral_spread_um', 'z_span_um',
               'ipl_depth_pct', 'sa_v_ratio', 'sphericity', 'aspect_ratio']

    results = []

    print("\n=== STATISTICAL ANALYSIS ===")

    # 1. S2 vs S4 (pooled)
    print("Running: S2 vs S4 (pooled)")
    for metric in metrics:
        res = run_comparison(
            df[df['cell_type']=='S2'],
            df[df['cell_type']=='S4'],
            'S2', 'S4', metric
        )
        if res:
            results.append(res)

    # 2. Dorsal vs Ventral (pooled)
    print("Running: Dorsal vs Ventral (pooled)")
    for metric in metrics:
        res = run_comparison(
            df[df['region']=='Dorsal'],
            df[df['region']=='Ventral'],
            'Dorsal', 'Ventral', metric
        )
        if res:
            results.append(res)

    # 3. Dorsal S4 vs Ventral S4
    print("Running: Dorsal S4 vs Ventral S4")
    for metric in metrics:
        res = run_comparison(
            df[(df['region']=='Dorsal') & (df['cell_type']=='S4')],
            df[(df['region']=='Ventral') & (df['cell_type']=='S4')],
            'Dorsal S4', 'Ventral S4', metric
        )
        if res:
            results.append(res)

    # 4. VL vs KM (inter-rater)
    print("Running: VL vs KM (inter-rater reliability)")
    for metric in metrics:
        res = run_comparison(
            df[df['annotator']=='VL'],
            df[df['annotator']=='KM'],
            'VL', 'KM', metric
        )
        if res:
            results.append(res)

    return pd.DataFrame(results)

# ============================================================================
# FIGURE GENERATION
# ============================================================================

def create_figure1_stratification(df: pd.DataFrame, output_dir: Path) -> None:
    """Figure 1: IPL Stratification Overview - Enhanced with violin + swarm plots"""
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for i, region in enumerate(['Dorsal', 'Ventral']):
        ax = axes[i]
        data = df[df['region'] == region]

        # Enhanced reference zones with borders
        s2_patch = ax.axvspan(17, 34, alpha=0.15, color=COLORS['S2_light'], zorder=0)
        s4_patch = ax.axvspan(50, 67, alpha=0.15, color=COLORS['S4_light'], zorder=0)

        # Add borders to reference zones
        ax.axvline(17, color=COLORS['S2'], linestyle='--', linewidth=1.5, alpha=0.6)
        ax.axvline(34, color=COLORS['S2'], linestyle='--', linewidth=1.5, alpha=0.6)
        ax.axvline(50, color=COLORS['S4'], linestyle='--', linewidth=1.5, alpha=0.6)
        ax.axvline(67, color=COLORS['S4'], linestyle='--', linewidth=1.5, alpha=0.6)

        # Add zone labels
        ax.text(25.5, 0.42, 'S2 Expected', ha='center', fontsize=9,
               color=COLORS['S2_dark'], fontweight='bold')
        ax.text(58.5, 0.42, 'S4 Expected', ha='center', fontsize=9,
               color=COLORS['S4_dark'], fontweight='bold')

        # Violin plots for density visualization
        for j, cell_type in enumerate(['S2', 'S4']):
            subset = data[data['cell_type'] == cell_type]
            if len(subset) > 0:
                # Create violin plot data
                y_pos = 0.15 if cell_type == 'S2' else -0.15
                parts = ax.violinplot([subset['ipl_depth_pct'].values],
                                     positions=[y_pos],
                                     vert=False, widths=0.25,
                                     showmeans=False, showmedians=True)

                # Color the violin
                for pc in parts['bodies']:
                    pc.set_facecolor(COLORS[cell_type])
                    pc.set_alpha(0.4)
                    pc.set_edgecolor(COLORS[cell_type])
                    pc.set_linewidth(1.5)

                # Style the median line
                parts['cmedians'].set_color(COLORS[f'{cell_type}_dark'])
                parts['cmedians'].set_linewidth(2.5)

                # Swarm plot overlay for individual cells
                y_jitter = np.random.normal(y_pos, 0.03, len(subset))
                ax.scatter(subset['ipl_depth_pct'], y_jitter,
                          c=COLORS[cell_type], s=35, alpha=0.7,
                          edgecolors=COLORS[f'{cell_type}_dark'], linewidths=0.8,
                          label=f'{cell_type} (n={len(subset)})', zorder=3)

        # Styling
        ax.set_ylabel(region, fontsize=14, fontweight='bold')
        ax.set_ylim(-0.5, 0.5)
        ax.set_yticks([])
        ax.grid(axis='x', alpha=0.4, linestyle=':', linewidth=1)
        ax.set_axisbelow(True)

        # Legend
        if i == 0:
            legend = ax.legend(loc='upper left', frameon=True, framealpha=0.95,
                             edgecolor='gray', fancybox=True, shadow=True)
            legend.get_frame().set_linewidth(1.5)

    # X-axis label and limits
    axes[1].set_xlabel('IPL Depth (%)', fontsize=13, fontweight='bold')
    axes[1].set_xlim(0, 100)

    # Add minor grid lines for better readability
    for ax in axes:
        ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
        ax.grid(which='minor', axis='x', alpha=0.2, linestyle=':', linewidth=0.5)

    plt.suptitle('IPL Stratification Depth of S2 and S4 Bipolar Cells',
                 fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()

    plt.savefig(output_dir / 'figure1_stratification.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure1_stratification.svg', bbox_inches='tight')
    plt.close()
    print("OK Figure 1 saved")

def create_figure2_regional_comparison(df: pd.DataFrame, stats_df: pd.DataFrame, output_dir: Path) -> None:
    """Figure 2: Regional Morphology Comparison - Enhanced main dissertation figure"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    metrics = [
        ('volume_um3', 'Volume (um3)', 'Volume\n(um3)'),
        ('lateral_spread_um', 'Lateral Spread (um)', 'Lateral Spread\n(um)'),
        ('surface_area_um2', 'Surface Area (um2)', 'Surface Area\n(um2)')
    ]

    # Panel letters
    panel_letters = ['A', 'B', 'C', 'D', 'E', 'F']
    letter_idx = 0

    for row, region in enumerate(['Dorsal', 'Ventral']):
        data = df[df['region'] == region]

        # Add subtle background shading for row distinction
        if row == 1:
            for col in range(3):
                axes[row, col].set_facecolor(COLORS['background'])

        for col, (metric, ylabel, title) in enumerate(metrics):
            ax = axes[row, col]

            # Prepare data for plotting
            plot_data = []
            sample_sizes = {}
            for cell_type in ['S2', 'S4']:
                subset = data[data['cell_type'] == cell_type]
                sample_sizes[cell_type] = len(subset)
                for val in subset[metric].dropna():
                    plot_data.append({'Cell Type': cell_type, ylabel: val})

            plot_df = pd.DataFrame(plot_data)

            # Enhanced box plot with thicker lines
            box_props = dict(linewidth=2.5)
            whisker_props = dict(linewidth=2)
            median_props = dict(linewidth=3, color='black')
            flier_props = dict(marker='o', markerfacecolor='gray', markersize=5,
                             linestyle='none', markeredgecolor='gray', alpha=0.5)

            bp = ax.boxplot([plot_df[plot_df['Cell Type']=='S2'][ylabel].values,
                            plot_df[plot_df['Cell Type']=='S4'][ylabel].values],
                           positions=[0, 1],
                           widths=0.5,
                           patch_artist=True,
                           boxprops=box_props,
                           whiskerprops=whisker_props,
                           medianprops=median_props,
                           flierprops=flier_props)

            # Color the boxes
            for patch, cell_type in zip(bp['boxes'], ['S2', 'S4']):
                patch.set_facecolor(COLORS[cell_type])
                patch.set_alpha(0.7)
                patch.set_edgecolor(COLORS[f'{cell_type}_dark'])

            # Swarm plot overlay for individual points
            for i, cell_type in enumerate(['S2', 'S4']):
                subset_vals = plot_df[plot_df['Cell Type']==cell_type][ylabel].values
                x_pos = np.random.normal(i, 0.04, len(subset_vals))
                ax.scatter(x_pos, subset_vals, c=COLORS[f'{cell_type}_dark'],
                          s=25, alpha=0.6, edgecolors='black', linewidths=0.5, zorder=3)

            # Sample sizes embedded in x-tick labels (avoids clipping)
            n_s2 = sample_sizes['S2']
            n_s4 = sample_sizes['S4']

            # Add significance annotation
            stat_row = stats_df[(stats_df['comparison'].str.contains('S2 vs S4')) &
                               (stats_df['metric'] == metric)]
            if not stat_row.empty and stat_row.iloc[0]['significant']:
                p = stat_row.iloc[0]['p_value']
                if p < 0.001:
                    sig_text = '***'
                elif p < 0.01:
                    sig_text = '**'
                else:
                    sig_text = '*'

                y_max = plot_df[ylabel].max()
                y_range = plot_df[ylabel].max() - plot_df[ylabel].min()
                y_bar = y_max + 0.1*y_range

                # Significance bar
                ax.plot([0, 1], [y_bar, y_bar], 'k-', linewidth=2)
                ax.plot([0, 0], [y_bar-0.02*y_range, y_bar], 'k-', linewidth=2)
                ax.plot([1, 1], [y_bar-0.02*y_range, y_bar], 'k-', linewidth=2)

                # Significance text
                ax.text(0.5, y_bar + 0.03*y_range, sig_text,
                       ha='center', fontsize=14, fontweight='bold')

            # Panel letter
            ax.text(-0.15, 1.05, panel_letters[letter_idx], transform=ax.transAxes,
                   fontsize=18, fontweight='bold', va='top')
            letter_idx += 1

            # Labels and styling
            if col == 0:
                ax.set_ylabel(ylabel.split('(')[0].strip(), fontsize=12, fontweight='bold')
            else:
                ax.set_ylabel('')

            if row == 0:
                ax.set_title(title, fontsize=13, fontweight='bold', pad=10)

            ax.set_xticks([0, 1])
            ax.set_xticklabels([f'S2\n(n={n_s2})', f'S4\n(n={n_s4})'],
                              fontsize=10, fontweight='bold')
            ax.set_xlabel('')

            # Grid and styling
            ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=1)
            ax.set_axisbelow(True)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['bottom'].set_linewidth(1.5)

            # Add region label on left with background box
            if col == 0:
                bbox_props = dict(boxstyle='round,pad=0.5', facecolor=COLORS[region],
                                alpha=0.8, edgecolor='black', linewidth=2)
                ax.text(-0.25, 0.5, region, transform=ax.transAxes,
                       fontsize=14, fontweight='bold', rotation=90,
                       verticalalignment='center', ha='center',
                       bbox=bbox_props, color='white')

    plt.suptitle('Regional Morphology Comparison: Dorsal vs Ventral',
                fontsize=18, fontweight='bold', y=1.02)

    plt.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.08,
                       hspace=0.35, wspace=0.25)

    plt.savefig(output_dir / 'figure2_regional_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure2_regional_comparison.svg', bbox_inches='tight')
    plt.close()
    print("OK Figure 2 saved")

def create_figure3_interrater(df: pd.DataFrame, output_dir: Path) -> None:
    """Figure 3: Inter-Rater Reliability - Enhanced with color-coding and statistics"""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.5))

    metrics = [
        ('volume_um3', 'Volume (um3)', 'Volume\n(um3)'),
        ('lateral_spread_um', 'Lateral Spread (um)', 'Lateral Spread\n(um)'),
        ('ipl_depth_pct', 'IPL Depth (%)', 'IPL Depth\n(%)')
    ]

    for i, (metric, label, title) in enumerate(metrics):
        ax = axes[i]

        # Get VL and KM data
        vl_data = df[df['annotator']=='VL'].groupby(['region', 'cell_type', 'cell_number'])[metric].mean()
        km_data = df[df['annotator']=='KM'].groupby(['region', 'cell_type', 'cell_number'])[metric].mean()

        # Find common cells
        common_keys = set(vl_data.index) & set(km_data.index)

        if len(common_keys) > 0:
            vl_vals = np.array([vl_data[k] for k in common_keys])
            km_vals = np.array([km_data[k] for k in common_keys])

            # Color code by cell type
            colors_list = []
            for k in common_keys:
                cell_type = k[1]  # Extract cell_type from tuple
                colors_list.append(COLORS[cell_type])

            # Scatter plot with color coding
            scatter = ax.scatter(vl_vals, km_vals, c=colors_list, s=70, alpha=0.7,
                               edgecolors='black', linewidths=1, zorder=3)

            # Add identity line (thicker and more visible)
            min_val = min(min(vl_vals), min(km_vals))
            max_val = max(max(vl_vals), max(km_vals))
            padding = (max_val - min_val) * 0.05
            ax.plot([min_val-padding, max_val+padding],
                   [min_val-padding, max_val+padding],
                   'k--', linewidth=2.5, alpha=0.7, label='Identity', zorder=1)

            # Calculate statistics
            r, p = stats.pearsonr(vl_vals, km_vals)
            r_squared = r**2

            # Statistics box
            stats_text = f'r = {r:.3f}\nR$^2$ = {r_squared:.3f}\np = {p:.3f}\nn = {len(vl_vals)}'
            ax.text(0.05, 0.95, stats_text,
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.6', facecolor='white',
                           edgecolor='black', linewidth=1.5, alpha=0.95),
                   fontsize=10, fontweight='bold')

            # Note: aspect ratio left free to prevent compression with mismatched ranges

        # Labels and styling
        ax.set_xlabel(f'VL Annotator - {label}', fontsize=11, fontweight='bold')
        ax.set_ylabel(f'KM Annotator - {label}', fontsize=11, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
        ax.grid(alpha=0.4, linestyle=':', linewidth=1, zorder=0)
        ax.set_axisbelow(True)

        # Spines styling
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

    # Add legend for cell types
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='S2',
                  markerfacecolor=COLORS['S2'], markersize=10,
                  markeredgecolor='black', markeredgewidth=1),
        plt.Line2D([0], [0], marker='o', color='w', label='S4',
                  markerfacecolor=COLORS['S4'], markersize=10,
                  markeredgecolor='black', markeredgewidth=1)
    ]
    axes[0].legend(handles=legend_elements, loc='lower right', frameon=True,
                  framealpha=0.95, edgecolor='black', fancybox=True, shadow=True)

    plt.suptitle('Inter-Rater Reliability: VL vs KM Annotators',
                fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    plt.savefig(output_dir / 'figure3_interrater.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure3_interrater.svg', bbox_inches='tight')
    plt.close()
    print("OK Figure 3 saved")

def create_figure4_shape_complexity(df: pd.DataFrame, output_dir: Path) -> None:
    """Figure 4: Shape Complexity - Enhanced with confidence intervals"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for i, region in enumerate(['Dorsal', 'Ventral']):
        ax = axes[i]
        data = df[df['region'] == region]

        for cell_type in ['S2', 'S4']:
            subset = data[data['cell_type'] == cell_type]
            x_data = subset['volume_um3'].dropna().values
            y_data = subset['surface_area_um2'].dropna().values

            # Scatter plot with enhanced styling
            ax.scatter(x_data, y_data,
                      c=COLORS[cell_type], s=80, alpha=0.65,
                      edgecolors=COLORS[f'{cell_type}_dark'], linewidths=1.2,
                      label=f'{cell_type} (n={len(x_data)})', zorder=3)

            # Calculate regression with confidence interval
            if len(x_data) > 2:
                # Linear regression
                slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)

                # Generate line
                x_line = np.linspace(x_data.min(), x_data.max(), 100)
                y_line = slope * x_line + intercept

                # Plot regression line
                ax.plot(x_line, y_line, color=COLORS[cell_type],
                       linestyle='-', linewidth=2.5, alpha=0.9, zorder=2)

                # Calculate 95% confidence interval
                predict_error = np.sqrt(np.sum((y_data - (slope * x_data + intercept))**2) / (len(x_data) - 2))
                margin = 1.96 * predict_error  # 95% CI
                ax.fill_between(x_line, y_line - margin, y_line + margin,
                               color=COLORS[cell_type], alpha=0.15, zorder=1)

                # Add R² value annotation (fixed position, below legend)
                r_squared = r_value**2
                label_y = 0.78 if cell_type == 'S2' else 0.68
                ax.text(0.05, label_y, f'{cell_type}: R$^2$ = {r_squared:.3f}',
                       fontsize=10, fontweight='bold',
                       color=COLORS[f'{cell_type}_dark'],
                       transform=ax.transAxes,
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                               edgecolor=COLORS[cell_type], linewidth=1.5, alpha=0.9))

        # Labels and styling
        ax.set_xlabel('Volume (um3)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Surface Area (um2)', fontsize=12, fontweight='bold')
        ax.set_title(f'{region} Region', fontsize=14, fontweight='bold', pad=10)

        # Legend with enhanced styling
        legend = ax.legend(loc='upper left', frameon=True, framealpha=0.95,
                         edgecolor='black', fancybox=True, shadow=True,
                         fontsize=10)
        legend.get_frame().set_linewidth(1.5)

        # Grid and spines
        ax.grid(alpha=0.4, linestyle=':', linewidth=1, zorder=0)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

    plt.suptitle('Shape Complexity: Volume vs Surface Area Relationship',
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()

    plt.savefig(output_dir / 'figure4_shape_complexity.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure4_shape_complexity.svg', bbox_inches='tight')
    plt.close()
    print("OK Figure 4 saved")

def create_figure5_effect_heatmap(stats_df: pd.DataFrame, output_dir: Path) -> None:
    """Figure 5: Effect Size Heatmap - Enhanced with readable labels and significance markers"""

    # Create readable metric labels
    metric_labels = {
        'volume_um3': 'Volume',
        'surface_area_um2': 'Surface Area',
        'lateral_spread_um': 'Lateral Spread',
        'z_span_um': 'Z-Span',
        'ipl_depth_pct': 'IPL Depth',
        'sa_v_ratio': 'Surface/Volume Ratio',
        'sphericity': 'Sphericity',
        'aspect_ratio': 'Aspect Ratio'
    }

    # Shorten comparison labels
    comparison_labels = {
        'S2 vs S4': 'S2 vs S4',
        'Dorsal vs Ventral': 'Dorsal vs\nVentral',
        'Dorsal S4 vs Ventral S4': 'Dorsal S4 vs\nVentral S4',
        'VL vs KM': 'VL vs KM'
    }

    # Pivot data
    pivot = stats_df.pivot(index='metric', columns='comparison', values='cohens_d')
    pivot_p = stats_df.pivot(index='metric', columns='comparison', values='p_value')

    # Rename axes
    pivot = pivot.rename(index=metric_labels, columns=comparison_labels)
    pivot_p = pivot_p.rename(index=metric_labels, columns=comparison_labels)

    fig, ax = plt.subplots(figsize=(11, 9))

    # Create heatmap with custom formatting
    sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
                cbar_kws={'label': "Cohen's d (Effect Size)", 'shrink': 0.8},
                ax=ax, vmin=-5, vmax=5, linewidths=1.5, linecolor='white',
                annot_kws={'size': 11, 'weight': 'bold'},
                cbar=True)

    # Add asterisks for significance
    for i, metric in enumerate(pivot.index):
        for j, comparison in enumerate(pivot.columns):
            p_val = pivot_p.loc[metric, comparison]
            if p_val < 0.05:
                d_val = pivot.loc[metric, comparison]
                # Add asterisks
                if p_val < 0.001:
                    sig_marker = '***'
                elif p_val < 0.01:
                    sig_marker = '**'
                else:
                    sig_marker = '*'

                # Position asterisks above the number
                ax.text(j + 0.5, i + 0.3, sig_marker,
                       ha='center', va='center', fontsize=14,
                       color='black', fontweight='bold')

    # Styling
    ax.set_xlabel('Comparison', fontsize=13, fontweight='bold', labelpad=10)
    ax.set_ylabel('Morphometric Feature', fontsize=13, fontweight='bold', labelpad=10)
    ax.set_title("Effect Size Summary (Cohen's d)\n* p<0.05, ** p<0.01, *** p<0.001",
                fontsize=15, fontweight='bold', pad=15)

    # Rotate labels for better readability
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center', fontsize=11)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11)

    # Add effect size interpretation guide
    guide_text = ('Effect Size Guide:\n'
                 'Small: |d| < 0.5\n'
                 'Medium: 0.5 <= |d| < 0.8\n'
                 'Large: |d| >= 0.8')
    ax.text(1.35, 0.5, guide_text, transform=ax.transAxes,
           fontsize=9, va='center',
           bbox=dict(boxstyle='round,pad=0.8', facecolor='white',
                    edgecolor='black', linewidth=1.5))

    plt.tight_layout()

    plt.savefig(output_dir / 'figure5_effect_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'figure5_effect_heatmap.svg', bbox_inches='tight')
    plt.close()
    print("OK Figure 5 saved")

# ============================================================================
# REPORT GENERATION
# ============================================================================

def generate_report(df_raw: pd.DataFrame, df_clean: pd.DataFrame, exclusions: list[dict],
                    stats_df: pd.DataFrame, calibration_info: dict, output_dir: Path) -> None:
    """Generate comprehensive analysis report."""

    report = []
    report.append("=" * 70)
    report.append("ZEBRAFISH BIPOLAR CELL MORPHOMETRICS REPORT")
    report.append("=" * 70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Dataset: ZF_retina_opl_ipl_model_v5_segmentation_v1")
    report.append("")

    # Dataset summary
    report.append("--- DATASET SUMMARY ---")
    report.append(f"Total STL files found: {len(df_raw)}")
    report.append(f"Successfully loaded: {len(df_raw[df_raw['n_triangles'] > 0])}")
    report.append(f"Excluded (QC): {len(exclusions)}")
    report.append(f"Final dataset: {len(df_clean)}")
    report.append("")

    report.append("Breakdown:")
    for annotator in ['VL', 'KM']:
        for region in ['Dorsal', 'Ventral']:
            for cell_type in ['S2', 'S4']:
                n = len(df_clean[(df_clean['annotator']==annotator) &
                                (df_clean['region']==region) &
                                (df_clean['cell_type']==cell_type)])
                report.append(f"  {annotator} {region} {cell_type}: {n}")
    report.append("")

    # IPL calibration
    report.append("--- IPL CALIBRATION ---")
    report.append(f"Empirical z_min (0%): {calibration_info['ipl_z_min_um']:.2f} µm")
    report.append(f"Empirical z_max (100%): {calibration_info['ipl_z_max_um']:.2f} µm")
    report.append(f"IPL span: {calibration_info['ipl_span_um']:.2f} µm")
    report.append(f"S2 mean depth: {calibration_info['s2_mean_depth_pct']:.1f}% (expected: 17–34%)")
    report.append(f"S4 mean depth: {calibration_info['s4_mean_depth_pct']:.1f}% (expected: 50–67%)")
    if calibration_info['z_axis_inverted']:
        report.append("NOTE: Z-axis was inverted during calibration")
    report.append(f"Calibration status: {'VALID' if calibration_info['calibration_valid'] else 'WARNING'}")
    report.append("")

    # Key findings
    report.append("--- KEY FINDINGS ---")
    significant_results = stats_df[stats_df['significant']]
    if len(significant_results) > 0:
        for _, row in significant_results.iterrows():
            direction = "greater" if row['mean_a'] > row['mean_b'] else "lower"
            report.append(
                f"{row['comparison']} showed significant difference in {row['metric']}: "
                f"{row['test']}, p={row['p_value']:.4f}, d={row['cohens_d']:.2f} ({row['effect_magnitude']} effect)"
            )
    else:
        report.append("No significant differences detected.")
    report.append("")

    # Morphometric summary
    report.append("--- MORPHOMETRIC SUMMARY (mean ± SD) ---")
    metrics_to_report = ['volume_um3', 'surface_area_um2', 'lateral_spread_um',
                        'z_span_um', 'ipl_depth_pct']

    for metric in metrics_to_report:
        report.append(f"\n{metric}:")
        for region in ['Dorsal', 'Ventral']:
            for cell_type in ['S2', 'S4']:
                subset = df_clean[(df_clean['region']==region) & (df_clean['cell_type']==cell_type)]
                if len(subset) > 0:
                    mean_val = subset[metric].mean()
                    sd_val = subset[metric].std()
                    report.append(f"  {region} {cell_type}: {mean_val:.2f} ± {sd_val:.2f}")
    report.append("")

    # Inter-rater reliability
    report.append("--- INTER-RATER RELIABILITY ---")
    vl_km_stats = stats_df[stats_df['comparison'].str.contains('VL vs KM')]
    for metric in ['volume_um3', 'lateral_spread_um']:
        row = vl_km_stats[vl_km_stats['metric']==metric]
        if not row.empty:
            # Report stats from comparison test
            mean_diff = abs(row.iloc[0]['mean_a'] - row.iloc[0]['mean_b'])
            report.append(f"VL vs KM for {metric}: p={row.iloc[0]['p_value']:.3f}, mean_diff={mean_diff:.2f}")
    report.append("Assessment: See statistical_results.csv for full inter-rater analysis")
    report.append("")

    # QC flags
    report.append("--- QC FLAGS ---")
    if len(exclusions) > 0:
        for exc in exclusions:
            report.append(f"  {exc['filename']}: {exc['reason']}")
    else:
        report.append("  No exclusions")
    report.append("")

    # Limitations
    report.append("--- LIMITATIONS ---")
    report.append(f"- Small sample sizes per group (n~{len(df_clean)//8}) limit statistical power")
    report.append("- Empirical IPL normalisation assumes terminal arbors span full IPL range")
    report.append("- Bounding box lateral spread is an overestimate of true arbor spread")
    report.append("- Volume measurements include all cell compartments visible in mesh, not terminal-only")
    report.append("")

    report.append("=" * 70)

    # Write to file
    with open(output_dir / 'ANALYSIS_REPORT.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report))

    print("OK Analysis report saved")

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main() -> None:
    """Main execution pipeline."""
    configure_plot_style()

    print("=" * 70)
    print("ZEBRAFISH BIPOLAR CELL 3D MORPHOMETRICS PIPELINE")
    print("=" * 70)
    print()

    # Get working directory (accepts CLI argument or defaults to script directory)
    work_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent.resolve()

    # Create output directory
    output_dir = work_dir / 'morphometrics_results'
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Log package versions for reproducibility
    log_package_versions(output_dir)
    print()

    # Find all STL files
    stl_files = sorted(work_dir.glob('*.stl'))
    print(f"Found {len(stl_files)} STL files")
    print()

    # Load all files and compute metrics
    print("=== LOADING FILES AND COMPUTING METRICS ===")
    data = []
    failed_loads = []

    for i, filepath in enumerate(stl_files, 1):
        filename = filepath.name
        print(f"[{i}/{len(stl_files)}] Processing: {filename}", end=" ... ")

        # Parse filename
        annotator, region, cell_type, cell_num = parse_filename(filename)

        if annotator is None:
            print("FAIL Invalid filename format")
            failed_loads.append({'filename': filename, 'reason': 'Invalid filename format'})
            continue

        # Read STL
        normals, triangles = read_binary_stl(filepath)

        if triangles is None:
            print("FAIL Corrupt or empty file")
            failed_loads.append({'filename': filename, 'reason': 'Corrupt or empty STL file'})
            continue

        # Compute all metrics
        metrics = compute_all_metrics(triangles)

        # Combine with metadata
        row = {
            'filename': filename,
            'annotator': annotator,
            'region': region,
            'cell_type': cell_type,
            'cell_number': cell_num,
            **metrics
        }

        data.append(row)
        print(f"OK ({metrics['n_triangles']} triangles, {metrics['volume_um3']:.2f} um3)")

    print()
    print(f"Successfully loaded: {len(data)}/{len(stl_files)}")
    print(f"Failed to load: {len(failed_loads)}")
    print()

    # Create DataFrame
    df = pd.DataFrame(data)

    # Save raw data
    df.to_csv(output_dir / 'morphometrics_raw.csv', index=False)
    print("OK Raw data saved")

    # Apply IPL normalization
    df, calibration_info = normalize_ipl_depth(df)

    # Save calibration info
    with open(output_dir / 'ipl_calibration.json', 'w') as f:
        json.dump(calibration_info, f, indent=2)
    print("OK IPL calibration saved")

    # Apply QC filters
    df_clean, exclusions = apply_qc_filters(df)

    # Save clean data and exclusions
    df_clean.to_csv(output_dir / 'morphometrics_clean.csv', index=False)
    pd.DataFrame(exclusions).to_csv(output_dir / 'exclusions_log.csv', index=False)
    print("OK Clean data and exclusions saved")

    # Metadata summary
    metadata_summary = df_clean.groupby(['annotator', 'region', 'cell_type']).size().reset_index(name='count')
    metadata_summary.to_csv(output_dir / 'metadata_summary.csv', index=False)
    print("OK Metadata summary saved")

    # Statistical analysis
    stats_df = run_all_comparisons(df_clean)
    stats_df.to_csv(output_dir / 'statistical_results.csv', index=False)
    print("OK Statistical results saved")

    # Descriptive statistics
    metrics = ['volume_um3', 'surface_area_um2', 'lateral_spread_um', 'z_span_um',
               'ipl_depth_pct', 'sa_v_ratio', 'sphericity', 'aspect_ratio']

    desc_stats = []
    for (region, cell_type), group in df_clean.groupby(['region', 'cell_type']):
        for metric in metrics:
            vals = group[metric].dropna()
            desc_stats.append({
                'region': region,
                'cell_type': cell_type,
                'metric': metric,
                'n': len(vals),
                'mean': vals.mean(),
                'sd': vals.std(),
                'median': vals.median(),
                'q25': vals.quantile(0.25),
                'q75': vals.quantile(0.75),
                'min': vals.min(),
                'max': vals.max()
            })

    desc_stats_df = pd.DataFrame(desc_stats)
    desc_stats_df.to_csv(output_dir / 'descriptive_stats.csv', index=False)
    print("OK Descriptive statistics saved")

    # Generate figures
    print("\n=== GENERATING FIGURES ===")
    create_figure1_stratification(df_clean, output_dir)
    create_figure2_regional_comparison(df_clean, stats_df, output_dir)
    create_figure3_interrater(df_clean, output_dir)
    create_figure4_shape_complexity(df_clean, output_dir)
    create_figure5_effect_heatmap(stats_df, output_dir)

    # Generate report
    print("\n=== GENERATING REPORT ===")
    generate_report(df, df_clean, exclusions, stats_df, calibration_info, output_dir)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"\nResults saved to: {output_dir}")
    print(f"Total files created: {len(list(output_dir.glob('*')))}")
    print("\nKey outputs:")
    print(f"  - {len(df_clean)} cells analyzed (from {len(df)} loaded)")
    print(f"  - {len(stats_df)} statistical comparisons")
    print(f"  - 5 figures (PNG + SVG)")
    print(f"  - Comprehensive analysis report")
    print()

if __name__ == '__main__':
    main()
