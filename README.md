# Regional Morphological Specialisation of Monostratifying Bipolar Cells in the Larval Zebrafish Retina

**Author:** Vanz Labitad
**Institution:** University of Sussex, School of Life Sciences
**Degree:** MSci Neuroscience, 2025
**Supervisor:** Katarina Moravkova

## Overview

The larval zebrafish retina is topographically organised: the ventral retina samples the bright dorsal visual field (sky), while the dorsal retina samples the darker ventral field (ground). This project used manual annotation of serial-section electron microscopy data in webKnossos (dataset: `ZF_retina_opl_ipl_model_v5_segmentation_v1`) to compare the 3D morphology of S2 (OFF-pathway, inner plexiform layer stratum 2) and S4 (ON-pathway, stratum 4) monostratifying bipolar cells between dorsal and ventral retinal regions. The principal finding is that ventral S4 bipolar cells have significantly wider terminal than their dorsal counterparts (Cohen's *d* = -3.12, *p* < 0.001), consistent with regional functional specialisation driven by the ecological demands of processing bright overhead stimuli.

## Dataset

All data files are in `data/`:

| File | Rows | Description |
|------|------|-------------|
| `morphometrics_clean.csv` | 185 | Post-Quality Check morphometric measurements (24 columns including volume, surface area, lateral spread, IPL depth %, sphericity, aspect ratio) |
| `morphometrics_raw.csv` | 222 | Pre-Quality Check measurements (23 columns, no IPL depth normalisation) |
| `exclusions_log.csv` | 37 | Excluded cells with reasons: volume outliers (>3 SD from group mean) or implausible IPL depth (<5% or >95%) |
| `metadata_summary.csv` | 8 | Cell counts per annotator x region x cell type |
| `descriptive_stats.csv` | 32 | Summary statistics (n, mean, SD, median, IQR) per group and metric |
| `statistical_results.csv` | 32 | Mann-Whitney U test results with Cohen's *d* effect sizes for all comparisons |

Coordinates were extracted from 3D mesh reconstructions in nanometres and converted to micrometres for analysis.

## Dual Annotator Design

All cells were independently annotated by two trained annotators to assess inter-rater reliability:

- **VL** -- Vanz Labitad (primary annotator)
- **KM** -- Katarina Moravkova (independent annotator/Supervisor)

No significant differences were found between annotators across any morphometric measure (all *p* > 0.05, all Cohen's *d* < 0.10), confirming excellent inter-rater agreement. See `figure3_interrater.png` and the VL vs KM rows in `statistical_results.csv`.

## Analysis Pipeline

1. **STL mesh export** from webKnossos (binary format)
2. **3D morphometric extraction**: volume (signed tetrahedron method), surface area, lateral spread, IPL depth, convex hull metrics, shape indices
3. **IPL stratification normalisation**: empirical percentile-based z-axis calibration (2nd-98th percentile of z-centroids)
4. **Quality control**: 37 cells excluded (see `exclusions_log.csv`)
5. **Statistical analysis**: Mann-Whitney U tests (non-parametric, appropriate for non-normal distributions), Cohen's *d* effect sizes
6. **Comparisons**: S2 vs S4 (cell type validation), Dorsal vs Ventral (regional differences), Dorsal S4 vs Ventral S4 (key regional test)

## Reproducing the Analysis

```bash
pip install pandas numpy scipy matplotlib seaborn
python analysis/morphometric_analysis.py
```

Requires Python 3.8+. All input data is in `data/`. Output figures are saved to `figures/`.

## Key Results

| Comparison | Metric | Cohen's *d* | *p*-value | Effect |
|------------|--------|-------------|-----------|--------|
| Dorsal S4 vs Ventral S4 | Lateral spread | -3.12 | < 0.001 | Large |
| Dorsal S4 vs Ventral S4 | IPL depth | -4.63 | < 0.001 | Large |
| S2 vs S4 | Surface area | 0.92 | < 0.001 | Large |
| VL vs KM | All metrics | < 0.10 | > 0.05 | None |

Morphological divergence between dorsal and ventral S4 cells supports retinal topographic tuning to visual ecology. S2 vs S4 differences in volume and surface area validate the cell type classification.

## Data Access

Raw EM image stacks and STL mesh files are proprietary to the University of Sussex and are not included in this repository. This repository contains only derived morphometric measurements and analysis code. Contact the Author for data access enquiries.

## Figures

| Figure | Description |
|--------|-------------|
| `figure1_stratification.png` | IPL depth distribution by cell type and region |
| `figure2_regional_comparison.png` | Regional morphology comparison (Dorsal vs Ventral) |
| `figure3_interrater.png` | Inter-rater reliability (VL vs KM) |
| `figure4_shape_complexity.png` | Volume-surface area relationship by region |
| `figure5_effect_heatmap.png` | Cohen's *d* effect size heatmap |

## License

Code: MIT License. Data: shared for academic transparency and reproducibility. For VL Dissertation!
