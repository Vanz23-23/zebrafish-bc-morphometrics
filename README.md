# Regional Morphological Specialisation of Monostratifying Bipolar Cells in the Larval Zebrafish Retina

**Author:** Vanz Labitad  
**Institution:** University of Sussex, School of Life Sciences  
**Degree:** MSci Neuroscience  
**Supervisor:** Katarina Moravkova

---

This repository contains the analysis code, morphometric data, & figures for my dissertation project investigating regional morphological differences between monostratifying bipolar cell subtypes (S2 and S4) across the dorsal and ventral retina of larval zebrafish.

Cells were manually reconstructed from serial-section electron microscopy data using WebKnossos. 3D morphometric features included the volume, surface area, lateral spread, IPL depth, convex hull, and shape indices these were extracted from STL mesh exports. Statistical comparisons use Mann-Whitney U tests and Cohen's *d* effect sizes throughout.

<div align="center">
<img src="assets/Ventral Vs Dorsal comparison.png"
     alt="3D mesh reconstructions of S2 and S4 monostratifying bipolar cells — dorsal vs ventral retinal regions"
     width="100%"/>
<br/>
<em>3D mesh reconstructions of S2 and S4 monostratifying bipolar cells — dorsal vs ventral retinal regions. Dataset: ZF_retina_opl_ipl_model_v5.</em>
</div>

---

## Repository Structure

```
zebrafish-bc-morphometrics/
├── analysis/
│   └── morphometric_analysis.py     # Main analysis pipeline
├── data/
│   ├── morphometrics_clean.csv      # 185 rows, 24 columns — post-QC dataset
│   ├── morphometrics_raw.csv        # 222 rows, 23 columns — pre-QC dataset
│   ├── exclusions_log.csv           # 37 excluded cells with reasons
│   ├── metadata_summary.csv         # 8 rows — cell counts per group
│   ├── descriptive_stats.csv        # 32 rows — summary statistics
│   └── statistical_results.csv      # 32 rows — Mann-Whitney U + Cohen's d
├── figures/
│   ├── figure1_stratification.png   # IPL depth by cell type and region
│   ├── figure2_regional_comparison.png  # Dorsal vs ventral morphology
│   ├── figure3_interrater.png       # Inter-rater reliability (VL vs KM)
│   ├── figure4_shape_complexity.png # Volume–surface area relationship
│   └── figure5_effect_heatmap.png  # Cohen's d effect size heatmap
└── assets/
    ├── Ventral Vs Dorsal comparison.png
    ├── S2 and S4 skeleotisations in Dorsal region.png
    ├── Dorsal render.png
    └── S4 Ventral skeletonisations.png
```

---

## Dataset

| File | Rows | Columns | Description |
|---|---|---|---|
| `morphometrics_raw.csv` | 222 | 23 | Pre-QC measurements from all annotated cells |
| `morphometrics_clean.csv` | 185 | 24 | Post-QC dataset used in all analyses |
| `exclusions_log.csv` | 37 | — | Excluded cells: volume outliers (>3 SD) or implausible IPL depth |
| `metadata_summary.csv` | 8 | — | Cell counts per annotator × cell type × region group |
| `descriptive_stats.csv` | 32 | — | Mean, SD, median, IQR per group and morphometric variable |
| `statistical_results.csv` | 32 | — | Mann-Whitney U statistics, *p*-values, Cohen's *d* |

**37 cells excluded** from the clean dataset: cells were removed if their volume fell outside ±3 SD of the group mean, or if their IPL depth fell outside the plausible 0–100% normalised range.

---

## Dual-Annotator Design

All cells were annotated in WebKnossos by two independent annotators to assess inter-rater reliability:

| Annotator | Role |
|---|---|
| **VL** — Vanz Labitad | Primary annotator |
| **KM** — Katarina Moravkova | Independent annotator (supervisor) |

Inter-rater agreement was assessed across all morphometric variables. Results confirmed high reliability (all Cohen's *d* < 0.10, all *p* > 0.05), supporting the validity of the primary annotation set.

---

## Analysis Pipeline

1. **STL mesh export** — Reconstructed cell meshes exported from WebKnossos in STL format
2. **3D morphometric extraction** — Volume, surface area, lateral spread, IPL depth, convex hull volume, and shape indices computed per cell
3. **IPL stratification normalisation** — IPL depth normalised to 0–100% range using 2nd–98th percentile z-axis calibration per dataset
4. **Quality control** — 37 cells excluded based on volume and IPL depth criteria (see the `exclusions_log.csv`)
5. **Statistical testing** — Mann-Whitney U tests with Cohen's *d* effect sizes for all pairwise comparisons
6. **Comparisons performed:**
   - S2 vs S4 (collapsed across region)
   - Dorsal vs Ventral (collapsed across cell type)
   - Dorsal S4 vs Ventral S4 (primary comparison of interest)

---

## Key Results

| Comparison | Variable | Cohen's *d* | *p*-value |
|---|---|---|---|
| Dorsal S4 vs Ventral S4 | Terminal width | −3.12 | < 0.001 |
| Dorsal S4 vs Ventral S4 | IPL depth | −4.63 | < 0.001 |
| S2 vs S4 | Surface area | 0.92 | < 0.001 |
| Inter-rater (VL vs KM) | All variables | < 0.10 | > 0.05 |

Ventral S4 cells show significantly wider terminals than their dorsal counterparts (Cohen's *d* = −3.12, *p* < 0.001). This is consistent with the functional demands of the ventral retina, which samples the bright overhead visual field and is expected to favour increased photon catch and spatial summation. The large effect size for IPL depth (Cohen's *d* = −4.63) confirms that S4 stratification position differs substantially between dorsal and ventral regions, suggesting topographic functional specialisation at the level of individual bipolar cell subtypes.

<table><tr>
<td width="33%" align="center">
<img src="assets/S2 and S4 skeleotisations in Dorsal region.png" width="100%"/>
<br/><em>Dorsal S2 and S4 skeletonisations</em></td>
<td width="33%" align="center">
<img src="assets/Dorsal render.png" width="100%"/>
<br/><em>Segmented dorsal region bipolar cells</em></td>
<td width="33%" align="center">
<img src="assets/S4 Ventral skeletonisations.png" width="100%"/>
<br/><em>Ventral S4 bipolar cell population</em></td>
</tr></table>

---

## Figures

| Figure | File | Description |
|---|---|---|
| Figure 1 | `figure1_stratification.png` | IPL depth distributions by cell type and region |
| Figure 2 | `figure2_regional_comparison.png` | Dorsal vs ventral morphometric comparison |
| Figure 3 | `figure3_interrater.png` | Inter-rater reliability: VL vs KM |
| Figure 4 | `figure4_shape_complexity.png` | Volume–surface area relationship |
| Figure 5 | `figure5_effect_heatmap.png` | Cohen's *d* effect size heatmap across all comparisons |

---

## Quickstart

**Requirements:** Python 3.8+

```bash
pip install pandas numpy scipy matplotlib seaborn
python analysis/morphometric_analysis.py
```

The script reads from `data/morphometrics_clean.csv` by default and writes outputs to `figures/` and `data/`.

---

## Data Access

Raw EM image stacks and STL mesh files are property to the University of Sussex & Baden Lab and are not included in this repository. This repository contains only morphometric measurements (CSV files) and the analysis code used to produce the reported results. Researchers wishing to access the raw data should contact the Baden Lab, University of Sussex.

---

## Limitations

- Sample sizes per subgroup are modest.
- Analyses are restricted to (6 days old)larval zebrafish; generalisation to adult or other species requires further study
- Morphometric extraction assumes mesh quality sufficient for accurate volume and surface area computation; cells failing QC thresholds were excluded
- Region boundaries (dorsal/ventral) were predefined can vary between Zebrafish.

---

## License

Code is released under the MIT License. These data files are made available for academic use. **Raw EM data remains the property of the University of Sussex.**

