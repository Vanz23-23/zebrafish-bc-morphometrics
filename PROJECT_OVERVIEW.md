<div align="center">

<!-- HERO SECTION -->

# 🧬 Zebrafish Retinal Bipolar Cell Morphometrics

### 3D Reconstruction · Quantitative Analysis · Computational Neuroscience

<br/>

<img src="assets/Ventral Vs Dorsal comparison.png"
     alt="3D mesh reconstructions of bipolar cells in dorsal and ventral zebrafish retina"
     width="90%"/>

<br/>

<sub>3D electron microscopy reconstructions of retinal bipolar cells — dorsal (left) vs ventral (right) regions</sub>

<br/><br/>

![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-Scientific_Computing-013243?style=for-the-badge&logo=numpy&logoColor=white)
![Pandas](https://img.shields.io/badge/Pandas-Data_Analysis-150458?style=for-the-badge&logo=pandas&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-Statistics-8CAAE6?style=for-the-badge&logo=scipy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-Visualisation-11557C?style=for-the-badge)
![WebKnossos](https://img.shields.io/badge/WebKnossos-3D_Annotation-FF6B6B?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)

</div>

<br/>

---

<br/>

## 🔬 The Problem

> **How does the shape of individual neurons adapt to the visual demands of different parts of the eye?**

The zebrafish retina contains specialised nerve cells called **bipolar cells** that relay visual information from photoreceptors to the brain. Different regions of the eye face dramatically different visual environments — the **ventral retina** looks upward at the bright sky, while the **dorsal retina** looks downward at darker terrain.

This project investigates whether these environmental pressures produce **measurable morphological differences** at the level of individual cell subtypes.

<br/>

---

<br/>

## 🛠️ What I Built

<table>
<tr>
<td width="50%">

### End-to-End Analysis Pipeline

A **1,300+ line Python pipeline** that takes raw 3D cell reconstructions and produces publication-quality statistical analysis and figures — fully automated and reproducible.

**Pipeline stages:**
- 🔹 Binary STL mesh parsing (custom parser)
- 🔹 8 morphometric feature extraction (volume, surface area, lateral spread, IPL depth, convex hull, shape indices)
- 🔹 Automated IPL depth normalisation with z-axis correction
- 🔹 Statistical QC filtering (±3 SD, depth validation)
- 🔹 Mann-Whitney U tests with Cohen's *d* effect sizes
- 🔹 5 publication-quality figure generators

</td>
<td width="50%">

### Technical Highlights

| Skill | Application |
|---|---|
| **3D Geometry** | Signed tetrahedron volume, convex hull computation, mesh surface area |
| **Statistics** | Non-parametric tests, effect sizes (Cohen's *d*, rank-biserial *r*) |
| **Data Engineering** | Automated QC pipeline, reproducible CSV outputs |
| **Visualisation** | Violin + swarm plots, heatmaps, regression with CI bands |
| **Scientific Rigour** | Dual-annotator validation, sensitivity analysis, documented limitations |
| **Domain Knowledge** | Retinal neuroscience, EM reconstruction, IPL stratification |

</td>
</tr>
</table>

<br/>

---

<br/>

## 📊 Key Results at a Glance

<div align="center">

```
                    ┌──────────────────────────────────────────────────┐
                    │         DORSAL S4  vs  VENTRAL S4               │
                    │                                                  │
                    │   Terminal Width    d = -3.12 ***   p < 0.001    │
                    │   IPL Depth        d = -4.63 ***   p < 0.001    │
                    │   Aspect Ratio     d =  2.35 ***   p < 0.001    │
                    │                                                  │
                    │   *** = highly significant (p < 0.001)           │
                    │   Cohen's d > 0.8 = large effect                 │
                    └──────────────────────────────────────────────────┘
```

</div>

<br/>

> **Key Finding:** Ventral S4 bipolar cells have **3× wider terminal arbors** than their dorsal counterparts (*d* = −3.12), consistent with functional adaptation to bright overhead illumination. Inter-rater reliability confirmed all *d* < 0.10 — the measurements are robust.

<br/>

<table>
<tr>
<td width="33%" align="center">
<img src="assets/S2 and S4 skeleotisations in Dorsal region.png" width="100%"/>
<br/><sub><b>Dorsal Region</b> — S2 & S4 skeletonisations</sub>
</td>
<td width="33%" align="center">
<img src="assets/Dorsal render.png" width="100%"/>
<br/><sub><b>Segmented Cells</b> — 3D mesh renders</sub>
</td>
<td width="33%" align="center">
<img src="assets/S4 Ventral skeletonisations.png" width="100%"/>
<br/><sub><b>Ventral S4</b> — wider terminal spread</sub>
</td>
</tr>
</table>

<br/>

---

<br/>

## 📈 Figures

<table>
<tr>
<td width="50%" align="center">
<img src="figures/figure1_stratification.png" width="100%"/>
<br/><sub><b>Figure 1</b> — IPL stratification confirms S2/S4 classification</sub>
</td>
<td width="50%" align="center">
<img src="figures/figure2_regional_comparison.png" width="100%"/>
<br/><sub><b>Figure 2</b> — Regional morphology comparison</sub>
</td>
</tr>
<tr>
<td width="50%" align="center">
<img src="figures/figure4_shape_complexity.png" width="100%"/>
<br/><sub><b>Figure 4</b> — Volume vs surface area relationship</sub>
</td>
<td width="50%" align="center">
<img src="figures/figure5_effect_heatmap.png" width="100%"/>
<br/><sub><b>Figure 5</b> — Effect size heatmap across all comparisons</sub>
</td>
</tr>
</table>

<div align="center">
<img src="figures/figure3_interrater.png" width="70%"/>
<br/><sub><b>Figure 3</b> — Inter-rater reliability: near-perfect agreement between annotators</sub>
</div>

<br/>

---

<br/>

## 🧪 Data & Methodology

<div align="center">

```mermaid
flowchart LR
    A["🔬 Serial-Section EM\n(electron microscopy)"] --> B["✏️ Manual 3D Tracing\n(WebKnossos)"]
    B --> C["📦 STL Mesh Export\n(222 cells)"]
    C --> D["📐 Morphometric\nExtraction"]
    D --> E["🧹 Quality Control\n(185 cells remain)"]
    E --> F["📊 Statistical\nAnalysis"]
    F --> G["📈 Publication\nFigures"]
```

</div>

<br/>

| Stage | Detail |
|---|---|
| **Source Data** | Serial-section electron microscopy (University of Sussex, Baden Lab) |
| **Annotation** | 222 cells manually traced in WebKnossos by 2 independent annotators |
| **QC Filtering** | 37 cells excluded (volume outliers > 3 SD, implausible IPL depths) |
| **Final Dataset** | 185 cells · 24 features · 4 subgroups (S2/S4 × Dorsal/Ventral) |
| **Statistics** | Mann-Whitney U, Cohen's *d*, rank-biserial *r* · 32 comparisons |
| **Validation** | Dual-annotator inter-rater reliability (all *d* < 0.10, all *p* > 0.05) |
| **Sensitivity** | QC-excluded vs full dataset comparison confirms findings hold |

<br/>

---

<br/>

## 🗂️ Repository Structure

```
zebrafish-bc-morphometrics/
│
├── 📄 README.md                        ← Technical documentation
├── 📄 PROJECT_OVERVIEW.md              ← You are here
│
├── 🐍 morphometrics_pipeline.py        ← Full pipeline (1,339 lines)
│
├── 📂 analysis/
│   └── morphometric_analysis.py        ← Standalone statistical analysis
│
├── 📂 sensitivity_analysis/
│   └── sensitivity_analysis.py         ← QC robustness check
│
├── 📂 data/
│   ├── morphometrics_clean.csv         ← 185 cells, 24 features (post-QC)
│   ├── morphometrics_raw.csv           ← 222 cells, 23 features (pre-QC)
│   ├── exclusions_log.csv              ← 37 excluded cells with reasons
│   ├── statistical_results.csv         ← 32 Mann-Whitney U test results
│   ├── descriptive_stats.csv           ← Summary statistics per group
│   └── metadata_summary.csv            ← Cell counts per subgroup
│
├── 📂 figures/                         ← 5 publication-quality figures (PNG + SVG)
└── 📂 assets/                          ← 3D reconstruction renders
```

<br/>

---

<br/>

## 💡 Skills Demonstrated

<table>
<tr>
<td align="center" width="20%">
<h3>🐍</h3>
<b>Python</b><br/>
<sub>1,300+ lines of production-quality scientific code with type hints, docstrings, and modular design</sub>
</td>
<td align="center" width="20%">
<h3>📐</h3>
<b>3D Computation</b><br/>
<sub>Custom STL parser, signed tetrahedron volumes, convex hull analysis, mesh operations</sub>
</td>
<td align="center" width="20%">
<h3>📊</h3>
<b>Statistics</b><br/>
<sub>Non-parametric hypothesis testing, effect size estimation, multiple comparison awareness</sub>
</td>
<td align="center" width="20%">
<h3>🎨</h3>
<b>Data Viz</b><br/>
<sub>Publication-quality figures with Matplotlib & Seaborn — violin plots, heatmaps, regression CI bands</sub>
</td>
<td align="center" width="20%">
<h3>🧠</h3>
<b>Neuroscience</b><br/>
<sub>Domain expertise in retinal circuitry, IPL stratification, and EM reconstruction methodology</sub>
</td>
</tr>
</table>

<br/>

---

<br/>

## ⚡ Quick Start

```bash
# Clone the repository
git clone https://github.com/Vanz23-23/zebrafish-bc-morphometrics.git
cd zebrafish-bc-morphometrics

# Install dependencies
pip install pandas numpy scipy matplotlib seaborn

# Run the analysis (reads from data/, outputs to figures/)
python analysis/morphometric_analysis.py

# Or regenerate all figures from the full pipeline
python morphometrics_pipeline.py --figures-only
```

<br/>

---

<br/>

<div align="center">

## 👤 Author

**Vanz Labitad**
<br/>
MSci Neuroscience · University of Sussex
<br/><br/>
Supervised by **Katarina Moravkova** · Baden Lab, School of Life Sciences

<br/>

---

<sub>Code released under MIT License · Data files for academic use · Raw EM data property of University of Sussex</sub>

</div>
