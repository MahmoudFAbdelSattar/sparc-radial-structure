```markdown
# median_discrepancy_mass

**Radial Structure in the Baryonic Mass Discrepancy of SPARC Spirals:
Beyond the Tully–Fisher and Mass–Radius Relations**

[![Python 3.13+](https://img.shields.io/badge/python-3.13%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## Overview

This repository contains the complete, reproducible analysis pipeline for the
manuscript:

> **"Radial Structure in the Baryonic Mass Discrepancy of SPARC Spirals:
> Beyond the Tully–Fisher and Mass–Radius Relations"**
> *Submitted to The Astrophysical Journal*

The pipeline analyzes the logarithmic baryonic mass discrepancy,

$$\phi(r) = \ln\!\left[\frac{v_{\mathrm{obs}}(r)}{v_{\mathrm{bar}}(r)}\right],$$

in the saturated kinematic regime ($r \geq 5$ kpc) of $N = 83$ SPARC disc
galaxies with $M_{\mathrm{bar}} < 10^{11}\,M_\odot$, and separates it into
three components:

1. **A median discrepancy–mass relation**
   $$\phi_{\mathrm{end}} = 4.04 - 0.351\,\log_{10}\!\left(\frac{M_{\mathrm{bar}}}{M_\odot}\right)$$
   with $r = -0.797$ and $R^2 = 0.635$.

2. **An algebraic identity** showing that the *slope* of the endpoint value
   $\phi(R_{\max})$ is fully accounted for by the combination of the
   Baryonic Tully–Fisher Relation (BTFR) and the empirical mass–radius
   relation $R_{\max} \propto M_{\mathrm{bar}}^{p}$. Approximately 77% of
   the observed slope is therefore known physics.

3. **A new observable — the radial structure**
   $$\Delta\phi \equiv \phi_{\mathrm{end}} - \phi(R_{\max}),$$
   which correlates strongly with galaxy size at fixed baryonic mass
   ($\rho_p = -0.42$, $p < 10^{-4}$) and constitutes the remaining 23% of
   the signal.

---

## Key results at a glance

| Quantity | Value | Source |
|---|---|---|
| Sample size | $N = 83$ galaxies, 923 radial points | §2 |
| Median relation | $\phi_{\mathrm{end}} = 4.04 - 0.351 \log_{10} M_{\mathrm{bar}}$ | §4.1 |
| Pearson $r$ | $-0.797$ | §4.1 |
| $R^2$ | $0.635$ | §4.1 |
| BTFR slope $\alpha$ | $0.2304 \pm 0.0111$ | §4.2 |
| Mass–radius slope $p$ | $0.3038 \pm 0.0262$ | §4.2 |
| BTFR+R(M) predicted slope | $-0.2709$ | §4.2 |
| Observed endpoint slope | $-0.2709$ | §4.2 |
| Residual (observed − predicted) | $< 0.0001$ | §4.2 |
| Radial structure $\langle\Delta\phi\rangle$ | $-0.1026$ | §4.3 |
| $\rho_p(\Delta\phi, \log R_{\max} \mid \log M_{\mathrm{bar}})$ | $-0.419$ ($p < 10^{-4}$) | §4.3 |
| Tertile stability $r$ | $-0.574, -0.375, -0.455$ | §4.4 |
| Total systematic uncertainty | $\pm 0.12$ | §5 |

---

## Installation

### Requirements

- Python **3.13.18** (tested)
- NumPy 2.4.6, Pandas 2.2.3, SciPy 1.16.3
- statsmodels 0.15.0, scikit-learn 1.6.1
- matplotlib, requests

### Install dependencies

```bash
pip install numpy==2.4.6 pandas==2.2.3 scipy==1.16.3 \
            statsmodels==0.15.0 scikit-learn==1.6.1 \
            matplotlib requests
```

Or use the provided requirements file:

```bash
pip install -r requirements.txt
```

Clone the repository

```bash
git clone https://github.com/[user]/median-discrepancy-mass.git
cd median-discrepancy-mass
```

---

Quick start

Run the complete pipeline

```bash
python median_discrepancy_mass.py
```

Expected runtime: ~30 seconds on a modern laptop (after the first
download of the SPARC catalogue, which takes 1–2 seconds).

Run within Python

```python
from median_discrepancy_mass import run_pipeline, Config

results = run_pipeline(Config(), verbose=True)

# Access key results
print(results["main_median"]["equation"])
print(results["btfr_consistency"]["predicted_slope"])
print(results["radial_structure"]["partial_correlations"]["log_R_max"])
```

Command-line options

```bash
python median_discrepancy_mass.py --help
python median_discrepancy_mass.py --force-refresh     # re-download SPARC
python median_discrepancy_mass.py --quiet             # suppress logging
python median_discrepancy_mass.py --output-dir ./out  # custom output path
python median_discrepancy_mass.py --seed 42           # master random seed
```

---

Output structure

All artefacts are written to median_discrepancy_mass_output/ by default:

```
median_discrepancy_mass_output/
├── sample_final.csv                    # 83 galaxies × 26 columns
├── table1_sample.tex                   # LaTeX tables 1–11
├── table2_correlations.tex
├── table3_partials.tex
├── table4_models.tex
├── table5_ml_sensitivity.tex
├── table6_cv.tex
├── table7_selection.tex
├── table8_r_sat.tex
├── table9_btfr_consistency.tex
├── table10_radial_structure.tex
├── table11_disappearing.tex
├── table2_correlations.csv             # Raw data for each table
├── table3_partials.csv
├── table4_models.csv
├── table5_ml_sensitivity.csv
├── table6_cv.csv
├── table7_selection.csv
├── table8_r_sat.csv
├── table9_btfr_consistency.csv
├── table10_radial_structure.csv
├── table11_disappearing.csv
├── rejection_log.csv                   # Per-galaxy rejection reasons
├── run_manifest.json                   # Full reproducibility manifest
├── sparc_massmodels.pkl                # Cached SPARC catalogue
├── figures/
│   ├── figure1_main_dual.pdf           # Median vs endpoint panels
│   ├── figure2_btfr_prediction.pdf     # BTFR+R(M) consistency
│   ├── figure3_radial_structure.pdf    # Δφ vs R_max
│   ├── figure4_tertile_stability.pdf   # Stability across mass tertiles
│   ├── figure5_residuals.pdf           # 4-panel residual diagnostics
│   ├── figure6_ml_sensitivity.pdf      # Robustness to Υ_disk
│   └── figure7_disappearing.pdf        # Galaxies that lose φ>0 at Υ=1.4
└── tables/
    └── (same .tex files as above)
```

run_manifest.json

The manifest records everything needed to reproduce the results:

· Versions of Python, NumPy, Pandas, SciPy, statsmodels, scikit-learn
· All Config parameters
· Sample fingerprint (N, saturated points, ranges, SHA-256 hash)
· Main results (slope, intercept, r, R², p-value, bootstrap CI)
· BTFR consistency (α, p, predicted, observed, residual)
· Radial structure summary (mean, median, partial correlation, tertile r)
· Disappearing galaxies count

---

Configuration

All analysis parameters are declared in Config (a Python dataclass):

Parameter Default Description
UPSILON_DISK 0.5 Stellar M/L of the disc at 3.6 µm
UPSILON_BULGE 0.7 Stellar M/L of the bulge at 3.6 µm
R_SAT_KPC 5.0 Saturated-regime threshold
V_BAR_MIN_KMS 3.0 Minimum baryonic velocity per point
PHI_MAX 1.5 Per-point outlier threshold
N_SAT_MIN 4 Minimum saturated points per galaxy
M_BAR_MAX_MSUN $10^{11}$ Maximum baryonic mass per galaxy
BOOTSTRAP_N 5000 Bootstrap resamples
CV_SPLITS 10 Cross-validation folds
RANDOM_SEED 42 Master random seed

To customize, create a modified Config:

```python
from median_discrepancy_mass import run_pipeline, Config

cfg = Config(
    UPSILON_DISK=0.7,
    M_BAR_MAX_MSUN=2e11,
    OUTPUT_DIR="custom_output",
)
results = run_pipeline(cfg)
```

---

Data source

The SPARC catalogue is downloaded automatically from:

· Primary: 
· Fallback: 

The catalogue is cached locally (sparc_massmodels.pkl) after the first
download. To force a fresh download, use --force-refresh or set
FORCE_REFRESH=True in Config.

Reference: Lelli, F., McGaugh, S. S., & Schombert, J. M. 2016,
AJ, 152, 157.

---

Reproducibility

Guarantees

· Single random seed (42) governs every stochastic procedure
  (bootstrap, cross-validation, tertile partitioning).
· Deterministic fits: OLS, AIC/BIC, residuals are all deterministic.
· Cached data: the SPARC catalogue is archived with the pipeline.
· Byte-identical outputs: running the pipeline with the same Config
  and the same dependency versions yields the same figures, tables, and
  numerical results.

Verify the installation

```python
from median_discrepancy_mass import run_pipeline, Config

results = run_pipeline(Config())

assert results["main_median"]["N"] == 83
assert abs(results["main_median"]["slope"] + 0.351) < 1e-3
assert abs(results["btfr_consistency"]["residual"]) < 1e-15
assert abs(
    results["radial_structure"]["partial_correlations"]["log_R_max"]["rho"] + 0.419
) < 1e-2

print("✓ All checks passed.")
```

Software environment

The exact versions used at the time of publication:

```
Python   3.13.18
NumPy    2.4.6
Pandas   2.2.3
SciPy    1.16.3
statsmodels 0.15.0
scikit-learn 1.6.1
matplotlib 3.x
```

Any deviation may produce small numerical differences, but the qualitative
results (slopes, correlations, partial correlations) are robust to
version changes.

---

Citation

If you use this code or the derived quantities in your research, please
cite:

```bibtex
@article{AbdelSattar2026radial,
  author  = {Abdel-Sattar, Mahmoud F.},
  title   = {Radial Structure in the Baryonic Mass Discrepancy of
             {SPARC} Spirals: Beyond the {T}ully--{F}isher and
             Mass--Radius Relations},
  journal = {The Astrophysical Journal},
  year    = {2026},
  note    = {Submitted},
  doi     = {10.5281/zenodo.XXXXXXX}
}
```

The analysis pipeline is archived at Zenodo:

DOI: 10.5281/zenodo.XXXXXXX

---

Repository structure

```
median-discrepancy-mass/
├── median_discrepancy_mass.py       # Complete analysis pipeline
├── README.md                        # This file
├── LICENSE                          # MIT License
├── requirements.txt                 # Pinned dependencies
├── CITATION.cff                     # Citation File Format
├── .gitignore                       # Ignore output/cache directories
├── tests/
│   └── test_pipeline.py             # Smoke tests and reproducibility checks
└── manuscript/                      # (Optional) LaTeX source
    ├── main.tex
    ├── references.bib
    ├── figures/
    └── tables/
```

---

Backward compatibility

The module was previously named concentration_law.py. If you have code
that imports the old name, update it as follows:

```python
# Before
from concentration_law import run_pipeline, Config, ConcentrationSample

# After
from median_discrepancy_mass import run_pipeline, Config, DiscrepancySample
```

Class renaming:

· ConcentrationSample → DiscrepancySample
· ConcentrationAnalyzer → DiscrepancyAnalyzer
· PaperVisualizer → unchanged
· SPARCFetcher → unchanged

Output directory: concentration_law_output_v4/ →
median_discrepancy_mass_output/

The LEGACY_CACHE_DIR = "concentration_law_output" field is preserved
only to reuse the cached SPARC catalogue (avoiding a re-download); no
analysis depends on the legacy path.

---

Terminology note

The term "concentration law" used in earlier drafts was replaced
throughout the manuscript and code with "median discrepancy–mass relation",
because the analysis shows that approximately 77% of the observed slope is
already encoded in the BTFR + R(M) identity. This is a purely descriptive
rename; no numerical results changed.

---

License

This project is licensed under the MIT License — see the LICENSE
file for details.

The SPARC catalogue is distributed under the terms of its original
publication (Lelli et al. 2016) and is subject to the policies of the
SPARC collaboration.

---

Contact

Author: Mahmoud F. Abdel-Sattar
Affiliation: Department of Physics, Al-Azhar University, Cairo, Egypt
Email: m.f.abdel-sattar@azhar.edu.eg

For questions about the pipeline, reproducibility issues, or potential
collaborations, please open an issue or contact the author directly.

---

Acknowledgements

This research made use of:

· The SPARC (Spitzer Photometry and Accurate Rotation Curves) catalogue
  (Lelli, McGaugh & Schombert 2016)
· NASA's Astrophysics Data System (ADS)
· VizieR catalogue service at CDS, Strasbourg
· Open-source Python libraries: NumPy, Pandas, SciPy, statsmodels,
  scikit-learn, matplotlib
· Zenodo for archiving the analysis pipeline

---

Changelog

v4.1.0 (current)

· Renamed module: concentration_law.py → median_discrepancy_mass.py
· Renamed classes: ConcentrationSample → DiscrepancySample,
  ConcentrationAnalyzer → DiscrepancyAnalyzer
· Renamed output directory: concentration_law_output_v4 →
  median_discrepancy_mass_output
· Updated CLI, User-Agent, and plot labels to match the manuscript's
  terminology
· No numerical results changed

v4.0.1

· Fixed LaTeX column count in table10_radial_structure (6 cells)
· Added guard against qcut failures in figure4_tertile_stability
· Added R² and N to the text box in figure2_btfr_prediction
· Replaced cryptic label in table10 with the word "partial"
· Added R² to both panels' legends in figure1_main_dual
· Refactored SPARCFetcher.load() into helper _load_from()
· Documented N difference vs v3.0.1 in table8_r_sat caption

v4.0.0

· Major reframing: BTFR + R(M) identity, radial structure $\Delta\phi$
· Added compute_btfr_rm_consistency() and compute_radial_structure()
· Added tertile stability analysis
· Added disappearing-galaxies test at Υ=1.4
· Retracted earlier claims of "disk dominance" and "gas inversion"
  (both were consequences of the BTFR + R(M) identity, not independent
  effects)

---

References

Key references for the underlying data and methods (full list in the
manuscript):

· Lelli, F., McGaugh, S. S., & Schombert, J. M. 2016, AJ, 152, 157 — SPARC
· McGaugh, S. S., Lelli, F., & Schombert, J. M. 2016, PRL, 117, 201101 — RAR
· McGaugh, S. S., Schombert, J. M., Bothun, G. D., & de Blok, W. J. G. 2000, ApJ, 533, L99 — BTFR
· Lelli, F., McGaugh, S. S., Schombert, J. M., & Pawlowski, M. S. 2017, ApJ, 836, 152 — RAR (detailed)
· McGaugh, S. S., & Schombert, J. M. 2014, AJ, 148, 77 — Color–M/L calibration

---

Contributing

Contributions are welcome. Please:

1. Fork the repository
2. Create a feature branch (git checkout -b feature/my-improvement)
3. Ensure tests pass (pytest tests/)
4. Ensure the pipeline reproduces the published numbers
5. Submit a pull request with a clear description

Important: any change that alters the numerical results must be
justified in the pull request, with a corresponding update to the
manuscript and/or the reproducibility manifest.

---

Disclaimer

This pipeline is provided for scientific reproducibility purposes. The
results are specific to the SPARC sample and the sample-selection cuts
described in the manuscript. Applying the pipeline to other datasets
requires careful re-validation of the cuts and configurations.

```
