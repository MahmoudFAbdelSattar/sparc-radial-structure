# CFT Concentration Law – Observational Analysis

This repository contains the numerical code that accompanies the paper

> **“Unifying the Sectors: How the Central Equation of Causal Field Theory Succeeds on the Galactic Plateau”**  
> M. F. Abdel‑Sattar, 2026

The code extracts the equilibrium saturation value \(\phi_{\text{end}}\) of the causal field from the SPARC galaxy rotation curves, computes the **Concentration Law**
\[
\phi_{\text{end}} = 5.3319 - 0.4935\,\log_{10}(M_{\text{tot}}/M_\odot),
\]
and produces the publication‑ready figures that appear in the manuscript.

---

## Description

The observational analysis proceeds in three steps:

1. **Load the SPARC data** – The publicly available SPARC database (Lelli et al. 2016) is downloaded and the galaxies that have been independently classified as *Saturated (Plateau)* in Paper I are selected.
2. **Extract \(\phi_{\text{end}}\)** – For each galaxy the observed causal field \(\phi_{\text{obs}}(r)=\ln(v_{\text{obs}}/v_{\text{bar}})\) is computed. The median value in the outer plateau region (\(r \ge 5\) kpc) defines \(\phi_{\text{end}}\). The total baryonic mass \(M_{\text{tot}}\) and the mean surface density \(\Sigma\) are evaluated at the outermost reliable point.
3. **Fit the Concentration Law** – A linear regression of \(\phi_{\text{end}}\) against \(\log_{10}M_{\text{tot}}\) is performed, and the statistical significance of the correlation is assessed. The same analysis is repeated for the surface density to show that mass is the primary driver.

All numerical values, figures, and CSV files reported in the paper are reproduced by executing the single script `concentration_law_analysis.py`.

---

## Requirements

- Python 3.8 or later
- NumPy
- SciPy
- Matplotlib
- pandas
- requests (included in the standard library)

Install the dependencies with

```bash
pip install numpy scipy matplotlib pandas
```

---

Usage

Clone the repository and run the analysis script:

```bash
git clone https://github.com/msattar/cft-concentration-law.git
cd cft-concentration-law
python concentration_law_analysis.py
```

The script automatically downloads the SPARC data and the regime classification from Paper I; no manual data preparation is needed. Execution time is about one minute on a modern laptop.

---

Output

After a successful run the following files are created in the working directory:

File Content
concentration_law_results.csv \phi_{\text{end}}, \log M_{\text{tot}}, \log\Sigma, \phi_{\text{std}} and number of outer points for each galaxy
fig1_example_phi.pdf Causal field \phi(r) for the example galaxy NGC 3198
fig2_phi_end_distribution.pdf Histogram of \phi_{\text{end}} across the 36 plateau galaxies
fig3_concentration_law.pdf Left: \phi_{\text{end}} vs. \log M_{\text{tot}} with the best‑fit line; Right: \phi_{\text{end}} vs. \log\Sigma
fig4_residuals.pdf Residuals of the linear Concentration Law fit

The terminal output contains the full statistical summary (correlation coefficients, p-values, R^2, and the coefficients of the linear and multiple regressions).

---

Reference

If you use this code in your research, please cite the accompanying paper:

M. F. Abdel‑Sattar, “Unifying the Sectors: How the Central Equation of Causal Field Theory Succeeds on the Galactic Plateau,” in preparation (2026).

---

Author

Mahmoud F. Abdel‑Sattar
Department of Astronomy and Meteorology
Faculty of Science, Al‑Azhar University, Cairo, Egypt