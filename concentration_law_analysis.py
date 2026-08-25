#!/usr/bin/env python3
"""
CFT Concentration Law – Observational Analysis
===============================================
This script extracts the saturation value phi_end as the mean of the
causal field phi(r) = ln(v_obs / v_bar) in the outer plateau region
(r >= R_START) for SPARC galaxies classified as 'Saturated (Plateau)'.

It then fits the Concentration Law:
    phi_end = a + b * log10(M_tot)

and produces four publication-quality figures for Paper III:
    1. Example of phi(r) and phi_end for NGC 3198
    2. Distribution of phi_end across the sample
    3. Concentration Law (phi_end vs logM) and weak correlation with Sigma
    4. Residuals of the linear fit

Author : Mahmoud F. Abdel‑Sattar
Date   : 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import requests
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# 1.  CONSTANTS
# =============================================================================
G = 4.302e-6          # kpc (km/s)^2 / Msun
R_START = 5.0         # kpc – start of the plateau region
R_MAX = 25.0          # kpc – maximum radius used
MIN_OUTER_PTS = 3     # minimum number of points in the outer region

# Saturation parameters (from Papers I & II)
PHI_MAX = 0.429       # maximum observed saturation value

# =============================================================================
# 2.  LOAD DATA
# =============================================================================
print("Loading SPARC rotation curves ...")
url_sparc = "https://astroweb.case.edu/SPARC/MassModels_Lelli2016c.mrt"
response = requests.get(url_sparc)
lines = response.text.strip().split('\n')

galaxy_data = {}
for line in lines:
    if not line.strip() or line.startswith(('URL', 'Title', '===')):
        continue
    if len(line) < 60:
        continue
    try:
        gal = line[0:11].strip()
        R   = float(line[20:26])
        Vobs= abs(float(line[27:33]))
        Vgas= abs(float(line[40:46])) if line[40:46].strip() else 0.0
        Vdsk= abs(float(line[47:53])) if line[47:53].strip() else 0.0
        Vbul= abs(float(line[54:60])) if line[54:60].strip() else 0.0
        if R > 0 and Vobs > 0:
            galaxy_data.setdefault(gal, []).append([R, Vobs, Vgas, Vdsk, Vbul])
    except (ValueError, IndexError):
        continue
print(f"Loaded {len(galaxy_data)} galaxies.")

# Load regime classification from Paper I
url_stats = ("https://raw.githubusercontent.com/MahmoudFAbdelSattar/"
             "causal-field-sparc/main/galaxy_statistics_20260612_042844.csv")
df_stats = pd.read_csv(url_stats)
plateau_gals = df_stats[df_stats['regime'] == 'Saturated (Plateau)']['galaxy'].tolist()
print(f"Plateau galaxies in catalogue: {len(plateau_gals)}")

# =============================================================================
# 3.  EXTRACT phi_end OBSERVATIONALLY
# =============================================================================
def analyze_galaxy(gal_name):
    """
    Extract phi_end as the median of phi_obs in the outer plateau region.
    Also compute total mass and mean surface density.
    """
    if gal_name not in galaxy_data:
        return None

    data = np.array(galaxy_data[gal_name])
    r_obs   = data[:, 0]
    v_obs   = data[:, 1]
    v_gas   = data[:, 2]
    v_disk  = data[:, 3]
    v_bulge = data[:, 4]
    v_bar   = np.sqrt(v_gas**2 + v_disk**2 + v_bulge**2)

    # Require Vobs > Vbar (physical regime where phi is positive)
    mask = (v_obs > v_bar) & (r_obs <= R_MAX)
    if np.sum(mask) < 5:
        return None
    r_obs, v_obs, v_bar = r_obs[mask], v_obs[mask], v_bar[mask]

    # Sort by radius
    idx = np.argsort(r_obs)
    r_obs, v_obs, v_bar = r_obs[idx], v_obs[idx], v_bar[idx]

    # Outer plateau region
    mask_outer = r_obs >= R_START
    if np.sum(mask_outer) < MIN_OUTER_PTS:
        return None

    r_outer = r_obs[mask_outer]
    v_obs_outer = v_obs[mask_outer]
    v_bar_outer = v_bar[mask_outer]

    # Causal field in the outer region
    phi_obs_outer = np.log(v_obs_outer / v_bar_outer)

    # Use median as robust estimator of phi_end
    phi_end = np.median(phi_obs_outer)
    phi_std = np.std(phi_obs_outer)

    # Total mass at the outermost point
    M_tot = v_bar_outer[-1]**2 * r_outer[-1] / G
    Sigma = M_tot / (np.pi * r_outer[-1]**2)

    return {
        "gal": gal_name,
        "phi_end": phi_end,
        "phi_std": phi_std,
        "logM": np.log10(M_tot),
        "logSigma": np.log10(Sigma),
        "n_outer": len(r_outer),
        "r_outer": r_outer,
        "phi_obs_outer": phi_obs_outer
    }

# =============================================================================
# 4.  RUN ANALYSIS ON PLATEAU GALAXIES
# =============================================================================
print("\nExtracting phi_end from observed rotation curves ...")
results = []
for gal in plateau_gals:
    res = analyze_galaxy(gal)
    if res is not None and res["phi_end"] > 0:
        results.append(res)

df = pd.DataFrame(results)
n_success = len(df)
print(f"Successfully extracted phi_end for {n_success} galaxies")

if n_success < 5:
    print("Not enough galaxies. Exiting.")
    exit()

# =============================================================================
# 5.  STATISTICAL ANALYSIS
# =============================================================================
phi_end = df["phi_end"].values
logM    = df["logM"].values
logS    = df["logSigma"].values

# Pearson correlations
r_M, p_M = stats.pearsonr(logM, phi_end)
r_S, p_S = stats.pearsonr(logS, phi_end)
print(f"\nCorrelation φ_end vs logM  : r = {r_M:.4f}, p = {p_M:.4e}")
print(f"Correlation φ_end vs logΣ  : r = {r_S:.4f}, p = {p_S:.4e}")

# Linear regression: phi_end = a + b * logM
slope, intercept, r_value, p_value, std_err = stats.linregress(logM, phi_end)
print(f"\nLinear regression (Concentration Law):")
print(f"phi_end = {intercept:.4f} + {slope:.4f} * logM")
print(f"R² = {r_value**2:.4f}, p = {p_value:.4e}")

# Multiple regression with Sigma
X = np.column_stack([np.ones(len(df)), logM, logS])
coeff, _, _, _ = np.linalg.lstsq(X, phi_end, rcond=None)
A, B, C = coeff
phi_pred = X @ coeff
SS_res = np.sum((phi_end - phi_pred)**2)
SS_tot = np.sum((phi_end - np.mean(phi_end))**2)
R2_mult = 1.0 - SS_res / SS_tot
print(f"\nMultiple regression: phi_end = {A:.4f} + {B:.4f}·logM + {C:.4f}·logΣ")
print(f"R² = {R2_mult:.4f}")

# =============================================================================
# 6.  PUBLICATION‑QUALITY FIGURES
# =============================================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "figure.dpi": 150
})

# -----------------------------------------------------------------
# Figure 1: Example galaxy (NGC 3198) – phi(r) and phi_end
# -----------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(8, 6))
example_gal = "NGC3198"
res_ex = df[df["gal"] == example_gal]
if not res_ex.empty:
    # Reload data for this galaxy
    data = np.array(galaxy_data[example_gal])
    r_all = data[:, 0]
    v_obs_all = data[:, 1]
    v_gas_all = data[:, 2]
    v_disk_all = data[:, 3]
    v_bulge_all = data[:, 4]
    v_bar_all = np.sqrt(v_gas_all**2 + v_disk_all**2 + v_bulge_all**2)
    mask = (r_all <= R_MAX) & (v_obs_all > 0)
    r_all = r_all[mask]
    v_obs_all = v_obs_all[mask]
    v_bar_all = v_bar_all[mask]
    phi_all = np.log(v_obs_all / v_bar_all)
    ax1.plot(r_all, phi_all, 'k-', lw=1.5, label=r'$\phi_{\rm obs}(r)$')
    # Highlight plateau region
    mask_plt = r_all >= R_START
    ax1.plot(r_all[mask_plt], phi_all[mask_plt], 'bo', markersize=4, alpha=0.7)
    phi_end_ex = res_ex["phi_end"].values[0]
    ax1.axhline(phi_end_ex, color='red', linestyle='--', lw=2,
                label=rf'$\phi_{{\rm end}} = {phi_end_ex:.3f}$')
    ax1.axvline(R_START, color='gray', linestyle=':', lw=1.5,
                label=r'$r_{\rm start}=5$ kpc')
    ax1.set_xlabel(r'$r$ [kpc]')
    ax1.set_ylabel(r'$\phi$')
    ax1.set_title(f'Causal field in {example_gal}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig("fig1_example_phi.pdf", bbox_inches="tight")
    print("Figure 1 saved.")
else:
    print("Figure 1 skipped (NGC3198 not in sample).")

# -----------------------------------------------------------------
# Figure 2: Distribution of phi_end
# -----------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(8, 6))
ax2.hist(phi_end, bins=15, color='steelblue', edgecolor='black', alpha=0.8)
ax2.axvline(PHI_MAX, color='red', linestyle='--', lw=2,
            label=rf'$\phi_{{\max}}={PHI_MAX}$ (upper bound)')
ax2.axvline(np.mean(phi_end), color='orange', linestyle='-.', lw=2,
            label=rf'Mean = {np.mean(phi_end):.3f}')
ax2.set_xlabel(r'$\phi_{\rm end}$')
ax2.set_ylabel('Number of galaxies')
ax2.set_title(f'Distribution of saturation values (N = {n_success})')
ax2.legend()
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig("fig2_phi_end_distribution.pdf", bbox_inches="tight")
print("Figure 2 saved.")

# -----------------------------------------------------------------
# Figure 3: Concentration Law (left) and weak vs Sigma (right)
# -----------------------------------------------------------------
fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(16, 6))

# Left: φ_end vs logM
ax3a.scatter(logM, phi_end, c='steelblue', edgecolors='black',
             alpha=0.8, s=80, zorder=5)
x_fit = np.linspace(logM.min(), logM.max(), 100)
y_fit = intercept + slope * x_fit
ax3a.plot(x_fit, y_fit, 'r-', lw=2.5, label=rf'$R^2 = {r_value**2:.3f}$')
ax3a.set_xlabel(r'$\log_{10}(M_{\rm tot}/M_\odot)$')
ax3a.set_ylabel(r'$\phi_{\rm end}$')
ax3a.set_title('Concentration Law (observed)')
ax3a.legend(loc='upper right')
ax3a.grid(True, alpha=0.3)

# Right: φ_end vs logSigma
ax3b.scatter(logS, phi_end, c='darkgreen', edgecolors='black',
             alpha=0.8, s=80)
s_fit = np.polyfit(logS, phi_end, 1)
s_poly = np.poly1d(s_fit)
x_s = np.linspace(logS.min(), logS.max(), 100)
ax3b.plot(x_s, s_poly(x_s), 'r--', lw=2, label=rf'$r = {r_S:.3f}$')
ax3b.set_xlabel(r'$\log_{10}(\Sigma / M_\odot\,{\rm kpc}^{-2})$')
ax3b.set_ylabel(r'$\phi_{\rm end}$')
ax3b.set_title('Weak correlation with surface density')
ax3b.legend(loc='upper right')
ax3b.grid(True, alpha=0.3)

fig3.tight_layout()
fig3.savefig("fig3_concentration_law.pdf", bbox_inches="tight")
print("Figure 3 saved.")

# -----------------------------------------------------------------
# Figure 4: Residuals of the linear fit
# -----------------------------------------------------------------
phi_pred_linear = intercept + slope * logM
residuals = phi_end - phi_pred_linear

fig4, ax4 = plt.subplots(figsize=(8, 6))
ax4.scatter(logM, residuals, c='purple', edgecolors='black',
            alpha=0.8, s=60)
ax4.axhline(0, color='k', linestyle='-', lw=1)
# Trend line to check randomness
z = np.polyfit(logM, residuals, 1)
p = np.poly1d(z)
x_res = np.linspace(logM.min(), logM.max(), 100)
ax4.plot(x_res, p(x_res), 'r--', lw=1.5, label='Trend (should be flat)')
ax4.set_xlabel(r'$\log_{10}(M_{\rm tot}/M_\odot)$')
ax4.set_ylabel(r'Residuals ($\phi_{\rm end,obs} - \phi_{\rm end,pred}$)')
ax4.set_title('Residuals of the Concentration Law')
ax4.legend()
ax4.grid(True, alpha=0.3)
fig4.tight_layout()
fig4.savefig("fig4_residuals.pdf", bbox_inches="tight")
print("Figure 4 saved.")

# =============================================================================
# 7.  EXPORT RESULTS
# =============================================================================
df_out = df[["gal", "logM", "logSigma", "phi_end", "phi_std", "n_outer"]]
df_out.to_csv("concentration_law_results.csv", index=False)
print("\nResults saved to concentration_law_results.csv")

# =============================================================================
# 8.  SUMMARY
# =============================================================================
print("\n" + "=" * 60)
print("SUMMARY OF OBSERVATIONAL CONCENTRATION LAW")
print("=" * 60)
print(f"Galaxies analysed            : {n_success}")
print(f"Mean φ_end                   : {np.mean(phi_end):.4f}")
print(f"Median φ_end                 : {np.median(phi_end):.4f}")
print(f"Std φ_end                    : {np.std(phi_end):.4f}")
print(f"Min φ_end                    : {np.min(phi_end):.4f}")
print(f"Max φ_end                    : {np.max(phi_end):.4f}")
print(f"Correlation with logM        : r = {r_M:.4f} (p = {p_M:.4e})")
print(f"Correlation with logΣ        : r = {r_S:.4f} (p = {p_S:.4e})")
print(f"Concentration law (linear)   : φ_end = {intercept:.4f} + {slope:.4f}·logM")
print(f"R² (linear)                  : {r_value**2:.4f}")
print(f"Multiple regression R²       : {R2_mult:.4f}")
print("=" * 60)