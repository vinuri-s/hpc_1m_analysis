"""
Figure: Thor-database basalt UCS distribution, with the 7 sampled K/UCS
points marked and the excluded Schmidt-hammer-saturation outliers
flagged, plus the UCS -> K conversion curve. Identical to the 15 m
pipeline's sens_07_ucs_distribution.py -- doesn't depend on grid
resolution, only on the Thor spreadsheet, included here for a
self-contained package.

Usage:
    python3 05_ucs_distribution.py
(expects "Search tool.xlsx" in this same directory -- copy it over from
the main project folder if it isn't here already)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = HERE
XLSX_PATH = os.path.join(HERE, "Search tool.xlsx")

df = pd.read_excel(XLSX_PATH, sheet_name="RN data")
basalt = df[df["Class III - Lito"] == "Basalt"].copy()
ucs = pd.to_numeric(basalt["UCS from RN"], errors="coerce").dropna()

UCS_REF = 37.741401877848745
K_REF = 1.0e-7

K_VALUES = {
    "K1_3.02e-09": 3.024e-09,
    "K2_1.19e-08": 1.191e-08,
    "K3_4.69e-08": 4.693e-08,
    "K4_1.85e-07": 1.849e-07,
    "K5_7.28e-07": 7.283e-07,
    "K6_2.87e-06": 2.869e-06,
    "K7_1.13e-05": 1.130e-05,
}
sample_ucs = {k: UCS_REF * np.sqrt(K_REF / v) for k, v in K_VALUES.items()}

excluded = ucs[ucs > 300]
pct99 = ucs.quantile(0.99)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
bins = np.logspace(np.log10(ucs.min()), np.log10(ucs.max()), 40)
ax.hist(ucs, bins=bins, color="#8c9fc4", edgecolor="white", alpha=0.85,
        label=f"All basalt UCS (n={len(ucs)})")
ax.set_xscale("log")
ax.axvline(pct99, color="k", ls="--", lw=1.3, label=f"99th pct = {pct99:.1f} MPa (used as max)")
ax.axvline(ucs.min(), color="tab:green", ls=":", lw=1.3, label=f"min = {ucs.min():.2f} MPa (used as min)")
for u in excluded:
    ax.axvline(u, color="tab:red", ls="-", lw=1.5, alpha=0.7)
ax.axvline(excluded.iloc[0], color="tab:red", ls="-", lw=1.5, alpha=0.7,
           label="Excluded outliers (Schmidt-hammer\nsaturation artefacts, Teymen &\nMenguec 2020)")
# Okabe-Ito colorblind-safe qualitative palette (Okabe & Ito, 2008), same
# K1-K7 assignment used across all figures in the paper (see
# 04_make_figures_1kyr.py), replacing viridis_r so K colours are
# consistent between the main-text and supplementary figures.
OKABE_ITO_7 = ["#0072B2", "#56B4E9", "#009E73", "#F0E442",
               "#E69F00", "#D55E00", "#CC79A7"]
for i, (k, u) in enumerate(sample_ucs.items()):
    ax.axvline(u, color=OKABE_ITO_7[i], lw=1.8, alpha=0.9)
    ax.text(u, ax.get_ylim()[1] * (0.92 - 0.06 * (i % 3)), f"K{i+1}",
            rotation=90, fontsize=8, color=OKABE_ITO_7[i], ha="right", va="top")
ax.set_xlabel("UCS (MPa, log scale)")
ax.set_ylabel("Count")
ax.set_title("(a) Basalt UCS distribution (Thor database) and the\n7 sampled erodibility values", fontsize=10)
ax.legend(fontsize=7, loc="upper right")

ax = axes[1]
u_curve = np.logspace(np.log10(1), np.log10(1000), 300)
k_curve = K_REF * (UCS_REF / u_curve) ** 2
ax.plot(u_curve, k_curve, color="grey", lw=1.5, label=r"$K=K_{ref}(UCS_{ref}/UCS)^2$")
ax.scatter([UCS_REF], [K_REF], color="k", zorder=5, marker="*", s=120,
           label=f"Anchor: Quye-Sawyer et al. (2020)\n$K_{{ref}}$={K_REF:.1e} yr$^{{-1}}$ @ UCS$_{{ref}}$={UCS_REF:.1f} MPa")
for i, (k, u) in enumerate(sample_ucs.items()):
    ax.scatter([u], [K_VALUES[k]], color=OKABE_ITO_7[i], zorder=5, s=50)
    ax.annotate(f"K{i+1}", (u, K_VALUES[k]), fontsize=8, color=OKABE_ITO_7[i],
                xytext=(5, 5), textcoords="offset points")
ax.scatter(excluded.values, [K_REF * (UCS_REF / u) ** 2 for u in excluded.values],
           color="tab:red", marker="x", s=70, zorder=5, label="Excluded outliers (would give\nK as low as 1.5e-10)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("UCS (MPa, log scale)")
ax.set_ylabel(r"$K_{br}$ (yr$^{-1}$, log scale)")
ax.set_title("(b) UCS -> erodibility conversion and anchor point", fontsize=10)
ax.legend(fontsize=7, loc="upper right")

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/fig_ucs_distribution_1m.png", dpi=230, bbox_inches="tight")
plt.savefig(f"{OUT_DIR}/fig_ucs_distribution_1m.pdf", bbox_inches="tight")
print("Saved fig_ucs_distribution_1m.png/pdf")
