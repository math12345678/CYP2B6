"""
KDE density plots for RMSD, per Prof. Bishop's request (July 29 meeting):
WT should show a sharp, tight peak with good replicate agreement. Mutants that
sample multiple conformations will show broader or multimodal density.

Run from ~/Desktop/Research/Research_Projects/CYP2B6 after the RMSD .xvg files exist.
Produces one grid figure (kde_rmsd_all_alleles.png) with all 11 alleles vs WT.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

base = "."

# (allele label, rep1 file, rep2 file)
alleles = [
    ("G99E",        "G99E/rmsd_G99E.xvg",               "G99E_2/rmsd_G99E_2.xvg"),
    ("K139E",       "K139E/rmsd_K139E.xvg",             "K139E_2/rmsd_K139E_2.xvg"),
    ("M46V",        "M46V/rmsd_M46V.xvg",               "M46V_2/rmsd_M46V_2.xvg"),
    ("I328T",       "I328T/rmsd_I328T.xvg",             "I328T_2/rmsd_I328T_2.xvg"),
    ("I391N",       "I391N/rmsd_I391N.xvg",             "I391N_2/rmsd_I391N_2.xvg"),
    ("K262R",       "K262R/rmsd_K262R.xvg",             "K262R_2/rmsd_K262R_2.xvg"),
    ("R140Q",       "R140Q/rmsd_R140Q.xvg",             "R140Q_2/rmsd_R140Q_2.xvg"),
    ("R487C",       "R487C/rmsd_R487C.xvg",             "R487C_2/rmsd_R487C_2.xvg"),
    ("P428T",       "P428T/rmsd_P428T.xvg",             "P428T_2/rmsd_P428T_2.xvg"),
    ("S259R",       "S259R/rmsd_S259R.xvg",             "S259R_2/rmsd_S259R_2.xvg"),
    ("T306S-R378K", "T306S-R378K/rmsd_T306S-R378K.xvg", "T306S-R378K_2/rmsd_T306S-R378K_2.xvg"),
]

wt1 = np.loadtxt(f"{base}/WT/rmsd_WT.xvg", comments=["#", "@"])[:, 1]
wt2 = np.loadtxt(f"{base}/WT_2/rmsd_WT.xvg", comments=["#", "@"])[:, 1]

def kde_curve(data, grid):
    kde = gaussian_kde(data)
    return kde(grid)

fig, axes = plt.subplots(4, 3, figsize=(16, 18))
axes = axes.flatten()

xmin = min(wt1.min(), wt2.min())
xmax = max(wt1.max(), wt2.max())

for i, (name, f1, f2) in enumerate(alleles):
    ax = axes[i]
    try:
        m1 = np.loadtxt(f"{base}/{f1}", comments=["#", "@"])[:, 1]
        m2 = np.loadtxt(f"{base}/{f2}", comments=["#", "@"])[:, 1]
    except OSError:
        ax.set_visible(False)
        continue

    lo = min(xmin, m1.min(), m2.min())
    hi = max(xmax, m1.max(), m2.max())
    grid = np.linspace(lo, hi, 400)

    ax.plot(grid, kde_curve(wt1, grid), color="black", label="WT rep1", linewidth=1.5)
    ax.plot(grid, kde_curve(wt2, grid), color="gray", label="WT rep2", linewidth=1.5)
    ax.plot(grid, kde_curve(m1, grid), color="tab:blue", label=f"{name} rep1", linewidth=1.5)
    ax.plot(grid, kde_curve(m2, grid), color="tab:orange", label=f"{name} rep2", linewidth=1.5)

    ax.set_title(f"WT vs. {name}", fontsize=14, fontweight="bold")
    ax.set_xlabel("RMSD (nm)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.tick_params(labelsize=10)
    ax.legend(fontsize=9, loc="upper right")

# hide any unused axes (12 slots, 10 alleles)
for j in range(len(alleles), len(axes)):
    axes[j].set_visible(False)

fig.suptitle("Backbone RMSD — KDE density, all alleles vs. WT (both replicates each)",
             fontsize=15, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("kde_rmsd_all_alleles.png", dpi=150)
plt.show()

print("Saved kde_rmsd_all_alleles.png")
