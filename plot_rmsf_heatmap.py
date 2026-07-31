"""
RMSF heatmap across all 11 alleles, per Prof. Bishop's request (July 29 meeting).
Rows = allele (averaged across both replicates), columns = residue number
(GROMACS numbering). Color = delta RMSF vs. WT average (mutant - WT_avg),
following the same delta convention as the CYP3A4 paper (red = more flexible
than WT, blue = more rigid than WT).

Run from ~/Desktop/Research/Research_Projects/CYP2B6.
"""
import numpy as np
import matplotlib.pyplot as plt

base = "."

alleles = [
    ("G99E",        "G99E/rmsf_G99E.xvg",               "G99E_2/rmsf_G99E_2.xvg"),
    ("K139E",       "K139E/rmsf_K139E.xvg",             "K139E_2/rmsf_K139E_2.xvg"),
    ("M46V",        "M46V/rmsf_M46V.xvg",               "M46V_2/rmsf_M46V_2.xvg"),
    ("I328T",       "I328T/rmsf_I328T.xvg",             "I328T_2/rmsf_I328T_2.xvg"),
    ("I391N",       "I391N/rmsf_I391N.xvg",             "I391N_2/rmsf_I391N_2.xvg"),
    ("K262R",       "K262R/rmsf_K262R.xvg",             "K262R_2/rmsf_K262R_2.xvg"),
    ("R140Q",       "R140Q/rmsf_R140Q.xvg",             "R140Q_2/rmsf_R140Q_2.xvg"),
    ("R487C",       "R487C/rmsf_R487C.xvg",             "R487C_2/rmsf_R487C_2.xvg"),
    ("P428T",       "P428T/rmsf_P428T.xvg",             "P428T_2/rmsf_P428T_2.xvg"),
    ("S259R",       "S259R/rmsf_S259R.xvg",             "S259R_2/rmsf_S259R_2.xvg"),
    ("T306S-R378K", "T306S-R378K/rmsf_T306S-R378K.xvg", "T306S-R378K_2/rmsf_T306S-R378K_2.xvg"),
]

wt1 = np.loadtxt(f"{base}/WT/rmsf_WT.xvg", comments=["#", "@"])
wt2 = np.loadtxt(f"{base}/WT_2/rmsf_WT_2.xvg", comments=["#", "@"])
residues = wt1[:, 0]
wt_avg = (wt1[:, 1] + wt2[:, 1]) / 2

rows = []
labels = []
for name, f1, f2 in alleles:
    m1 = np.loadtxt(f"{base}/{f1}", comments=["#", "@"])
    m2 = np.loadtxt(f"{base}/{f2}", comments=["#", "@"])
    m_avg = (m1[:, 1] + m2[:, 1]) / 2
    delta = m_avg - wt_avg
    rows.append(delta)
    labels.append(name)

matrix = np.array(rows)

fig, ax = plt.subplots(figsize=(20, 6))
vmax = np.max(np.abs(matrix))
im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
               extent=[residues[0], residues[-1], len(labels), 0])

ax.set_yticks(np.arange(len(labels)) + 0.5)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("Residue number (GROMACS numbering)", fontsize=11)
ax.set_title("RMSF delta vs. WT (mutant avg - WT avg), all alleles\n"
             "Red = more flexible than WT, Blue = more rigid than WT",
             fontsize=13, fontweight="bold")

# mark the shared hotspot loop and known mutation sites for reference
ax.axvline(108, color="black", linestyle=":", linewidth=1, alpha=0.6)
ax.axvline(112, color="black", linestyle=":", linewidth=1, alpha=0.6)

cbar = fig.colorbar(im, ax=ax, pad=0.01)
cbar.set_label("Delta RMSF (nm)", fontsize=10)

plt.tight_layout()
plt.savefig("rmsf_heatmap_all_alleles.png", dpi=150)
plt.show()

print("Saved rmsf_heatmap_all_alleles.png")
