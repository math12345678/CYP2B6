"""
DRN centrality heatmap across all 11 alleles, mirroring plot_rmsf_heatmap.py.
Three panels (BC, CC, EC). Rows = allele (averaged across both replicates),
columns = residue number (GROMACS numbering). Color = delta centrality vs.
WT average (mutant - WT_avg). Residue 408 (heme CM1, confirmed via .cif
auth_seq_id cross-reference in drn_significance_check.py) is excluded from
every row, same as the significance check.

Mutation-site markers are drawn per-row at each allele's own site, since
(unlike RMSF/H-bonds) DRN findings here are overwhelmingly local rather than
panel-wide, so a single shared vertical line isn't informative here.

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6, inside the cyp2b6 env
(needs numpy/pandas/matplotlib -- same env used for the other plot_*.py
scripts, not mdmtaskweb).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESID_EXCLUDE = {408}
METRICS = ["BC", "CC", "EC"]

alleles = [
    ("G99E", "G99E", "G99E_2", 71),
    ("K139E", "K139E", "K139E_2", 111),
    ("M46V", "M46V", "M46V_2", 18),
    ("I328T", "I328T", "I328T_2", 300),
    ("I391N", "I391N", "I391N_2", 363),
    ("K262R", "K262R", "K262R_2", 234),
    ("R140Q", "R140Q", "R140Q_2", 112),
    ("R487C", "R487C", "R487C_2", 459),
    ("P428T", "P428T", "P428T_2", 400),
    ("S259R", "S259R", "S259R_2", 231),
    ("T306S-R378K", "T306S-R378K", "T306S-R378K_2", None),  # two sites, handled separately
]

def load_metric(path):
    df = pd.read_csv(path)
    df.index = range(1, len(df) + 1)
    df = df.drop(index=[r for r in RESID_EXCLUDE if r in df.index])
    return df

wt1 = load_metric("WT/md_noWAT_mean.csv")
wt2 = load_metric("WT_2/md_noWAT_mean.csv")
all_resids = sorted(wt1.index)  # 1..463 minus 408, shared across all systems

fig, axes = plt.subplots(3, 1, figsize=(22, 16))

for ax, metric in zip(axes, METRICS):
    wt_avg = (wt1[metric].reindex(all_resids) + wt2[metric].reindex(all_resids)) / 2

    rows, labels, sites = [], [], []
    for name, d1, d2, site in alleles:
        m1 = load_metric(f"{d1}/md_noWAT_mean.csv")[metric].reindex(all_resids)
        m2 = load_metric(f"{d2}/md_noWAT_mean.csv")[metric].reindex(all_resids)
        m_avg = (m1 + m2) / 2
        delta = (m_avg - wt_avg).values
        rows.append(delta)
        labels.append(name)
        sites.append(site)

    matrix = np.array(rows)
    vmax = np.nanmax(np.abs(matrix))
    im = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   extent=[all_resids[0], all_resids[-1], len(labels), 0])

    ax.set_yticks(np.arange(len(labels)) + 0.5)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Residue number (GROMACS numbering)", fontsize=10)
    ax.set_title(f"{metric} delta vs. WT (mutant avg - WT avg)", fontsize=12, fontweight="bold")

    # mark each allele's own mutation site with a short tick on its row
    for i, site in enumerate(sites):
        if site is not None:
            ax.plot(site, i + 0.5, marker="|", color="black", markersize=14, markeredgewidth=1.5)
        else:
            # T306S-R378K: mark both sites
            ax.plot(278, i + 0.5, marker="|", color="black", markersize=14, markeredgewidth=1.5)
            ax.plot(350, i + 0.5, marker="|", color="black", markersize=14, markeredgewidth=1.5)

    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label(f"Delta {metric}", fontsize=9)

fig.suptitle("DRN centrality delta vs. WT, all alleles (black ticks = each allele's own mutation site)\n"
             "Residue 408 (heme CM1) excluded from all rows", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("drn_heatmap_all_alleles.png", dpi=150)
print("Saved drn_heatmap_all_alleles.png")
