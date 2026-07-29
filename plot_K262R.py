import numpy as np
import matplotlib.pyplot as plt

base = "."  # run this from ~/Desktop/Research/Research_Projects/CYP2B6

wt1_rmsd = np.loadtxt(f"{base}/WT/rmsd_WT.xvg", comments=["#", "@"])
wt2_rmsd = np.loadtxt(f"{base}/WT_2/rmsd_WT.xvg", comments=["#", "@"])
k1_rmsd = np.loadtxt(f"{base}/K262R/rmsd_K262R.xvg", comments=["#", "@"])
k2_rmsd = np.loadtxt(f"{base}/K262R_2/rmsd_K262R_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsd[:, 0], wt1_rmsd[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsd[:, 0], wt2_rmsd[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(k1_rmsd[:, 0], k1_rmsd[:, 1], label="K262R rep1", color="tab:blue", alpha=0.7)
plt.plot(k2_rmsd[:, 0], k2_rmsd[:, 1], label="K262R rep2", color="tab:orange", alpha=0.7)
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("WT vs. K262R — Backbone RMSD")
plt.legend()
plt.tight_layout()
plt.savefig("rmsd_WT_vs_K262R.png", dpi=150)
plt.show()

wt1_rmsf = np.loadtxt(f"{base}/WT/rmsf_WT.xvg", comments=["#", "@"])
wt2_rmsf = np.loadtxt(f"{base}/WT_2/rmsf_WT_2.xvg", comments=["#", "@"])
k1_rmsf = np.loadtxt(f"{base}/K262R/rmsf_K262R.xvg", comments=["#", "@"])
k2_rmsf = np.loadtxt(f"{base}/K262R_2/rmsf_K262R_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsf[:, 0], wt1_rmsf[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsf[:, 0], wt2_rmsf[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(k1_rmsf[:, 0], k1_rmsf[:, 1], label="K262R rep1", color="tab:blue", alpha=0.7)
plt.plot(k2_rmsf[:, 0], k2_rmsf[:, 1], label="K262R rep2", color="tab:orange", alpha=0.7)
plt.axvline(234, color="red", linestyle="--", linewidth=1,
            label="Residue 234 (GROMACS numbering) = true K262 mutation site")
plt.axvline(108, color="green", linestyle=":", linewidth=1, alpha=0.7,
            label="Residues 108-112 (GROMACS) = shared hotspot loop")
plt.axvline(112, color="green", linestyle=":", linewidth=1, alpha=0.7)
plt.xlabel("Residue number (GROMACS numbering)")
plt.ylabel("RMSF (nm)")
plt.title("WT vs. K262R — per-residue Backbone RMSF")
plt.legend()
plt.tight_layout()
plt.savefig("rmsf_WT_vs_K262R.png", dpi=150)
plt.show()
