import numpy as np
import matplotlib.pyplot as plt

base = "."  # run this from ~/Desktop/Research/Research_Projects/CYP2B6

wt1_rmsd = np.loadtxt(f"{base}/WT/rmsd_WT.xvg", comments=["#", "@"])
wt2_rmsd = np.loadtxt(f"{base}/WT_2/rmsd_WT.xvg", comments=["#", "@"])
g1_rmsd = np.loadtxt(f"{base}/G99E/rmsd_G99E.xvg", comments=["#", "@"])
g2_rmsd = np.loadtxt(f"{base}/G99E_2/rmsd_G99E_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsd[:, 0], wt1_rmsd[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsd[:, 0], wt2_rmsd[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(g1_rmsd[:, 0], g1_rmsd[:, 1], label="G99E rep1", color="tab:blue", alpha=0.7)
plt.plot(g2_rmsd[:, 0], g2_rmsd[:, 1], label="G99E rep2", color="tab:orange", alpha=0.7)
plt.xlabel("Time (ns)", fontsize=12)
plt.ylabel("RMSD (nm)", fontsize=12)
plt.title("WT vs. G99E — Backbone RMSD", fontsize=14, fontweight="bold")
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("rmsd_WT_vs_G99E.png", dpi=150)
plt.show()

wt1_rmsf = np.loadtxt(f"{base}/WT/rmsf_WT.xvg", comments=["#", "@"])
wt2_rmsf = np.loadtxt(f"{base}/WT_2/rmsf_WT_2.xvg", comments=["#", "@"])
g1_rmsf = np.loadtxt(f"{base}/G99E/rmsf_G99E.xvg", comments=["#", "@"])
g2_rmsf = np.loadtxt(f"{base}/G99E_2/rmsf_G99E_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsf[:, 0], wt1_rmsf[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsf[:, 0], wt2_rmsf[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(g1_rmsf[:, 0], g1_rmsf[:, 1], label="G99E rep1", color="tab:blue", alpha=0.7)
plt.plot(g2_rmsf[:, 0], g2_rmsf[:, 1], label="G99E rep2", color="tab:orange", alpha=0.7)
plt.axvline(71, color="red", linestyle="--", linewidth=1,
            label="Residue 71 (GROMACS numbering) = true G99 mutation site")
plt.axvline(108, color="green", linestyle=":", linewidth=1, alpha=0.7,
            label="Residues 108-112 (GROMACS) = shared hotspot loop")
plt.axvline(112, color="green", linestyle=":", linewidth=1, alpha=0.7)
plt.xlabel("Residue number (GROMACS numbering)", fontsize=12)
plt.ylabel("RMSF (nm)", fontsize=12)
plt.title("WT vs. G99E — per-residue Backbone RMSF", fontsize=14, fontweight="bold")
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("rmsf_WT_vs_G99E.png", dpi=150)
plt.show()
