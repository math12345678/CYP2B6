import numpy as np
import matplotlib.pyplot as plt

base = "."  # run this from ~/Desktop/Research/Research_Projects/CYP2B6

wt1_rmsd = np.loadtxt(f"{base}/WT/rmsd_WT.xvg", comments=["#", "@"])
wt2_rmsd = np.loadtxt(f"{base}/WT_2/rmsd_WT.xvg", comments=["#", "@"])
m1_rmsd = np.loadtxt(f"{base}/M46V/rmsd_M46V.xvg", comments=["#", "@"])
m2_rmsd = np.loadtxt(f"{base}/M46V_2/rmsd_M46V_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsd[:, 0], wt1_rmsd[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsd[:, 0], wt2_rmsd[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(m1_rmsd[:, 0], m1_rmsd[:, 1], label="M46V rep1", color="tab:blue", alpha=0.7)
plt.plot(m2_rmsd[:, 0], m2_rmsd[:, 1], label="M46V rep2", color="tab:orange", alpha=0.7)
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("WT vs. M46V — Backbone RMSD")
plt.legend()
plt.tight_layout()
plt.savefig("rmsd_WT_vs_M46V.png", dpi=150)
plt.show()

wt1_rmsf = np.loadtxt(f"{base}/WT/rmsf_WT.xvg", comments=["#", "@"])
wt2_rmsf = np.loadtxt(f"{base}/WT_2/rmsf_WT_2.xvg", comments=["#", "@"])
m1_rmsf = np.loadtxt(f"{base}/M46V/rmsf_M46V.xvg", comments=["#", "@"])
m2_rmsf = np.loadtxt(f"{base}/M46V_2/rmsf_M46V_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsf[:, 0], wt1_rmsf[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsf[:, 0], wt2_rmsf[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(m1_rmsf[:, 0], m1_rmsf[:, 1], label="M46V rep1", color="tab:blue", alpha=0.7)
plt.plot(m2_rmsf[:, 0], m2_rmsf[:, 1], label="M46V rep2", color="tab:orange", alpha=0.7)
plt.axvline(18, color="red", linestyle="--", linewidth=1,
            label="Residue 18 (GROMACS numbering) = true M46 mutation site")
plt.axvline(108, color="green", linestyle=":", linewidth=1, alpha=0.7,
            label="Residues 108-112 (GROMACS) = shared G99E/K139E hotspot loop")
plt.axvline(112, color="green", linestyle=":", linewidth=1, alpha=0.7)
plt.xlabel("Residue number (GROMACS numbering)")
plt.ylabel("RMSF (nm)")
plt.title("WT vs. M46V — per-residue Backbone RMSF")
plt.legend()
plt.tight_layout()
plt.savefig("rmsf_WT_vs_M46V.png", dpi=150)
plt.show()
