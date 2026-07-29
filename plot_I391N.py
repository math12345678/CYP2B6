import numpy as np
import matplotlib.pyplot as plt

base = "."  # run this from ~/Desktop/Research/Research_Projects/CYP2B6

wt1_rmsd = np.loadtxt(f"{base}/WT/rmsd_WT.xvg", comments=["#", "@"])
wt2_rmsd = np.loadtxt(f"{base}/WT_2/rmsd_WT.xvg", comments=["#", "@"])
i1_rmsd = np.loadtxt(f"{base}/I391N/rmsd_I391N.xvg", comments=["#", "@"])
i2_rmsd = np.loadtxt(f"{base}/I391N_2/rmsd_I391N_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsd[:, 0], wt1_rmsd[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsd[:, 0], wt2_rmsd[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(i1_rmsd[:, 0], i1_rmsd[:, 1], label="I391N rep1", color="tab:blue", alpha=0.7)
plt.plot(i2_rmsd[:, 0], i2_rmsd[:, 1], label="I391N rep2", color="tab:orange", alpha=0.7)
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("WT vs. I391N — Backbone RMSD")
plt.legend()
plt.tight_layout()
plt.savefig("rmsd_WT_vs_I391N.png", dpi=150)
plt.show()

wt1_rmsf = np.loadtxt(f"{base}/WT/rmsf_WT.xvg", comments=["#", "@"])
wt2_rmsf = np.loadtxt(f"{base}/WT_2/rmsf_WT_2.xvg", comments=["#", "@"])
i1_rmsf = np.loadtxt(f"{base}/I391N/rmsf_I391N.xvg", comments=["#", "@"])
i2_rmsf = np.loadtxt(f"{base}/I391N_2/rmsf_I391N_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsf[:, 0], wt1_rmsf[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsf[:, 0], wt2_rmsf[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(i1_rmsf[:, 0], i1_rmsf[:, 1], label="I391N rep1", color="tab:blue", alpha=0.7)
plt.plot(i2_rmsf[:, 0], i2_rmsf[:, 1], label="I391N rep2", color="tab:orange", alpha=0.7)
plt.axvline(363, color="red", linestyle="--", linewidth=1,
            label="Residue 363 (GROMACS numbering) = true I391 mutation site")
plt.axvline(108, color="green", linestyle=":", linewidth=1, alpha=0.7,
            label="Residues 108-112 (GROMACS) = shared hotspot loop")
plt.axvline(112, color="green", linestyle=":", linewidth=1, alpha=0.7)
plt.xlabel("Residue number (GROMACS numbering)")
plt.ylabel("RMSF (nm)")
plt.title("WT vs. I391N — per-residue Backbone RMSF")
plt.legend()
plt.tight_layout()
plt.savefig("rmsf_WT_vs_I391N.png", dpi=150)
plt.show()
