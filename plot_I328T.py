import numpy as np
import matplotlib.pyplot as plt

base = "."  # run this from ~/Desktop/Research/Research_Projects/CYP2B6

wt1_rmsd = np.loadtxt(f"{base}/WT/rmsd_WT.xvg", comments=["#", "@"])
wt2_rmsd = np.loadtxt(f"{base}/WT_2/rmsd_WT.xvg", comments=["#", "@"])
i1_rmsd = np.loadtxt(f"{base}/I328T/rmsd_I328T.xvg", comments=["#", "@"])
i2_rmsd = np.loadtxt(f"{base}/I328T_2/rmsd_I328T_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsd[:, 0], wt1_rmsd[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsd[:, 0], wt2_rmsd[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(i1_rmsd[:, 0], i1_rmsd[:, 1], label="I328T rep1", color="tab:blue", alpha=0.7)
plt.plot(i2_rmsd[:, 0], i2_rmsd[:, 1], label="I328T rep2", color="tab:orange", alpha=0.7)
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("WT vs. I328T — Backbone RMSD")
plt.legend()
plt.tight_layout()
plt.savefig("rmsd_WT_vs_I328T.png", dpi=150)
plt.show()

wt1_rmsf = np.loadtxt(f"{base}/WT/rmsf_WT.xvg", comments=["#", "@"])
wt2_rmsf = np.loadtxt(f"{base}/WT_2/rmsf_WT_2.xvg", comments=["#", "@"])
i1_rmsf = np.loadtxt(f"{base}/I328T/rmsf_I328T.xvg", comments=["#", "@"])
i2_rmsf = np.loadtxt(f"{base}/I328T_2/rmsf_I328T_2.xvg", comments=["#", "@"])

plt.figure(figsize=(10, 5))
plt.plot(wt1_rmsf[:, 0], wt1_rmsf[:, 1], label="WT rep1", color="black", alpha=0.6)
plt.plot(wt2_rmsf[:, 0], wt2_rmsf[:, 1], label="WT rep2", color="gray", alpha=0.6)
plt.plot(i1_rmsf[:, 0], i1_rmsf[:, 1], label="I328T rep1", color="tab:blue", alpha=0.7)
plt.plot(i2_rmsf[:, 0], i2_rmsf[:, 1], label="I328T rep2", color="tab:orange", alpha=0.7)
plt.axvline(300, color="red", linestyle="--", linewidth=1,
            label="Residue 300 (GROMACS numbering) = true I328 mutation site")
plt.axvline(108, color="green", linestyle=":", linewidth=1, alpha=0.7,
            label="Residues 108-112 (GROMACS) = shared G99E/K139E hotspot loop")
plt.axvline(112, color="green", linestyle=":", linewidth=1, alpha=0.7)
plt.xlabel("Residue number (GROMACS numbering)")
plt.ylabel("RMSF (nm)")
plt.title("WT vs. I328T — per-residue Backbone RMSF")
plt.legend()
plt.tight_layout()
plt.savefig("rmsf_WT_vs_I328T.png", dpi=150)
plt.show()
