"""
Quantitative significance check for RMSD and RMSF findings, using the
disagreement between the two WT replicates as an empirical noise floor
(since only 2 WT replicates exist here, not the reference triplicate used
in the CYP3A4 paper's 3-SD framework -- this is the best available proxy).

A mutant's effect is only treated as "real" if:
  (a) both mutant replicates agree in the direction of the delta vs. WT, AND
  (b) the replicate-averaged delta exceeds the WT rep1-vs-rep2 noise floor.

For RMSF mutation-site comparisons, a +/-3 residue window max is used
instead of the exact single-residue value, since RMSF peaks can shift by a
residue or two between systems and a single-point comparison is not robust
to that.

Run from ~/Desktop/Research/Research_Projects/CYP2B6.
"""
import numpy as np

def load(path):
    return np.loadtxt(path, comments=["#", "@"])

def window_max(residues, arr, target, w=3):
    mask = (residues >= target - w) & (residues <= target + w)
    return arr[mask].max()

# ===== RMSD significance =====
wt1_rmsd = load("WT/rmsd_WT.xvg")[:, 1]
wt2_rmsd = load("WT_2/rmsd_WT.xvg")[:, 1]
rmsd_noise_floor = abs(wt1_rmsd.mean() - wt2_rmsd.mean())

alleles_rmsd = {
    "G99E": ("G99E/rmsd_G99E.xvg", "G99E_2/rmsd_G99E_2.xvg"),
    "K139E": ("K139E/rmsd_K139E.xvg", "K139E_2/rmsd_K139E_2.xvg"),
    "M46V": ("M46V/rmsd_M46V.xvg", "M46V_2/rmsd_M46V_2.xvg"),
    "I328T": ("I328T/rmsd_I328T.xvg", "I328T_2/rmsd_I328T_2.xvg"),
    "I391N": ("I391N/rmsd_I391N.xvg", "I391N_2/rmsd_I391N_2.xvg"),
    "K262R": ("K262R/rmsd_K262R.xvg", "K262R_2/rmsd_K262R_2.xvg"),
    "R140Q": ("R140Q/rmsd_R140Q.xvg", "R140Q_2/rmsd_R140Q_2.xvg"),
    "R487C": ("R487C/rmsd_R487C.xvg", "R487C_2/rmsd_R487C_2.xvg"),
    "P428T": ("P428T/rmsd_P428T.xvg", "P428T_2/rmsd_P428T_2.xvg"),
    "S259R": ("S259R/rmsd_S259R.xvg", "S259R_2/rmsd_S259R_2.xvg"),
    "T306S-R378K": ("T306S-R378K/rmsd_T306S-R378K.xvg", "T306S-R378K_2/rmsd_T306S-R378K_2.xvg"),
}

print(f"RMSD noise floor (|WT rep1 mean - WT rep2 mean|): {rmsd_noise_floor:.4f} nm\n")
print(f"{'Allele':15s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
rmsd_results = {}
for name, (f1, f2) in alleles_rmsd.items():
    m1 = load(f1)[:, 1]
    m2 = load(f2)[:, 1]
    d1 = m1.mean() - wt1_rmsd.mean()
    d2 = m2.mean() - wt2_rmsd.mean()
    d_avg = ((m1.mean() + m2.mean()) / 2) - ((wt1_rmsd.mean() + wt2_rmsd.mean()) / 2)
    agree = (d1 > 0) == (d2 > 0)
    robust = agree and abs(d_avg) > rmsd_noise_floor
    rmsd_results[name] = robust
    print(f"{name:15s} {d_avg:8.4f} {d1:8.4f} {d2:8.4f} {str(agree):>7s} {str(robust):>8s}")

# ===== RMSF mutation-site significance =====
wt1_rmsf = load("WT/rmsf_WT.xvg")
wt2_rmsf = load("WT_2/rmsf_WT_2.xvg")
residues = wt1_rmsf[:, 0].astype(int)

sites = {
    "G99E": ("G99E/rmsf_G99E.xvg", "G99E_2/rmsf_G99E_2.xvg", 71),
    "K139E": ("K139E/rmsf_K139E.xvg", "K139E_2/rmsf_K139E_2.xvg", 111),
    "M46V": ("M46V/rmsf_M46V.xvg", "M46V_2/rmsf_M46V_2.xvg", 18),
    "I328T": ("I328T/rmsf_I328T.xvg", "I328T_2/rmsf_I328T_2.xvg", 300),
    "I391N": ("I391N/rmsf_I391N.xvg", "I391N_2/rmsf_I391N_2.xvg", 363),
    "K262R": ("K262R/rmsf_K262R.xvg", "K262R_2/rmsf_K262R_2.xvg", 234),
    "R140Q": ("R140Q/rmsf_R140Q.xvg", "R140Q_2/rmsf_R140Q_2.xvg", 112),
    "R487C": ("R487C/rmsf_R487C.xvg", "R487C_2/rmsf_R487C_2.xvg", 459),
    "P428T": ("P428T/rmsf_P428T.xvg", "P428T_2/rmsf_P428T_2.xvg", 400),
    "S259R": ("S259R/rmsf_S259R.xvg", "S259R_2/rmsf_S259R_2.xvg", 231),
}

print("\nRMSF mutation-site significance (+/-3 residue window max)")
print(f"{'Allele':10s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (f1, f2, site) in sites.items():
    m1 = load(f1); m2 = load(f2)
    wt_v1 = window_max(residues, wt1_rmsf[:, 1], site)
    wt_v2 = window_max(residues, wt2_rmsf[:, 1], site)
    m_v1 = window_max(residues, m1[:, 1], site)
    m_v2 = window_max(residues, m2[:, 1], site)
    noise_floor = abs(wt_v1 - wt_v2)
    d1 = m_v1 - wt_v1
    d2 = m_v2 - wt_v2
    d_avg = ((m_v1 + m_v2) / 2) - ((wt_v1 + wt_v2) / 2)
    agree = (d1 > 0) == (d2 > 0)
    robust = agree and abs(d_avg) > noise_floor
    print(f"{name:10s} {d_avg:8.4f} {d1:8.4f} {d2:8.4f} {str(agree):>7s} {str(robust):>8s}")

# T306S-R378K, two sites
t1 = load("T306S-R378K/rmsf_T306S-R378K.xvg")
t2 = load("T306S-R378K_2/rmsf_T306S-R378K_2.xvg")
for site, label in [(278, "T306"), (350, "R378")]:
    wt_v1 = window_max(residues, wt1_rmsf[:, 1], site)
    wt_v2 = window_max(residues, wt2_rmsf[:, 1], site)
    m_v1 = window_max(residues, t1[:, 1], site)
    m_v2 = window_max(residues, t2[:, 1], site)
    noise_floor = abs(wt_v1 - wt_v2)
    d1 = m_v1 - wt_v1
    d2 = m_v2 - wt_v2
    d_avg = ((m_v1 + m_v2) / 2) - ((wt_v1 + wt_v2) / 2)
    agree = (d1 > 0) == (d2 > 0)
    robust = agree and abs(d_avg) > noise_floor
    print(f"T306S-R378K({label}):{'':2s} {d_avg:8.4f} {d1:8.4f} {d2:8.4f} {str(agree):>7s} {str(robust):>8s}")
