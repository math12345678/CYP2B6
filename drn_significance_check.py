"""
Quantitative robustness check for DRN (dynamic residue network) centrality
metrics -- betweenness (BC), closeness (CC), eigenvector centrality (EC) --
mirroring the same framework used for RMSD/RMSF (significance_check.py) and
H-bonds (hbond_significance_check.py): the disagreement between the two WT
replicates is used as an empirical noise floor, and a mutant effect only
counts as "robust" if both of its own replicates agree in direction and the
replicate-averaged delta exceeds that noise floor.

Row numbering: confirmed directly against md_noWAT_mean_BC.cif's auth_seq_id
column (which carries the real residue number) -- CSV row N (1-indexed)
corresponds exactly to GROMACS resid N, the same numbering used throughout
this project (true residue = GROMACS residue + 28). Resid 408 is NOT a real
amino acid: it is the heme cofactor's CM1 component, which survived the
CA/CB-based network reduction because this force field happens to name one
of its carbons "CB". It is excluded explicitly below (RESID_EXCLUDE) rather
than relying on no mutation-site window happening to reach it, following the
same precaution added to hbond_significance_check.py after a near-miss there.

Two metrics per centrality type, analogous to RMSD/RMSF and the H-bond check:
  - Global: mean value of the metric across all (real) residues.
  - Local: window-max value (+/-3 residues) at each allele's own mutation
    site, analogous to the windowed RMSF/H-bond checks.

Run from ~/Desktop/Research/Research_Projects/CYP2B6 after run_all_drn.py.
"""
import numpy as np
import pandas as pd

RESID_EXCLUDE = {408}  # heme CM1, confirmed via .cif auth_seq_id cross-check

METRICS = ["BC", "CC", "EC"]

def load_metric(path):
    df = pd.read_csv(path)
    df.index = range(1, len(df) + 1)  # row N -> GROMACS resid N
    df = df.drop(index=[r for r in RESID_EXCLUDE if r in df.index])
    return df

def window_max(series, target, w=3):
    mask = (series.index >= target - w) & (series.index <= target + w)
    return series[mask].max()

alleles = {
    "G99E": ("G99E", "G99E_2", 71),
    "K139E": ("K139E", "K139E_2", 111),
    "M46V": ("M46V", "M46V_2", 18),
    "I328T": ("I328T", "I328T_2", 300),
    "I391N": ("I391N", "I391N_2", 363),
    "K262R": ("K262R", "K262R_2", 234),
    "R140Q": ("R140Q", "R140Q_2", 112),
    "R487C": ("R487C", "R487C_2", 459),
    "P428T": ("P428T", "P428T_2", 400),
    "S259R": ("S259R", "S259R_2", 231),
}

wt1 = load_metric("WT/md_noWAT_mean.csv")
wt2 = load_metric("WT_2/md_noWAT_mean.csv")

for metric in METRICS:
    print(f"\n{'='*70}\nMetric: {metric}\n{'='*70}")

    wt1_m, wt2_m = wt1[metric], wt2[metric]
    global_noise_floor = abs(wt1_m.mean() - wt2_m.mean())
    print(f"Global noise floor (|WT rep1 mean - WT rep2 mean|): {global_noise_floor:.6f}\n")

    print(f"{'Allele':10s} {'d_avg':>12s} {'rep1 d':>12s} {'rep2 d':>12s} {'agree?':>7s} {'ROBUST?':>8s}")
    for name, (d1, d2, _) in alleles.items():
        m1 = load_metric(f"{d1}/md_noWAT_mean.csv")[metric]
        m2 = load_metric(f"{d2}/md_noWAT_mean.csv")[metric]
        delta1 = m1.mean() - wt1_m.mean()
        delta2 = m2.mean() - wt2_m.mean()
        d_avg = ((m1.mean() + m2.mean()) / 2) - ((wt1_m.mean() + wt2_m.mean()) / 2)
        agree = (delta1 > 0) == (delta2 > 0)
        robust = agree and abs(d_avg) > global_noise_floor
        print(f"{name:10s} {d_avg:12.6f} {delta1:12.6f} {delta2:12.6f} {str(agree):>7s} {str(robust):>8s}")

    print(f"\nLocal {metric} significance at each allele's own mutation site (+/-3 residue window max)")
    print(f"{'Allele':10s} {'d_avg':>12s} {'rep1 d':>12s} {'rep2 d':>12s} {'agree?':>7s} {'ROBUST?':>8s}")
    for name, (d1, d2, site) in alleles.items():
        wt_v1 = window_max(wt1_m, site)
        wt_v2 = window_max(wt2_m, site)
        noise_floor = abs(wt_v1 - wt_v2)

        m1 = load_metric(f"{d1}/md_noWAT_mean.csv")[metric]
        m2 = load_metric(f"{d2}/md_noWAT_mean.csv")[metric]
        m_v1 = window_max(m1, site)
        m_v2 = window_max(m2, site)

        delta1 = m_v1 - wt_v1
        delta2 = m_v2 - wt_v2
        d_avg = ((m_v1 + m_v2) / 2) - ((wt_v1 + wt_v2) / 2)
        agree = (delta1 > 0) == (delta2 > 0)
        robust = agree and abs(d_avg) > noise_floor
        print(f"{name:10s} {d_avg:12.6f} {delta1:12.6f} {delta2:12.6f} {str(agree):>7s} {str(robust):>8s}")

    # T306S-R378K: two sites
    t1 = load_metric("T306S-R378K/md_noWAT_mean.csv")[metric]
    t2 = load_metric("T306S-R378K_2/md_noWAT_mean.csv")[metric]

    t_delta1 = t1.mean() - wt1_m.mean()
    t_delta2 = t2.mean() - wt2_m.mean()
    t_d_avg = ((t1.mean() + t2.mean()) / 2) - ((wt1_m.mean() + wt2_m.mean()) / 2)
    t_agree = (t_delta1 > 0) == (t_delta2 > 0)
    t_robust = t_agree and abs(t_d_avg) > global_noise_floor
    print(f"\n{'T306S-R378K (global)':22s} {t_d_avg:12.6f} {t_delta1:12.6f} {t_delta2:12.6f} {str(t_agree):>7s} {str(t_robust):>8s}")

    for site, label in [(278, "T306"), (350, "R378")]:
        wt_v1 = window_max(wt1_m, site)
        wt_v2 = window_max(wt2_m, site)
        noise_floor = abs(wt_v1 - wt_v2)
        m_v1 = window_max(t1, site)
        m_v2 = window_max(t2, site)
        delta1 = m_v1 - wt_v1
        delta2 = m_v2 - wt_v2
        d_avg = ((m_v1 + m_v2) / 2) - ((wt_v1 + wt_v2) / 2)
        agree = (delta1 > 0) == (delta2 > 0)
        robust = agree and abs(d_avg) > noise_floor
        print(f"T306S-R378K({label}):{'':2s} {d_avg:12.6f} {delta1:12.6f} {delta2:12.6f} {str(agree):>7s} {str(robust):>8s}")
