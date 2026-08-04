"""
Summarizes cpptraj conformational clustering results (run_all_clustering.sh)
across all 22 systems, analogous in spirit to Shaylyn Govender's thesis Table
4.1 (number of clusters and per-cluster frame fractions), so results can be
compared directly against her prior findings on this same system.

Metric used for cross-replicate/WT comparison: the largest cluster's frame
fraction ("dominant cluster fraction"). A value close to 1.0 means the
trajectory is dominated by one conformational state (low conformational
diversity); a lower value spread across multiple substantial clusters means
the system samples more than one distinct conformation. This is the same
quantity Shaylyn's Table 4.1 reports per system, just re-derived here from
our own (Backbone RMSD-only, no L-helix/heme-Fe metric -- see
run_all_clustering.sh for why) clustering run.

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6, inside the cyp2b6
env (only needs the standard library).
"""
import re

SYSTEMS = [
    "WT", "WT_2",
    "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
    "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2",
    "R140Q", "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2",
    "S259R", "S259R_2", "T306S-R378K", "T306S-R378K_2",
]

def parse_summary(path):
    """Returns list of (cluster_id, n_frames, frac) sorted by frac descending."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            # #Cluster Frames Frac AvgDist Stdev Centroid AvgCDist
            cluster_id, n_frames, frac = int(parts[0]), int(parts[1]), float(parts[2])
            rows.append((cluster_id, n_frames, frac))
    rows.sort(key=lambda r: -r[2])
    return rows

print(f"{'System':16s} {'#Clusters':>9s} {'Dominant frac':>14s} {'Fracs (desc)':>30s}")
results = {}
for sys in SYSTEMS:
    path = f"{sys}/cluster_{sys}_summary.dat"
    rows = parse_summary(path)
    n_clusters = len(rows)
    dominant_frac = rows[0][2]
    fracs_str = ", ".join(f"{r[2]:.3f}" for r in rows)
    results[sys] = (n_clusters, dominant_frac, rows)
    print(f"{sys:16s} {n_clusters:9d} {dominant_frac:14.3f} {fracs_str:>30s}")

# Same WT-replicate-noise-floor robustness framework used throughout this
# project (significance_check.py, hbond/drn/rg_sasa_significance_check.py):
# WT rep1 vs rep2 disagreement sets the noise floor, and an allele's effect
# only counts as robust if both of its replicates move in the same direction
# relative to their respective WT replicate AND the replicate-averaged delta
# exceeds that noise floor. This matters a lot here: the two WT replicates
# disagree sharply on dominant-cluster-fraction (0.454 vs 0.808), so this is
# a noisy metric and most raw deltas below should NOT be read as real effects
# without clearing that bar.
alleles = ["G99E", "K139E", "M46V", "I328T", "I391N", "K262R",
           "R140Q", "R487C", "P428T", "S259R"]

wt1, wt2 = results["WT"][1], results["WT_2"][1]
noise_floor = abs(wt1 - wt2)
print(f"\nWT rep1={wt1:.3f}, WT rep2={wt2:.3f}, noise floor={noise_floor:.3f}\n")

print(f"{'Allele':16s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
for a in alleles:
    r1, r2 = results[a][1], results[f"{a}_2"][1]
    delta1, delta2 = r1 - wt1, r2 - wt2
    d_avg = ((r1 + r2) / 2) - ((wt1 + wt2) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor
    print(f"{a:16s} {d_avg:8.3f} {delta1:8.3f} {delta2:8.3f} {str(agree):>7s} {str(robust):>8s}")

t1, t2 = results["T306S-R378K"][1], results["T306S-R378K_2"][1]
delta1, delta2 = t1 - wt1, t2 - wt2
d_avg = ((t1 + t2) / 2) - ((wt1 + wt2) / 2)
agree = (delta1 > 0) == (delta2 > 0)
robust = agree and abs(d_avg) > noise_floor
print(f"{'T306S-R378K':16s} {d_avg:8.3f} {delta1:8.3f} {delta2:8.3f} {str(agree):>7s} {str(robust):>8s}")

print("\nFor reference, Shaylyn's thesis (Table 4.1, cpptraj RMSD + heme Fe +")
print("L-helix metric, not directly reproduced here) flagged I328T, K262R,")
print("P428T, and R140Q as showing the most cluster/structural deviation.")
