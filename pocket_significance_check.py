#!/usr/bin/env python
"""
Part 8 binding-pocket significance check -- same WT-replicate noise-floor
framework as significance_check.py (RMSD/RMSF), applied to the active-site
pocket metrics aggregated by pocket_summary.py:

  - mean_vol_all       mean active-site pocket volume over all snapshots (A^3)
                       (closed snapshots count as 0 -- combines "how big"
                       with "how often open")
  - mean_vol_open      mean active-site volume over open snapshots only (A^3)
  - open_frac          fraction of snapshots in which the active-site pocket
                       exists (is accessible)
  - heme_hbond_sum     total trajectory frames with a protein-heme H-bond
  - heme_drift_mean    mean heme COM <-> protein COM distance (nm)
  - active_site_rmsf   mean RMSF of the WT active-site lining residues (nm)

As in the established framework, a mutant's effect is only called "robust"
if (a) both mutant replicates agree in sign vs. their paired WT replicate,
(b) the replicate-averaged delta exceeds the WT rep1-vs-rep2 noise floor,
and (c) the delta exceeds a metric-specific minimum effect size (prevents
the class of false positive documented in significance_check.py where a
tiny noise floor lets a physically negligible delta "pass").

Metric-specific MIN_EFFECT values:
  - volumes: 2-3 A^3 -- typical open active-site volumes here are 25-60 A^3,
    so this is ~5-10% of the signal, far above numerical noise (volumes are
    reported to 0.01 A^3).
  - open_frac: 0.03 -- 3% of snapshots (~18 of 601), comparable to the
    sampling granularity of the open/closed dynamics.
  - heme_hbond_sum: 3 frames (of 30001) -- H-bond counts are per-frame
    integers; 3 is a small but non-trivial contact frequency.
  - drift/RMSF: 0.01 nm -- the established MIN_EFFECT_NM from the RMSD/RMSF
    framework.

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.
"""
import csv

import numpy as np

ALLELES = ["G99E", "K139E", "M46V", "I328T", "I391N", "K262R",
           "R140Q", "R487C", "P428T", "S259R", "T306S-R378K"]

METRICS = {
    "mean_vol_all": 2.0,      # A^3
    "mean_vol_open": 3.0,     # A^3
    "open_frac": 0.03,        # fraction
    "heme_hbond_sum": 3.0,    # frames (of 30001)
    "heme_drift_mean": 0.01,  # nm
    "active_site_rmsf": 0.01, # nm
}


def load_summary():
    rows = {}
    with open("pocket_summary_all.csv") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows[r["system"]] = {k: float(v) for k, v in r.items()
                                 if k != "system" and v not in ("", "nan")}
    return rows


def main():
    data = load_summary()
    wt1, wt2 = data["WT"], data["WT_2"]

    for metric, min_effect in METRICS.items():
        if metric not in wt1 or metric not in wt2:
            print(f"\n{metric}: missing in summary -- skipping")
            continue
        noise_floor = abs(wt1[metric] - wt2[metric])
        print(f"\n{metric}  (noise floor |WT1-WT2| = {noise_floor:.4f}, "
              f"min effect = {min_effect})")
        print(f"{'Allele':12s} {'d1':>8s} {'d2':>8s} {'d_avg':>8s} "
              f"{'agree?':>7s} {'ROBUST?':>9s}")
        for a in ALLELES:
            r1, r2 = data.get(a), data.get(f"{a}_2")
            if not r1 or not r2 or metric not in r1 or metric not in r2:
                print(f"{a:12s}  (incomplete data)")
                continue
            d1 = r1[metric] - wt1[metric]
            d2 = r2[metric] - wt2[metric]
            d_avg = (d1 + d2) / 2
            agree = (d1 > 0) == (d2 > 0)
            robust = agree and abs(d_avg) > noise_floor and abs(d_avg) > min_effect
            print(f"{a:12s} {d1:8.4f} {d2:8.4f} {d_avg:8.4f} "
                  f"{str(agree):>7s} {str(robust):>9s}")


if __name__ == "__main__":
    main()
