#!/usr/bin/env python
"""
Part 8 binding-pocket figures (analogs of the reference paper's active-site
volume and heme-contact panels, Rehema et al. JMB 2025):

  - pocket_volume_all_alleles.png      mean active-site volume + open
                                       fraction per allele, with the
                                       WT-replicate noise floor band
  - pocket_volume_timeseries_WT.png    per-snapshot active-site volume for
                                       both WT replicates (Fig 4D analog)
  - pocket_heme_contacts_all_alleles.png  protein-heme H-bond contact
                                       counts per allele (Fig 4E analog)
  - active_site_rmsf_all_alleles.png   mean RMSF of the active-site lining
                                       residues per allele

All inputs come from pocket_summary_all.csv and the per-system
pocket_<SYS>.csv time series written by pocket_summary.py.

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ALLELES = ["G99E", "K139E", "M46V", "I328T", "I391N", "K262R",
           "R140Q", "R487C", "P428T", "S259R", "T306S-R378K"]
WT_BLUE = "#4477AA"
MUT_RED = "#BB5566"


def load_summary():
    rows = {}
    with open("pocket_summary_all.csv") as fh:
        for r in csv.DictReader(fh):
            rows[r["system"]] = {k: (float(v) if v not in ("", "nan") else np.nan)
                                 for k, v in r.items() if k != "system"}
    return rows


def pair_values(data, metric):
    """Replicate-paired values: (WT1, WT2), then [(allele rep1, rep2), ...]"""
    w1, w2 = data["WT"][metric], data["WT_2"][metric]
    pairs = [(data[a][metric], data[f"{a}_2"][metric]) for a in ALLELES]
    return (w1, w2), pairs


def bar_with_noise(ax, w_pair, pairs, ylabel, title, min_effect=0.0):
    """Bars = replicate means per allele; error = rep1..rep2 range.
    Shaded band = WT-replicate disagreement (noise floor)."""
    noise_floor = abs(w_pair[0] - w_pair[1])
    means = np.array([np.nanmean(p) for p in pairs])
    lows = np.array([min(p) for p in pairs])
    highs = np.array([max(p) for p in pairs])
    x = np.arange(len(ALLELES))
    ax.bar(x, means, color=MUT_RED, alpha=0.85, width=0.6,
           yerr=[means - lows, highs - means], capsize=3,
           error_kw={"elinewidth": 1, "capthick": 1})
    ax.axhline(np.mean([w_pair[0], w_pair[1]]), color=WT_BLUE, ls="--", lw=1.2,
               label="WT mean")
    ax.fill_between(x, np.mean([w_pair[0], w_pair[1]]) - noise_floor,
                    np.mean([w_pair[0], w_pair[1]]) + noise_floor,
                    color=WT_BLUE, alpha=0.15, label="WT rep1/rep2 spread")
    ax.set_xticks(x)
    ax.set_xticklabels(ALLELES, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="best")


def main():
    data = load_summary()

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
    (w1, w2), pairs = pair_values(data, "mean_vol_all")
    bar_with_noise(axes[0], (w1, w2), pairs,
                   "mean active-site volume, all snapshots (A^3)",
                   "Active-site pocket volume (closed snapshots = 0)")
    (w1, w2), pairs = pair_values(data, "open_frac")
    bar_with_noise(axes[1], (w1, w2), pairs,
                   "fraction of snapshots with open pocket",
                   "Active-site open fraction")
    fig.tight_layout()
    fig.savefig("pocket_volume_all_alleles.png", dpi=150)
    plt.close(fig)

    # WT volume time series (both replicates), rolling mean -- Fig 4D analog
    fig, ax = plt.subplots(figsize=(11, 3.6))
    for sys, col in [("WT", "#4477AA"), ("WT_2", "#88AACC")]:
        snap, vol = np.loadtxt(f"{sys}/pocket_{sys}.csv", delimiter=",", skiprows=1).T
        ax.plot(snap, vol, color=col, lw=0.5, alpha=0.35, label=f"{sys}")
        win = 25
        roll = np.convolve(vol, np.ones(win) / win, mode="same")
        ax.plot(snap, roll, color=col, lw=1.8, label=f"{sys} rolling mean (n={win})")
    ax.set_xlabel("snapshot (0.5 ns per frame)")
    ax.set_ylabel("active-site pocket volume (A^3)")
    ax.set_title("WT active-site pocket volume dynamics")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("pocket_volume_timeseries_WT.png", dpi=150)
    plt.close(fig)

    # Heme contacts -- Fig 4E analog
    fig, ax = plt.subplots(figsize=(9, 3.8))
    (w1, w2), pairs = pair_values(data, "heme_hbond_sum")
    bar_with_noise(ax, (w1, w2), pairs,
                   "total protein-heme H-bond frames (of 30001)",
                   "Heme H-bond contacts")
    fig.tight_layout()
    fig.savefig("pocket_heme_contacts_all_alleles.png", dpi=150)
    plt.close(fig)

    # Active-site RMSF
    fig, ax = plt.subplots(figsize=(9, 3.8))
    (w1, w2), pairs = pair_values(data, "active_site_rmsf")
    bar_with_noise(ax, (w1, w2), pairs,
                   "mean RMSF of active-site lining residues (nm)",
                   "Active-site flexibility")
    fig.tight_layout()
    fig.savefig("active_site_rmsf_all_alleles.png", dpi=150)
    plt.close(fig)

    print("Wrote pocket_volume_all_alleles.png, pocket_volume_timeseries_WT.png,")
    print("      pocket_heme_contacts_all_alleles.png, active_site_rmsf_all_alleles.png")


if __name__ == "__main__":
    main()
