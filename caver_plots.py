#!/usr/bin/env python
"""
Part 8 substrate-access-channel figures from CAVER 3.0 (analog of the
reference paper's substrate-channel-dynamics panels, Rehema et al. JMB 2025):

  - caver_open_frac_all_alleles.png        fraction of snapshots with an open
                                           channel (max bottleneck >= 1.3 A,
                                           water VDW radius)
  - caver_max_bottleneck_all_alleles.png   mean widest tunnel bottleneck
                                           radius per snapshot
  - caver_open_bottleneck_all_alleles.png  mean bottleneck radius over open
                                           snapshots only

Inputs come from caver_analysis/caver_summary_all.csv written by
caver_analysis/analyze_caver_results.py.

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
    with open("caver_analysis/caver_summary_all.csv") as fh:
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

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.2))
    (w1, w2), pairs = pair_values(data, "open_frac")
    bar_with_noise(axes[0], (w1, w2), pairs,
                   "fraction of snapshots with channel >= 1.3 A",
                   "Channel open fraction (water-passable)")
    (w1, w2), pairs = pair_values(data, "mean_max_br")
    bar_with_noise(axes[1], (w1, w2), pairs,
                   "mean max bottleneck radius, all snapshots (A)",
                   "Widest channel bottleneck (all snapshots)")
    (w1, w2), pairs = pair_values(data, "mean_open_br")
    bar_with_noise(axes[2], (w1, w2), pairs,
                   "mean bottleneck radius, open snapshots only (A)",
                   "Open-channel bottleneck")
    fig.tight_layout()
    fig.savefig("caver_open_frac_all_alleles.png", dpi=150)
    plt.close(fig)

    print("Wrote caver_open_frac_all_alleles.png (3 panels: open fraction,")
    print("      max bottleneck, open-channel bottleneck)")


if __name__ == "__main__":
    main()
