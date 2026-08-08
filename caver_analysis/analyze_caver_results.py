#!/usr/bin/env python
"""Analyze CAVER 3.0 tunnel results across all CYP2B6 systems.

Reads caver_results/<SYS>/analysis/tunnel_characteristics.csv (one row per
tunnel per snapshot) and computes per-snapshot tunnel metrics:

  - n_tunnels        : number of tunnels detected in the snapshot
  - max_br           : maximum bottleneck radius over tunnels (A)
  - max_len          : maximum tunnel length (A)
  - open(>=1.3)      : fraction of snapshots with a tunnel wider than the
                       water VDW radius (1.3 A) - 'open' channel

These are aggregated per system (mean over snapshots), then compared across
alleles with the same WT-replicate noise-floor framework used elsewhere:
a metric is 'robust' for an allele if both replicates move the same direction
relative to WT and the replicate-averaged delta exceeds the WT noise floor
(|delta| > max(abs(WT_rep1 - WT_rep2))).

Also reports the top tunnel cluster per system (highest No_snaps) and its
bottleneck-lining residues from bottlenecks.csv (GROMACS numbering).
"""
import csv
import os
import sys
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(BASE, "caver_results")

SYSTEMS = ["WT", "WT_2", "G99E", "G99E_2", "K139E", "K139E_2", "M46V",
           "M46V_2", "I328T", "I328T_2", "I391N", "I391N_2", "K262R",
           "K262R_2", "R140Q", "R140Q_2", "R487C", "R487C_2", "P428T",
           "P428T_2", "S259R", "S259R_2", "T306S-R378K", "T306S-R378K_2"]

OPEN_BR = 1.3  # water van der Waals radius (A)


def parse_characteristics(sysname):
    path = os.path.join(RESULTS, sysname, "analysis", "tunnel_characteristics.csv")
    if not os.path.isfile(path):
        return None
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [fn.strip() for fn in reader.fieldnames]
        for r in reader:
            try:
                rows.append({
                    "snap": int(r["Snapshot"].split("_")[1].split(".")[0]),
                    "cluster": int(r["Tunnel cluster"]),
                    "br": float(r["Bottleneck radius"]),
                    "length": float(r["Length"]),
                })
            except (ValueError, KeyError):
                continue
    return rows


def system_metrics(sysname):
    rows = parse_characteristics(sysname)
    if not rows:
        return None
    by_snap = defaultdict(list)
    for r in rows:
        by_snap[r["snap"]].append(r)
    snap_ids = sorted(set(range(1, 602, 5)) | set(by_snap.keys()))
    n_frames = len(snap_ids)
    n_tunnels = []
    max_br = []
    max_len = []
    for s in snap_ids:
        snaps = by_snap.get(s, [])
        n_tunnels.append(len(snaps))
        max_br.append(max((x["br"] for x in snaps), default=0.0))
        max_len.append(max((x["length"] for x in snaps), default=0.0))
    open_brs = [v for v in max_br if v >= OPEN_BR]
    return {
        "mean_n_tunnels": sum(n_tunnels) / n_frames,
        "mean_max_br": sum(max_br) / n_frames,
        "mean_open_br": (sum(open_brs) / len(open_brs)) if open_brs else 0.0,
        "mean_max_len": sum(max_len) / n_frames,
        "open_frac": len(open_brs) / n_frames,
        "open_frac_17": sum(1 for v in max_br if v >= 1.7) / n_frames,
    }


def top_cluster(sysname):
    """Return (cluster_id, No_snaps, Avg_BR, lining_residues) for the top cluster."""
    path = os.path.join(RESULTS, sysname, "analysis", "bottlenecks.csv")
    residues = {}
    if os.path.isfile(path):
        with open(path) as f:
            for line in f:
                if line.startswith("Snapshot,"):
                    continue
                parts = line.rstrip().split(",")
                if len(parts) < 2:
                    continue
                try:
                    cluster = int(parts[1])
                except ValueError:
                    continue
                if len(parts) > 9 and parts[9]:
                    res = [x for x in parts[9].split(":")[1:] if x]
                    residues.setdefault(cluster, set()).update(int(x) for x in res)
    # top cluster by tunnel count from characteristics
    rows = parse_characteristics(sysname) or []
    counts = defaultdict(int)
    for r in rows:
        counts[r["cluster"]] += 1
    if not counts:
        return None
    top = max(counts, key=lambda c: (counts[c], c))
    return {
        "cluster": top,
        "n_tunnels": counts[top],
        "residues": sorted(residues.get(top, [])),
    }


def main():
    metrics = {}
    tops = {}
    for s in SYSTEMS:
        m = system_metrics(s)
        if m is None:
            print(f"MISSING {s}")
            continue
        metrics[s] = m
        tops[s] = top_cluster(s)

    wt_rep1 = metrics["WT"]
    wt_rep2 = metrics["WT_2"]
    noise = {}
    for k in wt_rep1:
        noise[k] = abs(wt_rep1[k] - wt_rep2[k])

    out_csv = os.path.join(BASE, "caver_summary_all.csv")
    with open(out_csv, "w", newline="") as f:
        fieldnames = ["system"] + list(metrics["WT"].keys())
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for s in SYSTEMS:
            if s in metrics:
                w.writerow({"system": s, **metrics[s]})
    print(f"Wrote {out_csv}")

    print("=" * 96)
    print("CAVER 3.0 tunnel metrics per system (121 snapshots each)")
    print("=" * 96)
    hdr = (f"{'system':<16}{'n_tun/snap':>12}{'mean_maxBR':>12}"
           f"{'mean_openBR':>13}{'max_len':>9}{'open>=1.3':>10}{'open>=1.7':>10}")
    print(hdr)
    print("-" * 96)
    for s in SYSTEMS:
        m = metrics.get(s)
        if not m:
            continue
        print(f"{s:<16}{m['mean_n_tunnels']:>12.2f}{m['mean_max_br']:>12.3f}"
              f"{m['mean_open_br']:>13.3f}{m['mean_max_len']:>9.1f}"
              f"{m['open_frac']:>10.3f}{m['open_frac_17']:>10.3f}")

    print()
    print("WT noise floor (|WT_rep1 - WT_rep2|):")
    for k, v in noise.items():
        print(f"  {k}: {v:.4f}")

    print()
    print("=" * 88)
    print("Robust allele deltas vs WT (both replicates same sign, "
          "|delta| > WT noise floor)")
    print("=" * 88)
    alleles = [s for s in SYSTEMS if not s.endswith("_2") and s != "WT"]
    for metric in ["mean_n_tunnels", "mean_max_br", "mean_open_br",
                   "open_frac", "open_frac_17"]:
        print(f"\n--- {metric} (WT noise floor {noise[metric]:.4f}) ---")
        for a in alleles:
            a1 = metrics.get(a)
            a2 = metrics.get(a + "_2")
            if not a1 or not a2:
                print(f"  {a:<14} incomplete")
                continue
            d1 = a1[metric] - wt_rep1[metric]
            d2 = a2[metric] - wt_rep2[metric]
            mean_d = (d1 + d2) / 2
            same_sign = (d1 * d2) > 0
            robust = same_sign and abs(mean_d) > noise[metric]
            flag = "ROBUST" if robust else ("same-sign" if same_sign else "opposite-sign")
            print(f"  {a:<14} rep1={d1:+.4f} rep2={d2:+.4f} mean={mean_d:+.4f} -> {flag}")

    print()
    print("=" * 88)
    print("Top tunnel cluster per system (highest tunnel count) - "
          "lining residues (GROMACS numbering)")
    print("=" * 88)
    for s in SYSTEMS:
        t = tops.get(s)
        if not t:
            continue
        res = ",".join(str(r) for r in t["residues"][:20])
        print(f"{s:<16}cluster {t['cluster']:>3}  n={t['n_tunnels']:>3}  "
              f"residues: {res}")


if __name__ == "__main__":
    main()
