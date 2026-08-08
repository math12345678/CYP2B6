#!/usr/bin/env python3
"""
DRN network-structure deep-dive for the three single-metric DRN-only
findings (README "Next steps" item 9): K139E (local BC), R140Q (local BC),
R487C (local EC + global EC). 

drn_significance_check.py showed each of these alleles picks up exactly one
robust DRN result at its own mutation site (window-max, +/-3), but NOTHING
robust at the global level. Two competing interpretations of such a
"local-only" finding:

  (A) A genuine, localized network reorganization: the mutation site and its
      immediate neighborhood genuinely change centrality, and this is where
      the mutation's dynamical effect concentrates. If the affected residues
      fall in a literature-defined functional region (C/D-loop / proximal
      face / CPR-binding surface for K139E-R140Q; beta-4 / heme-distal
      region for R487C), the finding is structurally meaningful.

  (B) A single-residue artifact: the window-max passes because ONE residue in
      the +/-3 window (often the mutated residue itself) spikes, while the
      rest of the window and protein are unchanged -- the same failure mode
      that s259r_sasa_followup.py diagnosed for the S259R window-max SASA.

This script resolves (A) vs (B) per allele:

  1. Per-residue centrality delta profiles (mutant avg - WT avg) for the
     robust metric of each allele, built from the per-timepoint .dat traces
     (301 timepoints/system, same data md_noWAT_mean.csv averages).
  2. Per-residue WT-replicate disagreement as an empirical noise floor
     (|WT rep1 - WT rep2| at each residue), same convention as the S259R
     SASA follow-up and the framework throughout this project.
  3. Window-max vs window-mean at each mutation site, to test whether the
     finding is a coherent neighborhood effect or a single-residue outlier.
  4. Top-mover census: rank residues by |delta| above the noise floor, and
     annotate each against literature-defined CYP2B6 regions:
       - C/D loop / proximal face (true 136-140; 3IBD HELIX records put the
         C helix at 117-135 and D helix at 141-160; Zhang et al. 2011 locate
         K139 there and show K139E disrupts the P450-CPR complex).
       - SRS-1..SRS-6 (Gotoh 1992; CYP2B6 SRS-2=199-209, SRS-4=290-304,
         SRS-5=360-369 per Lin et al. 2016, DMD 44:1431-1440; SRS-6 around
         beta-4, G478 near SRS-6/beta4-2 per Kaneko et al. 2024).
       - Active-site 5-A set (Angle & Cox 2023, DMD 51:369): I101, I114,
         F115, F206, F297, A298, T302, L363, V367, V477.
  5. Localization summary: fraction of total |delta| mass within +/-3, +/-10,
     +/-20 of the mutation site, and whether robust movers cluster on the
     proximal face / CPR-binding surface.

Numbering: GROMACS resid = true residue - 28 (verified against each system's
md_noWAT_mean_*.cif auth_seq_id and against the mutation placements). Resid
408 is the heme cofactor's CM1 node, excluded (config.RESID_EXCLUDE).

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6 inside the cyp2b6 env:
    python drn_network_deepdive.py
"""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from config import ALLELES_MUTATION_SITES_SIMPLE, RESID_EXCLUDE

# (allele, robust metric) -- taken from drn_significance_check.py output:
#   K139E: local BC robust (+), R140Q: local BC robust (+),
#   R487C: local EC robust (-) AND global EC robust (-).
TARGETS = {
    "K139E": "BC",
    "R140Q": "BC",
    "R487C": "EC",
}

W = 3  # window half-width, same as drn_significance_check.py
N_TOP = 12  # top movers to report per allele


# ---- literature-defined CYP2B6 regions, in TRUE numbering ----
# Convert to GROMACS on the fly (g = t - 28).
TRUE_C_D_LOOP = (136, 140)          # C/D loop / proximal face (Zhang 2011; 3IBD)
TRUE_SRS = {                        # Gotoh 1992; CYP2B6 ranges
    "SRS-1": (96, 117),             #   B' helix + flanking
    "SRS-2": (199, 209),            #   Lin et al. 2016 (F helix)
    "SRS-3": (240, 245),            #   G helix N-terminus
    "SRS-4": (290, 304),            #   Lin et al. 2016 (I helix)
    "SRS-5": (360, 369),            #   Lin et al. 2016 (beta1-4)
    "SRS-6": (469, 478),            #   beta-4 (Kaneko et al. 2024: G478)
}
TRUE_ACTIVE_5A = {101, 114, 115, 206, 297, 298, 302, 363, 367, 477}  # Angle & Cox 2023


def load_trace(sysname, metric):
    """Stack the 301 per-timepoint centrality .dat files -> (time, resid)."""
    files = sorted(glob.glob(os.path.join(sysname, f"md_noWAT_*_{metric}.dat")))
    if not files:
        raise FileNotFoundError(f"no {metric} traces in {sysname}/")
    return np.array([np.loadtxt(f) for f in files])


def true_to_g(t):
    return t - 28


def annotate(g, resnames):
    """Map a GROMACS residue to functional-region labels (true numbering)."""
    t = g + 28
    labels = []
    if TRUE_C_D_LOOP[0] <= t <= TRUE_C_D_LOOP[1]:
        labels.append("C/D loop")
    for name, (lo, hi) in TRUE_SRS.items():
        if lo <= t <= hi:
            labels.append(name)
    if t in TRUE_ACTIVE_5A:
        labels.append("active-5A")
    return ", ".join(labels) if labels else "—"


def analyze(allele, metric):
    rep1, rep2, site = ALLELES_MUTATION_SITES_SIMPLE[allele]
    w1 = load_trace("WT", metric).mean(axis=0)
    w2 = load_trace("WT_2", metric).mean(axis=0)
    m1 = load_trace(rep1, metric).mean(axis=0)
    m2 = load_trace(rep2, metric).mean(axis=0)

    resids = np.arange(1, len(w1) + 1)
    valid = np.array([r not in RESID_EXCLUDE for r in resids])
    resids, w1, w2, m1, m2 = (x[valid] for x in (resids, w1, w2, m1, m2))

    delta = 0.5 * (m1 + m2) - 0.5 * (w1 + w2)
    noise = np.abs(w1 - w2)
    # per-residue robust: both mutant reps agree in sign AND exceed WT noise
    agree = np.sign(m1 - w1) == np.sign(m2 - w2)
    robust = agree & (np.abs(delta) > noise)

    # window max vs window mean at the mutation site
    def wmax(a):
        m = (resids >= site - W) & (resids <= site + W)
        return a[m].max()
    def wmean(a):
        m = (resids >= site - W) & (resids <= site + W)
        return a[m].mean()

    # noise floor for the window (same convention as significance check)
    nf_max = abs(wmax(w1) - wmax(w2))
    nf_mean = abs(wmean(w1) - wmean(w2))
    d_max = 0.5 * (wmax(m1) + wmax(m2)) - 0.5 * (wmax(w1) + wmax(w2))
    d_mean = 0.5 * (wmean(m1) + wmean(m2)) - 0.5 * (wmean(w1) + wmean(w2))
    win_max_agree = (wmax(m1) - wmax(w1)) * (wmax(m2) - wmax(w2)) > 0
    win_mean_agree = (wmean(m1) - wmean(w1)) * (wmean(m2) - wmean(w2)) > 0
    win_max_robust = win_max_agree and abs(d_max) > nf_max
    win_mean_robust = win_mean_agree and abs(d_mean) > nf_mean

    # localization mass
    for r in (W, 10, 20):
        pass
    mass = np.abs(delta).sum()
    fracs = {}
    for r in (W, 10, 20):
        in_win = (resids >= site - r) & (resids <= site + r)
        fracs[r] = float(np.abs(delta)[in_win].sum() / mass)

    # top movers
    idx = np.argsort(-np.abs(delta))[:N_TOP]
    movers = []
    for i in idx:
        movers.append({
            "resid": int(resids[i]),
            "true": int(resids[i] + 28),
            "delta": float(delta[i]),
            "noise": float(noise[i]),
            "ratio": float(delta[i] / noise[i]) if noise[i] > 0 else float("inf"),
            "robust": bool(robust[i]),
            "in_site_window": abs(int(resids[i]) - site) <= W,
            "region": annotate(int(resids[i]), None),
        })

    n_robust = int(robust.sum())
    n_robust_in_site = int((robust & (np.abs(resids - site) <= W)).sum())
    # robust movers on the C/D-loop proximal face (108-112) or SRS-6/beta-4
    cd_lo, cd_hi = true_to_g(TRUE_C_D_LOOP[0]), true_to_g(TRUE_C_D_LOOP[1])
    s6_lo, s6_hi = true_to_g(TRUE_SRS["SRS-6"][0]), true_to_g(TRUE_SRS["SRS-6"][1])
    on_cd = (robust & (resids >= cd_lo) & (resids <= cd_hi)).sum()
    on_s6 = (robust & (resids >= s6_lo) & (resids <= s6_hi)).sum()

    return {
        "allele": allele, "metric": metric, "site": site,
        "resids": resids, "delta": delta, "noise": noise,
        "robust": robust, "agree": agree,
        "d_max": d_max, "nf_max": nf_max, "win_max_robust": win_max_robust,
        "d_mean": d_mean, "nf_mean": nf_mean, "win_mean_robust": win_mean_robust,
        "fracs": fracs, "movers": movers, "n_robust": n_robust,
        "n_robust_in_site": n_robust_in_site, "on_cd": int(on_cd),
        "on_s6": int(on_s6), "n_res": len(resids),
    }


def plot(results, out_png):
    fig, axes = plt.subplots(len(results), 1, figsize=(13, 4.2 * len(results)),
                             sharex=True)
    if len(results) == 1:
        axes = [axes]
    for ax, r in zip(axes, results):
        resids = r["resids"]
        ax.axhline(0, color="0.55", lw=0.8)
        ax.fill_between(resids, -r["noise"], r["noise"], color="0.8", lw=0,
                        label="WT replicate noise floor")
        col = np.where(r["robust"], "#d62728", "#1f77b4")
        ax.bar(resids, r["delta"], color=col, width=1.0)
        # mutation site window
        ax.axvspan(r["site"] - W, r["site"] + W, color="#ffd700", alpha=0.18,
                   label=f"mutation site ±{W}")
        # C/D loop (proximal face) span
        cd_lo, cd_hi = true_to_g(136), true_to_g(140)
        ax.axvspan(cd_lo, cd_hi, color="#2ca02c", alpha=0.12,
                   label="C/D loop (proximal face)")
        # SRS-6 / beta-4 span
        s6_lo, s6_hi = true_to_g(469), true_to_g(478)
        ax.axvspan(s6_lo, s6_hi, color="#9467bd", alpha=0.10, label="SRS-6/β4")
        # label top movers
        for m in r["movers"][:8]:
            ax.annotate(f"{m['true']}", (m["resid"], m["delta"]),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=7, color="#333")
        ax.set_title(f"{r['allele']} — {r['metric']} per-residue Δ "
                     f"(mutant − WT), GROMACS numbering "
                     f"(true = resid + 28); site = {r['site']}",
                     fontsize=11)
        ax.set_ylabel(f"Δ {r['metric']}")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel("GROMACS residue")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved {out_png}")


def main():
    results = []
    for allele, metric in TARGETS.items():
        r = analyze(allele, metric)
        results.append(r)

        print("=" * 100)
        print(f"{r['allele']} ({r['metric']}) — mutation site GROMACS {r['site']} "
              f"(true {r['site'] + 28})")
        print("=" * 100)
        print(f"\nWindow ±{W} at mutation site: "
              f"max Δ = {r['d_max']:+.5f} (noise floor {r['nf_max']:.5f}, "
              f"ROBUST={r['win_max_robust']}); "
              f"mean Δ = {r['d_mean']:+.5f} (noise floor {r['nf_mean']:.5f}, "
              f"ROBUST={r['win_mean_robust']})")
        print(f"  -> window-mean {'PASSES' if r['win_mean_robust'] else 'fails'}: "
              f"{'coherent neighborhood effect' if r['win_mean_robust'] else 'single-residue spike drives the window-max'}")

        print(f"\nRobust per-residue movers: {r['n_robust']}/{r['n_res']} residues "
              f"({100 * r['n_robust'] / r['n_res']:.1f}%), "
              f"of which {r['n_robust_in_site']} in the site ±{W} window; "
              f"{r['on_cd']} on the C/D-loop proximal face (108-112); "
              f"{r['on_s6']} in SRS-6/β4 (441-450)")
        print("Localization of |Δ| mass: "
              + ", ".join(f"±{rr} = {100 * r['fracs'][rr]:.0f}%"
                          for rr in (W, 10, 20)))

        print(f"\nTop {N_TOP} per-residue movers "
              f"(GROMACS / true residue, Δ, Δ/noise, robust, site-window, region):")
        for m in r["movers"]:
            flag = "ROBUST" if m["robust"] else "      "
            sw = "*" if m["in_site_window"] else " "
            print(f"  {m['resid']:3d}/{m['true']:3d}  {m['delta']:+10.5f}  "
                  f"{m['ratio']:8.1f}x  {flag}  {sw}  {m['region']}")
        print()

    plot(results, "drn_network_deepdive.png")


if __name__ == "__main__":
    main()
