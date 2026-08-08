#!/usr/bin/env python3
"""
S259R SASA divergence — targeted follow-up (README "Next steps" item 10).

S259R already shows a *robust* global SASA decrease (whole protein more
compact) while its own mutation site shows a robust *increase* in local SASA
(site 231 more exposed). This script digs into that apparent contradiction:

  1. Global total-SASA time series (S259R vs WT, both replicates).
  2. Per-residue mean-SASA delta profile (S259R avg − WT avg) across all
     residues, with the S259R site window (+/-3 around GROMACS 231) shaded.
  3. Local robustness check at site 231 (window-max, the framework metric)
     AND window-mean (less outlier-sensitive), each against the
     WT-replicate noise floor, same convention as rg_sasa_significance_check.py.
  4. Characterization: which residues are most/least exposed in S259R
     relative to WT, to say whether the global compaction is distributed or
     localized, and whether the site-231 exposure change is a real
     site-specific effect rather than part of a large-scale re-arrangement.

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6 inside the cyp2b6 env:
    python s259r_sasa_followup.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import ALLELES_MUTATION_SITES_SIMPLE

ALLELE = "S259R"
REP1, REP2, SITE = ALLELES_MUTATION_SITES_SIMPLE[ALLELE]
TRUE_RESIDUE = SITE + 28  # GROMACS = true - 28
W = 3


def load_sasa_total(path):
    """columns: time, total area (nm^2)"""
    return np.loadtxt(path, comments=["#", "@"])[:, 1]


def load_sasa_res(path):
    """columns: resid (1-based, Protein group, no heme row), mean area (nm^2), std.

    Returns {resid: mean area}. Resid is the file's own column, NOT array
    position -- required because sasa_res_* skips the heme resid 408 (the
    SASA group was "Protein"), so positional indexing would silently shift
    every residue above 408. (Same bug fix as rg_sasa_significance_check.py.)
    """
    data = np.loadtxt(path, comments=["#", "@"])
    return {int(row[0]): row[1] for row in data}


def window_max(res_dict, target, w=W):
    vals = [v for r, v in res_dict.items() if target - w <= r <= target + w]
    return max(vals)


def window_mean(res_dict, target, w=W):
    vals = [v for r, v in res_dict.items() if target - w <= r <= target + w]
    return float(np.mean(vals))


def robust_line(label, d_avg, delta1, delta2, noise_floor, fmt="10.5f"):
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor
    print(f"{label:12s} {d_avg:{fmt}} {delta1:{fmt}} {delta2:{fmt}} "
          f"{str(agree):>7s} {str(robust):>8s}")
    return robust


print("=" * 70)
print(f"{ALLELE} (GROMACS site {SITE} = true {TRUE_RESIDUE}) SASA follow-up")
print("=" * 70)

# ─── 1. Global total SASA ────────────────────────────────────────────────────
s259r1 = load_sasa_total(f"{REP1}/sasa_{REP1}.xvg")
s259r2 = load_sasa_total(f"{REP2}/sasa_{REP2}.xvg")
wt1 = load_sasa_total("WT/sasa_WT.xvg")
wt2 = load_sasa_total("WT_2/sasa_WT_2.xvg")

wt1_m, wt2_m = wt1.mean(), wt2.mean()
global_noise = abs(wt1_m - wt2_m)
d1 = s259r1.mean() - wt1_m
d2 = s259r2.mean() - wt2_m
d_avg = ((s259r1.mean() + s259r2.mean()) / 2) - ((wt1_m + wt2_m) / 2)
print(f"\nGLOBAL total SASA (nm^2). WT noise floor: {global_noise:.4f}")
print(f"{'metric':12s} {'d_avg':>10s} {'rep1 d':>10s} {'rep2 d':>10s} {'agree?':>7s} {'ROBUST?':>8s}")
robust_line("global", d_avg, d1, d2, global_noise)

# ─── 2. Local per-residue SASA at the site window ────────────────────────────
r1 = load_sasa_res(f"{REP1}/sasa_res_{REP1}.xvg")
r2 = load_sasa_res(f"{REP2}/sasa_res_{REP2}.xvg")
w1 = load_sasa_res("WT/sasa_res_WT.xvg")
w2 = load_sasa_res("WT_2/sasa_res_WT_2.xvg")

print(f"\nLOCAL per-residue SASA (nm^2) in the {2 * W + 1}-residue window around site {SITE}")
print(f"{'resid':>6s} {'S259R rep1':>11s} {'S259R rep2':>11s} {'WT rep1':>9s} {'WT rep2':>9s} {'delta':>9s}")
for r in range(SITE - W, SITE + W + 1):
    m = r1.get(r, np.nan)
    m2_ = r2.get(r, np.nan)
    wt1_v = w1.get(r, np.nan)
    wt2_v = w2.get(r, np.nan)
    d = ((m + m2_) / 2) - ((wt1_v + wt2_v) / 2) if not np.isnan(m) else np.nan
    print(f"{r:6d} {m:11.4f} {m2_:11.4f} {wt1_v:9.4f} {wt2_v:9.4f} {d:9.4f}")

for metric, fn in [("window-max", window_max), ("window-mean", window_mean)]:
    wt1_v, wt2_v = fn(w1, SITE), fn(w2, SITE)
    m1_v, m2_v = fn(r1, SITE), fn(r2, SITE)
    noise = abs(wt1_v - wt2_v)
    d1, d2 = m1_v - wt1_v, m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    print(f"\nLOCAL {metric} at site {SITE}. WT noise floor: {noise:.4f}")
    print(f"{'metric':12s} {'d_avg':>10s} {'rep1 d':>10s} {'rep2 d':>10s} {'agree?':>7s} {'ROBUST?':>8s}")
    robust_line(metric, d_avg, d1, d2, noise)

# ─── 3. Per-residue delta profile + residue characterization ────────────────
res_all = sorted(set(r1) & set(r2) & set(w1) & set(w2))
delta_prof = {r: ((r1[r] + r2[r]) / 2) - ((w1[r] + w2[r]) / 2) for r in res_all}
x = np.array(list(delta_prof.keys()))
y = np.array(list(delta_prof.values()))

print("\nTop 10 residues MOST exposed in S259R vs WT (nm^2):")
for r in sorted(delta_prof, key=delta_prof.get, reverse=True)[:10]:
    print(f"  resid {r:4d} (true {r + 28:4d}): +{delta_prof[r]:.4f}")

print("\nTop 10 residues MOST buried in S259R vs WT (nm^2):")
for r in sorted(delta_prof, key=delta_prof.get)[:10]:
    print(f"  resid {r:4d} (true {r + 28:4d}): {delta_prof[r]:.4f}")

site_win = [r for r in res_all if SITE - W <= r <= SITE + W]
win_vals = np.array([delta_prof[r] for r in site_win])
print(f"\nSite-{SITE} window: mean delta {win_vals.mean():.4f}, "
      f"max delta {win_vals.max():.4f} nm^2 "
      f"({'positive (more exposed)' if win_vals.mean() > 0 else 'negative (more buried)'})")

# ─── 4. Figure ───────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

ax1.plot(wt1, color="black", alpha=0.4, linewidth=0.5, label="WT rep1")
ax1.plot(wt2, color="gray", alpha=0.5, linewidth=0.5, label="WT rep2")
ax1.plot(s259r1, color="crimson", alpha=0.7, linewidth=0.6, label=f"{ALLELE} rep1")
ax1.plot(s259r2, color="darkorange", alpha=0.6, linewidth=0.6, label=f"{ALLELE} rep2")
ax1.set_xlabel("Time (ns)", fontsize=10)
ax1.set_ylabel("Total SASA (nm$^2$)", fontsize=10)
ax1.set_title("Global SASA over 300 ns: S259R vs WT (both replicates)",
              fontsize=12, fontweight="bold")
ax1.legend(fontsize=8, loc="lower left")
ax1.grid(True, alpha=0.3)

ax2.axhline(0, color="black", linewidth=0.8)
ax2.plot(x, y, color="steelblue", linewidth=0.8)
ax2.axvspan(SITE - W, SITE + W, color="crimson", alpha=0.15, label=f"site {SITE} window (+/-{W})")
ax2.axvline(SITE, color="crimson", linestyle="--", linewidth=1.2)
ax2.annotate(f"{ALLELE}\nsite {SITE} (true {TRUE_RESIDUE})",
             xy=(SITE, delta_prof.get(SITE, 0)), xytext=(SITE + 12, delta_prof.get(SITE, 0)),
             arrowprops=dict(arrowstyle="->", color="crimson", lw=1), color="crimson", fontsize=9)
ax2.set_xlabel("Residue (GROMACS numbering)", fontsize=10)
ax2.set_ylabel("SASA delta (nm$^2$)", fontsize=10)
ax2.set_title("Per-residue mean SASA delta (S259R avg − WT avg)", fontsize=12, fontweight="bold")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig("s259r_sasa_followup.png", dpi=200, bbox_inches="tight")
print("\nSaved s259r_sasa_followup.png")
