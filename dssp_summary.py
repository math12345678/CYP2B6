"""
Summarizes gmx dssp per-frame secondary structure matrices (run_all_dssp_gmx.sh)
into a compact per-residue "ordered secondary structure" fraction, and applies
the same WT-replicate-noise-floor robustness framework used throughout this
project (significance_check.py, hbond/drn/rg_sasa_significance_check.py,
cluster_summary.py). This replaces an earlier cpptraj-based attempt whose
output was unusable (near-zero helix/sheet content everywhere -- see
run_all_dssp.sh / .gitignore comments for why that was abandoned).

gmx dssp -o writes one line per frame, one character per residue, in GROMACS
DSSP one-letter codes:
    H = alpha helix       G = 3-10 helix        I = pi helix
    E = beta strand       B = beta bridge
    T = turn              S = bend              P = kappa/PPII
    ~ = coil/loop (no assigned structure)
    = = gap (residue not assignable, e.g. heme CM1 at resid 408, which lacks
        real backbone atoms)

Column position i (0-indexed) corresponds directly to GROMACS resid i+1 --
the same numbering used for every mutation site elsewhere in this project
(rg_sasa_significance_check.py, drn_significance_check.py, etc).

Metric: "ordered secondary structure" = fraction of frames where a residue
is in {H, G, I, E, B} (any defined helix or strand/bridge element). This
excludes Turn/Bend/PPII/coil, which are not stable fold elements -- a drop
signals local unfolding, a rise signals new/stabilized structure.

  - Global: mean ordered-SS fraction across all 463 real residues (excludes
    resid 408).
  - Local: window-max |delta| at each allele's own mutation site (+/-3
    residues), same windowed approach as RMSF/H-bond/DRN/SASA local checks.

Writes small per-system per-residue summary files (dssp_<SYS>_orderedss.csv,
463 rows, ~5KB each) so the underlying per-residue numbers are trackable in
git without committing the full ~14MB per-frame matrices.

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6 (needs numpy).
"""
import numpy as np

RESID_EXCLUDE = {408}  # heme CM1, appears as '=' (gap) in gmx dssp output
ORDERED_CHARS = set("HGIEB")

SYSTEMS_ALL = [
    "WT", "WT_2",
    "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
    "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2",
    "R140Q", "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2",
    "S259R", "S259R_2", "T306S-R378K", "T306S-R378K_2",
]

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

def ordered_ss_per_residue(sys_name):
    """Returns 1D float array, index i -> resid i+1, ordered-SS fraction
    (NaN for excluded residues)."""
    path = f"{sys_name}/dssp_{sys_name}.dat"
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    arr = np.array([list(l) for l in lines])  # (n_frames, n_res)
    n_res = arr.shape[1]
    ordered = np.isin(arr, list(ORDERED_CHARS))
    frac = ordered.mean(axis=0)
    for excl in RESID_EXCLUDE:
        if 1 <= excl <= n_res:
            frac[excl - 1] = np.nan
    return frac

def window_max_abs_delta(wt_arr, mut_arr, target, w=3):
    lo, hi = max(target - w, 1), min(target + w, len(wt_arr))
    deltas = mut_arr[lo - 1:hi] - wt_arr[lo - 1:hi]
    deltas = deltas[~np.isnan(deltas)]
    if len(deltas) == 0:
        return 0.0
    return deltas[np.argmax(np.abs(deltas))]

print("Computing per-residue ordered-SS fractions for all 24 systems (this reads ~14MB/system, may take a bit)...")
results = {}
for sys in SYSTEMS_ALL:
    frac = ordered_ss_per_residue(sys)
    results[sys] = frac
    # write small per-residue summary (trackable in git)
    out = f"{sys}/dssp_{sys}_orderedss.csv"
    with open(out, "w") as f:
        f.write("resid,ordered_ss_fraction\n")
        for i, v in enumerate(frac, start=1):
            f.write(f"{i},{'' if np.isnan(v) else f'{v:.4f}'}\n")
    print(f"  {sys}: mean={np.nanmean(frac):.4f}  (wrote {out})")

wt1, wt2 = results["WT"], results["WT_2"]
wt1_mean, wt2_mean = np.nanmean(wt1), np.nanmean(wt2)
noise_floor_global = abs(wt1_mean - wt2_mean)

print("\n" + "=" * 70)
print("Ordered secondary structure (H+G+I+E+B) -- global (mean across 463 residues)")
print("=" * 70)
print(f"WT rep1 mean={wt1_mean:.4f}, WT rep2 mean={wt2_mean:.4f}, noise floor={noise_floor_global:.4f}\n")

print(f"{'Allele':16s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
global_robust = {}
for name, (d1, d2, site) in alleles.items():
    m1, m2 = np.nanmean(results[d1]), np.nanmean(results[d2])
    delta1, delta2 = m1 - wt1_mean, m2 - wt2_mean
    d_avg = ((m1 + m2) / 2) - ((wt1_mean + wt2_mean) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor_global
    global_robust[name] = robust
    print(f"{name:16s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}")

t1m, t2m = np.nanmean(results["T306S-R378K"]), np.nanmean(results["T306S-R378K_2"])
delta1, delta2 = t1m - wt1_mean, t2m - wt2_mean
d_avg = ((t1m + t2m) / 2) - ((wt1_mean + wt2_mean) / 2)
agree = (delta1 > 0) == (delta2 > 0)
robust = agree and abs(d_avg) > noise_floor_global
global_robust["T306S-R378K"] = robust
print(f"{'T306S-R378K':16s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}")

print("\n" + "=" * 70)
print("Ordered secondary structure -- local (+/-3 window, max |delta| at own site)")
print("=" * 70)
print(f"{'Allele':16s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
local_robust = {}
for name, (d1, d2, site) in alleles.items():
    delta1 = window_max_abs_delta(wt1, results[d1], site)
    delta2 = window_max_abs_delta(wt2, results[d2], site)
    d_avg = (delta1 + delta2) / 2
    lo, hi = max(site - 3, 1), min(site + 3, len(wt1))
    wt_diffs = wt1[lo - 1:hi] - wt2[lo - 1:hi]
    wt_diffs = wt_diffs[~np.isnan(wt_diffs)]
    noise_floor_local = np.abs(wt_diffs).max() if len(wt_diffs) else 0.0
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor_local
    local_robust[name] = robust
    print(f"{name:16s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}  (noise_floor={noise_floor_local:.4f})")

for site, label in [(278, "T306"), (350, "R378")]:
    delta1 = window_max_abs_delta(wt1, results["T306S-R378K"], site)
    delta2 = window_max_abs_delta(wt2, results["T306S-R378K_2"], site)
    d_avg = (delta1 + delta2) / 2
    lo, hi = max(site - 3, 1), min(site + 3, len(wt1))
    wt_diffs = wt1[lo - 1:hi] - wt2[lo - 1:hi]
    wt_diffs = wt_diffs[~np.isnan(wt_diffs)]
    noise_floor_local = np.abs(wt_diffs).max() if len(wt_diffs) else 0.0
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor_local
    print(f"T306S-R378K({label}):{'':2s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}  (noise_floor={noise_floor_local:.4f})")

print("\nSummary -- alleles with a ROBUST global ordered-SS shift:",
      [a for a, r in global_robust.items() if r] or "none")
print("Summary -- alleles with a ROBUST local (own-site) ordered-SS shift:",
      [a for a, r in local_robust.items() if r] or "none")
