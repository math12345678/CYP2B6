"""
Quantitative robustness check for H-bond results, mirroring
significance_check.py's framework for RMSD/RMSF (same rationale: only 2 WT
replicates exist, so their disagreement is used as an empirical noise floor,
and a mutant effect only counts as "robust" if both of its own replicates
agree in direction, the replicate-averaged delta exceeds the noise floor, and
the delta exceeds a minimum absolute effect-size floor).

Two metrics, analogous to RMSD/RMSF:
  - Global: mean H-bonds/frame for the whole system (from hbond_count_*.csv).
  - Local: sum of donor-acceptor pair frequencies where the donor or acceptor
    residue falls within +/-3 residues (GROMACS numbering) of the allele's
    own mutation site (from hbond_pairs_*.csv). Frequency already integrates
    over the whole trajectory (count / n_frames), so summing frequencies
    gives the expected number of site-window H-bonds present in an average
    frame -- the H-bond analogue of a local RMSF magnitude.

Run from ~/Desktop/Research/Research_Projects/CYP2B6 after run_all_hbonds.py.
"""
import numpy as np
import csv

MIN_EFFECT_GLOBAL = 5.0   # H-bonds/frame; ~2% of the ~283 average, well
                          # above typical frame-to-frame noise but below the
                          # WT-replicate noise floor scale seen here (~11)
MIN_EFFECT_LOCAL = 0.15   # summed frequency units at a site window

def load_count(path):
    counts = np.loadtxt(path, delimiter=",", skiprows=1)
    return counts[:, 1]

def load_pairs(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append((int(r["donor_resid"]), int(r["acceptor_resid"]), float(r["frequency"])))
    return rows

# Heme cofactor residues (CM1/HM1/FE1), confirmed via verify_hbond_numbering.py
# to sit at resids 408, 464, 465 in this topology -- close enough to the
# protein's last residue (463) that a mutation-site window near the
# C-terminus (e.g. R487C at site 459, window 456-462) could in principle
# pull in heme atoms. Manual check showed no actual overlap for any current
# allele, but site_window_sum() excludes these resids explicitly rather than
# relying on that margin holding by chance for future sites.
HEME_RESIDS = {408, 464, 465}

def site_window_sum(pairs, target, w=3):
    """A pair counts if either residue is a real (non-heme) residue inside
    the window. Excluding a heme resid from ever counting as "the window
    residue" guards against the near-miss found in verify_hbond_numbering.py
    (R487C's window sits 2 residues from heme resid 464) -- but a window
    residue's H-bond to a heme atom (e.g. a nearby residue contacting a heme
    propionate) is real signal and must still be counted; only the heme side
    of that pair is excluded from *being* the window match, not the pair
    itself."""
    total = 0.0
    for d_resid, a_resid, freq in pairs:
        d_in_window = d_resid not in HEME_RESIDS and (target - w <= d_resid <= target + w)
        a_in_window = a_resid not in HEME_RESIDS and (target - w <= a_resid <= target + w)
        if d_in_window or a_in_window:
            total += freq
    return total

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

# ===== Global: mean H-bonds/frame =====
wt1_counts = load_count("WT/hbond_count_WT.csv")
wt2_counts = load_count("WT_2/hbond_count_WT_2.csv")
wt1_mean, wt2_mean = wt1_counts.mean(), wt2_counts.mean()
global_noise_floor = abs(wt1_mean - wt2_mean)

print(f"Global H-bond noise floor (|WT rep1 mean - WT rep2 mean|): {global_noise_floor:.2f} bonds/frame\n")
print(f"{'Allele':10s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, _) in alleles.items():
    m1 = load_count(f"{d1}/hbond_count_{d1}.csv").mean()
    m2 = load_count(f"{d2}/hbond_count_{d2}.csv").mean()
    delta1 = m1 - wt1_mean
    delta2 = m2 - wt2_mean
    d_avg = ((m1 + m2) / 2) - ((wt1_mean + wt2_mean) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > global_noise_floor and abs(d_avg) > MIN_EFFECT_GLOBAL
    print(f"{name:10s} {d_avg:8.2f} {delta1:8.2f} {delta2:8.2f} {str(agree):>7s} {str(robust):>8s}")

# ===== Local: site-window H-bond frequency sum =====
wt1_pairs = load_pairs("WT/hbond_pairs_WT.csv")
wt2_pairs = load_pairs("WT_2/hbond_pairs_WT_2.csv")

print("\nLocal H-bond significance at each allele's own mutation site (+/-3 residue window)")
print(f"{'Allele':10s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, site) in alleles.items():
    wt1_v = site_window_sum(wt1_pairs, site)
    wt2_v = site_window_sum(wt2_pairs, site)
    noise_floor = abs(wt1_v - wt2_v)

    m1_pairs = load_pairs(f"{d1}/hbond_pairs_{d1}.csv")
    m2_pairs = load_pairs(f"{d2}/hbond_pairs_{d2}.csv")
    m1_v = site_window_sum(m1_pairs, site)
    m2_v = site_window_sum(m2_pairs, site)

    delta1 = m1_v - wt1_v
    delta2 = m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor and abs(d_avg) > MIN_EFFECT_LOCAL
    print(f"{name:10s} {d_avg:8.3f} {delta1:8.3f} {delta2:8.3f} {str(agree):>7s} {str(robust):>8s}")

# T306S-R378K: two sites
t1_pairs = load_pairs("T306S-R378K/hbond_pairs_T306S-R378K.csv")
t2_pairs = load_pairs("T306S-R378K_2/hbond_pairs_T306S-R378K_2.csv")
t1_counts = load_count("T306S-R378K/hbond_count_T306S-R378K.csv")
t2_counts = load_count("T306S-R378K_2/hbond_count_T306S-R378K_2.csv")

t_m1, t_m2 = t1_counts.mean(), t2_counts.mean()
delta1 = t_m1 - wt1_mean
delta2 = t_m2 - wt2_mean
d_avg = ((t_m1 + t_m2) / 2) - ((wt1_mean + wt2_mean) / 2)
agree = (delta1 > 0) == (delta2 > 0)
robust = agree and abs(d_avg) > global_noise_floor and abs(d_avg) > MIN_EFFECT_GLOBAL
print(f"\n{'T306S-R378K (global)':22s} {d_avg:8.2f} {delta1:8.2f} {delta2:8.2f} {str(agree):>7s} {str(robust):>8s}")

for site, label in [(278, "T306"), (350, "R378")]:
    wt1_v = site_window_sum(wt1_pairs, site)
    wt2_v = site_window_sum(wt2_pairs, site)
    noise_floor = abs(wt1_v - wt2_v)
    m1_v = site_window_sum(t1_pairs, site)
    m2_v = site_window_sum(t2_pairs, site)
    delta1 = m1_v - wt1_v
    delta2 = m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor and abs(d_avg) > MIN_EFFECT_LOCAL
    print(f"T306S-R378K({label}):{'':2s} {d_avg:8.3f} {delta1:8.3f} {delta2:8.3f} {str(agree):>7s} {str(robust):>8s}")
