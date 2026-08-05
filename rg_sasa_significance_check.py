"""
Quantitative robustness check for radius of gyration (Rg) and solvent-
accessible surface area (SASA), mirroring the same WT-replicate-noise-floor
framework used throughout (significance_check.py, hbond_significance_check.py,
drn_significance_check.py).

Rg is a whole-molecule scalar (no per-residue breakdown), so only a global
metric applies: mean total Rg over the trajectory (gyrate_<SYSTEM>.xvg column
2, computed on the Protein group -- see run_all_rg_sasa.sh for rationale).

SASA has both:
  - Global: mean total SASA over the trajectory (sasa_<SYSTEM>.xvg column 2).
  - Local: per-residue average SASA (sasa_res_<SYSTEM>.xvg), window-max at
    each allele's own mutation site (+/-3 residues), the SASA analogue of the
    windowed RMSF/H-bond/DRN checks. sasa_res_* has exactly 462 rows (one per
    real amino acid; unlike the DRN CSVs, no heme row here since the SASA
    calculation group was "Protein", which excludes the heme cofactor
    entirely -- confirmed via row count, no exclusion needed).

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6 after
run_all_rg_sasa.sh, inside the cyp2b6 env.
"""
import numpy as np

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

def load_gyrate(path):
    # columns: time, Rg(total), Rg(X), Rg(Y), Rg(Z)
    data = np.loadtxt(path, comments=["#", "@"])
    return data[:, 1]

def load_sasa_total(path):
    # columns: time, total area (nm^2)
    data = np.loadtxt(path, comments=["#", "@"])
    return data[:, 1]

def load_sasa_res(path):
    # columns: residue index (1-based, Protein group only -- this is the
    # real GROMACS resid, and it SKIPS 408 (heme, not part of "Protein"),
    # e.g. row order goes ...406 407 409 410... with no row for 408.
    # BUG FIX (found during Part 7/PCA-DCCM work): originally this returned
    # a plain positional array (data[:, 1]) and window_max() assumed index
    # i == resid i+1, silently ignoring the file's own resid column. That
    # is wrong for any mutation site with resid > 408, since the missing
    # heme row shifts every later residue's *array position* one slot
    # earlier than its *resid* -- confirmed directly by inspecting
    # sasa_res_WT.xvg, which jumps "407 ... 409" with no 408 row. Of this
    # project's mutation sites, only R487C (resid 459) is affected. Fixed
    # by returning a {resid: value} dict and looking up by actual resid
    # instead of by array position.
    data = np.loadtxt(path, comments=["#", "@"])
    return {int(row[0]): row[1] for row in data}

def window_max(res_dict, target, w=3):
    vals = [v for r, v in res_dict.items() if target - w <= r <= target + w]
    return max(vals)

print("=" * 70)
print("Radius of gyration (Rg) -- global only (whole-molecule scalar)")
print("=" * 70)

wt1_rg = load_gyrate("WT/gyrate_WT.xvg")
wt2_rg = load_gyrate("WT_2/gyrate_WT_2.xvg")
wt1_mean, wt2_mean = wt1_rg.mean(), wt2_rg.mean()
rg_noise_floor = abs(wt1_mean - wt2_mean)
print(f"WT noise floor: {rg_noise_floor:.5f} nm\n")
print(f"{'Allele':10s} {'d_avg':>10s} {'rep1 d':>10s} {'rep2 d':>10s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, _) in alleles.items():
    m1 = load_gyrate(f"{d1}/gyrate_{d1}.xvg").mean()
    m2 = load_gyrate(f"{d2}/gyrate_{d2}.xvg").mean()
    delta1 = m1 - wt1_mean
    delta2 = m2 - wt2_mean
    d_avg = ((m1 + m2) / 2) - ((wt1_mean + wt2_mean) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > rg_noise_floor
    print(f"{name:10s} {d_avg:10.5f} {delta1:10.5f} {delta2:10.5f} {str(agree):>7s} {str(robust):>8s}")

t1 = load_gyrate("T306S-R378K/gyrate_T306S-R378K.xvg").mean()
t2 = load_gyrate("T306S-R378K_2/gyrate_T306S-R378K_2.xvg").mean()
delta1 = t1 - wt1_mean
delta2 = t2 - wt2_mean
d_avg = ((t1 + t2) / 2) - ((wt1_mean + wt2_mean) / 2)
agree = (delta1 > 0) == (delta2 > 0)
robust = agree and abs(d_avg) > rg_noise_floor
print(f"{'T306S-R378K':10s} {d_avg:10.5f} {delta1:10.5f} {delta2:10.5f} {str(agree):>7s} {str(robust):>8s}")

print("\n" + "=" * 70)
print("SASA -- global (total area)")
print("=" * 70)

wt1_sasa = load_sasa_total("WT/sasa_WT.xvg")
wt2_sasa = load_sasa_total("WT_2/sasa_WT_2.xvg")
wt1_mean, wt2_mean = wt1_sasa.mean(), wt2_sasa.mean()
sasa_noise_floor = abs(wt1_mean - wt2_mean)
print(f"WT noise floor: {sasa_noise_floor:.3f} nm^2\n")
print(f"{'Allele':10s} {'d_avg':>10s} {'rep1 d':>10s} {'rep2 d':>10s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, _) in alleles.items():
    m1 = load_sasa_total(f"{d1}/sasa_{d1}.xvg").mean()
    m2 = load_sasa_total(f"{d2}/sasa_{d2}.xvg").mean()
    delta1 = m1 - wt1_mean
    delta2 = m2 - wt2_mean
    d_avg = ((m1 + m2) / 2) - ((wt1_mean + wt2_mean) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > sasa_noise_floor
    print(f"{name:10s} {d_avg:10.3f} {delta1:10.3f} {delta2:10.3f} {str(agree):>7s} {str(robust):>8s}")

t1 = load_sasa_total("T306S-R378K/sasa_T306S-R378K.xvg").mean()
t2 = load_sasa_total("T306S-R378K_2/sasa_T306S-R378K_2.xvg").mean()
delta1 = t1 - wt1_mean
delta2 = t2 - wt2_mean
d_avg = ((t1 + t2) / 2) - ((wt1_mean + wt2_mean) / 2)
agree = (delta1 > 0) == (delta2 > 0)
robust = agree and abs(d_avg) > sasa_noise_floor
print(f"{'T306S-R378K':10s} {d_avg:10.3f} {delta1:10.3f} {delta2:10.3f} {str(agree):>7s} {str(robust):>8s}")

print("\n" + "=" * 70)
print("SASA -- local (per-residue, +/-3 window max at each allele's own site)")
print("=" * 70)

wt1_res = load_sasa_res("WT/sasa_res_WT.xvg")
wt2_res = load_sasa_res("WT_2/sasa_res_WT_2.xvg")

print(f"{'Allele':10s} {'d_avg':>10s} {'rep1 d':>10s} {'rep2 d':>10s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, site) in alleles.items():
    wt1_v = window_max(wt1_res, site)
    wt2_v = window_max(wt2_res, site)
    noise_floor = abs(wt1_v - wt2_v)

    m1_res = load_sasa_res(f"{d1}/sasa_res_{d1}.xvg")
    m2_res = load_sasa_res(f"{d2}/sasa_res_{d2}.xvg")
    m1_v = window_max(m1_res, site)
    m2_v = window_max(m2_res, site)

    delta1 = m1_v - wt1_v
    delta2 = m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor
    print(f"{name:10s} {d_avg:10.4f} {delta1:10.4f} {delta2:10.4f} {str(agree):>7s} {str(robust):>8s}")

t1_res = load_sasa_res("T306S-R378K/sasa_res_T306S-R378K.xvg")
t2_res = load_sasa_res("T306S-R378K_2/sasa_res_T306S-R378K_2.xvg")
for site, label in [(278, "T306"), (350, "R378")]:
    wt1_v = window_max(wt1_res, site)
    wt2_v = window_max(wt2_res, site)
    noise_floor = abs(wt1_v - wt2_v)
    m1_v = window_max(t1_res, site)
    m2_v = window_max(t2_res, site)
    delta1 = m1_v - wt1_v
    delta2 = m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > noise_floor
    print(f"T306S-R378K({label}):{'':2s} {d_avg:10.4f} {delta1:10.4f} {delta2:10.4f} {str(agree):>7s} {str(robust):>8s}")
