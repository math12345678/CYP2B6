"""
Summarizes PCA (essential dynamics) and DCCM (dynamic cross-correlation
matrix) output from run_all_pca_dccm.sh, applying the same WT-replicate-
noise-floor robustness framework used throughout this project.

PCA: eigenval_<SYSTEM>.xvg holds the 4158 eigenvalues of the Backbone
covariance matrix (nm^2), sorted descending. Two global metrics are used:
  - PC1 eigenvalue (nm^2): the variance captured by the single dominant
    motion -- a larger value means more of the trajectory's flexibility is
    concentrated into one collective motion rather than spread out.
  - PC1 fraction: PC1 eigenvalue / sum of all eigenvalues -- how dominant
    that single mode is relative to everything else, independent of the
    system's total flexibility.
PCA has no natural per-residue local metric in this analysis (that would
require examining eigenvector *components*, a further step not attempted
here), so only global comparisons are made.

DCCM: covar_<SYSTEM>.dat is the raw ascii covariance matrix (nm^2) for the
C-alpha selection (462 atoms = 462 real residues; heme has no atom named
"CA" so it is automatically and completely excluded from this group, unlike
some earlier per-residue datasets in this project that included a heme row).
The file is a flattened (3*462) x (3*462) matrix, 3 values/line. This is
reduced to a standard 462x462 residue-residue correlation matrix:
    C_ij = trace(Cov_block_ij) / sqrt(trace(Cov_ii) * trace(Cov_jj))
where Cov_block_ij is the 3x3 sub-block of x/y/z covariances between the
C-alpha atoms of residues i and j. This matrix is NOT written to disk (it
would be a large, awkward file to track); only derived summary numbers are
kept.

  - Global: mean |correlation| across all off-diagonal residue pairs --
    a measure of how tightly coupled the whole protein's motion is overall.
  - Local: mean |correlation| between each allele's own mutation-site
    residue and every OTHER residue more than 3 positions away (excluding
    immediate neighbors, whose high correlation is a trivial consequence of
    backbone connectivity, not a meaningful long-range coupling signal) --
    a proxy for how strongly that site's motion is coupled to the rest of
    the protein (an "allosteric coupling strength" metric).

IMPORTANT residue-numbering note: the C-alpha group (462 atoms, no heme) is
numbered contiguously 1-462 with NO gap at the heme's position (unlike the
463-464-row datasets elsewhere in this project, e.g. DRN/DSSP, whose
numbering runs 1-463 THROUGH a heme placeholder at resid 408). Any mutation
site at a GROMACS resid greater than 408 must have 1 subtracted before
indexing into this C-alpha array. Of this project's mutation sites, only
R487C's own site (resid 459) is affected; T306S-R378K's sites (278, 350)
are both below 408 and unaffected. This adjustment is applied below via
`ca_index()`. (This is worth a quick sanity check against the earlier
Rg/SASA local-metric scripts, which may not have made this same adjustment
for R487C specifically -- flagged in Next steps, not re-litigated here.)

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6 (needs numpy).
"""
import numpy as np

N_RES = 462  # C-alpha group size, no heme gap
HEME_RESID = 408  # combined (backbone/DRN/DSSP) numbering position of heme

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

def ca_index(resid):
    """Maps a GROMACS resid (combined/heme-inclusive numbering, 1-indexed)
    to a 0-indexed position in the 462-atom, heme-free C-alpha array."""
    idx1 = resid - 1 if resid > HEME_RESID else resid
    return idx1 - 1  # to 0-indexed

# ---------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------
def load_eigenvalues(path):
    data = np.loadtxt(path, comments=["#", "@"])
    return data[:, 1]  # nm^2, descending

print("=" * 70)
print("PCA -- PC1 eigenvalue and PC1 variance fraction (global only)")
print("=" * 70)

pca = {}
for sys in SYSTEMS_ALL:
    eig = load_eigenvalues(f"{sys}/eigenval_{sys}.xvg")
    pc1 = eig[0]
    frac = eig[0] / eig.sum()
    pca[sys] = (pc1, frac)

wt1_pc1, wt1_frac = pca["WT"]
wt2_pc1, wt2_frac = pca["WT_2"]
nf_pc1 = abs(wt1_pc1 - wt2_pc1)
nf_frac = abs(wt1_frac - wt2_frac)
print(f"WT rep1 PC1={wt1_pc1:.4f} nm^2 ({wt1_frac:.4f}), "
      f"WT rep2 PC1={wt2_pc1:.4f} nm^2 ({wt2_frac:.4f})")
print(f"Noise floor: PC1 eigenvalue={nf_pc1:.4f}, PC1 fraction={nf_frac:.4f}\n")

print(f"{'Allele':16s} {'PC1 d_avg':>10s} {'agree?':>7s} {'ROBUST?':>8s}  "
      f"{'Frac d_avg':>10s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, site) in alleles.items():
    m1, f1 = pca[d1]
    m2, f2 = pca[d2]
    d1_pc1, d2_pc1 = m1 - wt1_pc1, m2 - wt2_pc1
    davg_pc1 = ((m1 + m2) / 2) - ((wt1_pc1 + wt2_pc1) / 2)
    agree_pc1 = (d1_pc1 > 0) == (d2_pc1 > 0)
    robust_pc1 = agree_pc1 and abs(davg_pc1) > nf_pc1

    d1_f, d2_f = f1 - wt1_frac, f2 - wt2_frac
    davg_f = ((f1 + f2) / 2) - ((wt1_frac + wt2_frac) / 2)
    agree_f = (d1_f > 0) == (d2_f > 0)
    robust_f = agree_f and abs(davg_f) > nf_frac

    print(f"{name:16s} {davg_pc1:10.4f} {str(agree_pc1):>7s} {str(robust_pc1):>8s}  "
          f"{davg_f:10.4f} {str(agree_f):>7s} {str(robust_f):>8s}")

t1m, t1f = pca["T306S-R378K"]
t2m, t2f = pca["T306S-R378K_2"]
d1_pc1, d2_pc1 = t1m - wt1_pc1, t2m - wt2_pc1
davg_pc1 = ((t1m + t2m) / 2) - ((wt1_pc1 + wt2_pc1) / 2)
agree_pc1 = (d1_pc1 > 0) == (d2_pc1 > 0)
robust_pc1 = agree_pc1 and abs(davg_pc1) > nf_pc1
d1_f, d2_f = t1f - wt1_frac, t2f - wt2_frac
davg_f = ((t1f + t2f) / 2) - ((wt1_frac + wt2_frac) / 2)
agree_f = (d1_f > 0) == (d2_f > 0)
robust_f = agree_f and abs(davg_f) > nf_frac
print(f"{'T306S-R378K':16s} {davg_pc1:10.4f} {str(agree_pc1):>7s} {str(robust_pc1):>8s}  "
      f"{davg_f:10.4f} {str(agree_f):>7s} {str(robust_f):>8s}")

# ---------------------------------------------------------------------
# DCCM
# ---------------------------------------------------------------------
def load_correlation_matrix(sys):
    """Reads covar_<sys>.dat (3*462 x 3*462 ascii covariance), reduces to
    a 462x462 correlation matrix."""
    path = f"{sys}/covar_{sys}.dat"
    flat = np.loadtxt(path)
    n3 = 3 * N_RES
    cov = flat.reshape(n3, n3)
    # 3x3 block traces: sum of the 3 diagonal-ish terms per block
    cov3 = cov.reshape(N_RES, 3, N_RES, 3)
    block_trace = np.einsum("iaja->ij", cov3)  # (N_RES, N_RES)
    diag = np.diag(block_trace).copy()
    diag[diag <= 0] = np.nan  # guard against numerical noise
    denom = np.sqrt(np.outer(diag, diag))
    corr = block_trace / denom
    return corr

print("\n" + "=" * 70)
print("DCCM -- global (mean |correlation| across all off-diagonal residue pairs)")
print("=" * 70)

print("Computing correlation matrices for all 24 systems (this reads/reshapes a "
      "~1.9M-element matrix per system, may take a bit)...")
dccm = {}
mask_offdiag = ~np.eye(N_RES, dtype=bool)
for sys in SYSTEMS_ALL:
    corr = load_correlation_matrix(sys)
    global_mean_abs = np.nanmean(np.abs(corr[mask_offdiag]))
    dccm[sys] = (corr, global_mean_abs)
    print(f"  {sys}: global mean|corr|={global_mean_abs:.4f}")

wt1_corr, wt1_g = dccm["WT"]
wt2_corr, wt2_g = dccm["WT_2"]
nf_g = abs(wt1_g - wt2_g)
print(f"\nWT rep1={wt1_g:.4f}, WT rep2={wt2_g:.4f}, noise floor={nf_g:.4f}\n")

print(f"{'Allele':16s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, site) in alleles.items():
    m1, m2 = dccm[d1][1], dccm[d2][1]
    delta1, delta2 = m1 - wt1_g, m2 - wt2_g
    d_avg = ((m1 + m2) / 2) - ((wt1_g + wt2_g) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > nf_g
    print(f"{name:16s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}")

t1g, t2g = dccm["T306S-R378K"][1], dccm["T306S-R378K_2"][1]
delta1, delta2 = t1g - wt1_g, t2g - wt2_g
d_avg = ((t1g + t2g) / 2) - ((wt1_g + wt2_g) / 2)
agree = (delta1 > 0) == (delta2 > 0)
robust = agree and abs(d_avg) > nf_g
print(f"{'T306S-R378K':16s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}")

def site_long_range_coupling(corr, resid, exclude_window=3):
    idx = ca_index(resid)
    others = np.ones(N_RES, dtype=bool)
    lo, hi = max(idx - exclude_window, 0), min(idx + exclude_window + 1, N_RES)
    others[lo:hi] = False
    row = corr[idx, others]
    return np.nanmean(np.abs(row))

print("\n" + "=" * 70)
print("DCCM -- local (mean |correlation| of mutation site vs all residues >3 away)")
print("=" * 70)
print(f"{'Allele':16s} {'d_avg':>8s} {'rep1 d':>8s} {'rep2 d':>8s} {'agree?':>7s} {'ROBUST?':>8s}")
for name, (d1, d2, site) in alleles.items():
    wt1_v = site_long_range_coupling(wt1_corr, site)
    wt2_v = site_long_range_coupling(wt2_corr, site)
    nf_local = abs(wt1_v - wt2_v)
    m1_v = site_long_range_coupling(dccm[d1][0], site)
    m2_v = site_long_range_coupling(dccm[d2][0], site)
    delta1, delta2 = m1_v - wt1_v, m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > nf_local
    print(f"{name:16s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}  (noise_floor={nf_local:.4f})")

for site, label in [(278, "T306"), (350, "R378")]:
    wt1_v = site_long_range_coupling(wt1_corr, site)
    wt2_v = site_long_range_coupling(wt2_corr, site)
    nf_local = abs(wt1_v - wt2_v)
    m1_v = site_long_range_coupling(dccm["T306S-R378K"][0], site)
    m2_v = site_long_range_coupling(dccm["T306S-R378K_2"][0], site)
    delta1, delta2 = m1_v - wt1_v, m2_v - wt2_v
    d_avg = ((m1_v + m2_v) / 2) - ((wt1_v + wt2_v) / 2)
    agree = (delta1 > 0) == (delta2 > 0)
    robust = agree and abs(d_avg) > nf_local
    print(f"T306S-R378K({label}):{'':2s} {d_avg:8.4f} {delta1:8.4f} {delta2:8.4f} {str(agree):>7s} {str(robust):>8s}  (noise_floor={nf_local:.4f})")
