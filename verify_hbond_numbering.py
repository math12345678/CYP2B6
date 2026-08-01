"""
Independent verification of the H-bond pipeline's residue numbering, run
after the fact to check assumptions that were never explicitly confirmed:

1. Does MDAnalysis's `resid` for the mutation-site residue actually match
   what Part 1's RMSD/RMSF analysis established (GROMACS residue = true
   residue - 28), i.e. is G99E's true residue 99 really GROMACS/MDAnalysis
   residue 71, and is it really glutamate (mutated from glycine)?
2. Where do the heme cofactor's residues (CM1/HM1/FE1) land relative to the
   protein's last residue number? This matters because R487C's mutation
   site is GROMACS residue 459, close to the protein's C-terminus (~462
   residues) -- if heme residues were numbered immediately afterward in a
   way that could fall inside a +/-3 window near the end of the chain, the
   site_window_sum() function in hbond_significance_check.py could
   inadvertently include heme atoms in a "mutation site" sum. This checks
   whether that's actually a risk or not.

Run from ~/Desktop/Research/Research_Projects/CYP2B6, inside the cyp2b6 env.
"""
import MDAnalysis as mda

print("=== Check 1: mutation-site residue identity in the H-bond topology ===")
sites_to_check = {
    "G99E": ("G99E/md_protein_ref.gro", 71, "true G99 -> expect GLU/GLH"),
    "M46V": ("M46V/md_protein_ref.gro", 18, "true M46 -> expect VAL (mutated from MET)"),
    "P428T": ("P428T/md_protein_ref.gro", 400, "true P428 -> expect THR (mutated from PRO)"),
    "R487C": ("R487C/md_protein_ref.gro", 459, "true R487 -> expect CYS (mutated from ARG)"),
}
for name, (gro, site, expectation) in sites_to_check.items():
    u = mda.Universe(gro)
    res = u.residues[u.residues.resids == site]
    if len(res) == 0:
        print(f"{name}: NO RESIDUE FOUND at resid {site} -- MISMATCH")
        continue
    resname = res.resnames[0]
    print(f"{name}: resid {site} -> resname {resname}  ({expectation})")

print("\n=== Check 2: protein max resid vs heme resid(s) ===")
for name, gro in [("WT", "WT/md_protein_ref.gro"), ("R487C", "R487C/md_protein_ref.gro")]:
    u = mda.Universe(gro)
    protein = u.select_atoms("protein")
    heme = u.select_atoms("resname CM1 or resname HM1 or resname FE1")
    protein_max_resid = protein.resids.max() if len(protein) else None
    heme_resids = sorted(set(heme.resids)) if len(heme) else []
    print(f"{name}: protein max resid = {protein_max_resid}, heme resid(s) = {heme_resids}")
    if heme_resids and protein_max_resid is not None:
        too_close = [r for r in heme_resids if abs(r - protein_max_resid) <= 6]
        if too_close:
            print(f"  WARNING: heme resid(s) {too_close} within 6 of protein max resid "
                  f"{protein_max_resid} -- a +/-3 window near the C-terminus could pull in heme atoms")
        else:
            print(f"  OK: heme resids are far enough from protein max resid; no window overlap risk")
