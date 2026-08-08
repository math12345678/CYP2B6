#!/usr/bin/env python
"""
Select the active-site pocket from fpocket output as the pocket whose alpha-
sphere centroid lies closest to the heme iron (FE1, GROMACS resid 408).

Rationale (grounded in the CYP2B6 reference structure):
  - CYP2B6 is a heme-containing P450; its active site is the pocket directly
    above the heme iron (distal pocket), where substrates bind and oxygenate.
  - fpocket returns 20-30 candidate pockets per structure; the heme-distal
    pocket is unambiguously identified by proximity to the Fe atom.
  - Validation (WT): the automatically selected pocket is lined by residues
    297, 301, 302, 363, 477, 114, 115 (true numbering) -- exactly the
    canonical CYP2B6 active-site residues from site-directed mutagenesis
    literature (Phe297, Glu301, Thr302, Val363, Val477, Ile114, Phe115).

Usage:
  select_active_pocket.py <SYS>

Writes <SYS>/selected_pocket_<SYS>.pdb (copy of the winning pocket's atom
PDB) and prints the selection summary. Requires <SYS>/md_protein_ref.pdb
and <SYS>/md_protein_ref_out/ (fpocket output from run_one_mdpocket.sh).

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.
"""
import glob
import os
import shutil
import sys

import numpy as np


def fe_iron(pdb_path):
    """Coordinates of the heme iron (FE1)."""
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "FE1":
            return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return None


def pocket_centroid(pqr_path):
    """Centroid of the pocket's alpha-sphere centers (vertices)."""
    xs, ys, zs = [], [], []
    for line in open(pqr_path):
        if line.startswith(("ATOM", "HETATM")):
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
    return np.array([np.mean(xs), np.mean(ys), np.mean(zs)])


def pocket_residues(pdb_path):
    return sorted({int(line[22:26]) for line in open(pdb_path) if line.startswith(("ATOM", "HETATM"))})


def main():
    sysdir = sys.argv[1]
    ref = os.path.join(sysdir, "md_protein_ref.pdb")
    outdir = os.path.join(sysdir, "md_protein_ref_out")
    pockets_dir = os.path.join(outdir, "pockets")

    fe = fe_iron(ref)
    if fe is None:
        print(f"ERROR: no FE1 heme iron found in {ref}")
        sys.exit(1)

    pqrs = sorted(glob.glob(os.path.join(pockets_dir, "pocket*_vert.pqr")))
    if not pqrs:
        print(f"ERROR: no fpocket pockets found in {pockets_dir}")
        sys.exit(1)

    # Pick the pocket whose alpha-sphere centroid is nearest the heme iron.
    best_n, best_d, best_c = None, np.inf, None
    for pqr in pqrs:
        n = int(os.path.basename(pqr).split("_")[0].replace("pocket", ""))
        c = pocket_centroid(pqr)
        d = np.linalg.norm(c - fe)
        if d < best_d:
            best_n, best_d, best_c = n, d, c

    atm_pdb = os.path.join(pockets_dir, f"pocket{best_n}_atm.pdb")
    out_pdb = os.path.join(sysdir, f"selected_pocket_{os.path.basename(sysdir)}.pdb")
    shutil.copy(atm_pdb, out_pdb)

    res = pocket_residues(atm_pdb)
    print(f"{sysdir}: selected pocket {best_n} (alpha-sphere centroid {best_d:.2f} A from heme Fe)")
    print(f"  wrote {out_pdb}")
    print(f"  lining residues (GROMACS numbering): {res}")
    print(f"  true numbering (+28): {sorted(r + 28 for r in res)}")


if __name__ == "__main__":
    main()
