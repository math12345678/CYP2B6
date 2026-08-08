#!/usr/bin/env python
"""Extract aligned PDB snapshots for CAVER 3.0 from each system's md_pocket.xtc.

Reads each system's protein+heme trajectory (601 frames), RMSD-fits every frame
to the reference PDB on protein backbone atoms (rotation+translation applied to
the whole system), and writes a subsampled set of PDB files into
caver_analysis/snapshots/<SYS>/snapshots/ so CAVER 3.0 receives identically
numbered structures in a common coordinate frame.

The CAVER starting point (substrate-cavity void above the heme) is computed as
the WT cavity centroid translated by each system's heme-Fe offset from WT, and
written to caver_analysis/starting_points.tsv.

Usage: source activate cyp2b6; python prep_caver_snapshots.py [SYS...]
"""
import os
import sys
import MDAnalysis as mda
from MDAnalysis.analysis import align

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_ROOT = os.path.dirname(os.path.abspath(__file__))

SYSTEMS = ["WT", "WT_2", "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
           "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2",
           "R140Q", "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2",
           "S259R", "S259R_2", "T306S-R378K", "T306S-R378K_2"]

# Subsampling: every Nth frame of the 601-frame md_pocket.xtc (0.5 ns cadence).
FRAME_STEP = 5   # -> 121 snapshots (~2.5 ns cadence)
FRAME_START = 0  # 0-indexed

BACKBONE = "name CA C N O"
FE_SEL = "resname FE1 and name FE"

# WT substrate-cavity centroid (verified void, clearance ~4.2 A at probe 0.9).
WT_FE = [42.42, 49.12, 40.05]
WT_CAVITY = [48.27, 43.87, 42.05]


def main():
    targets = sys.argv[1:] or SYSTEMS
    start_points = {}
    for sysname in targets:
        ref_pdb = os.path.join(ROOT, sysname, "md_protein_ref.pdb")
        xtc = os.path.join(ROOT, sysname, "md_pocket.xtc")
        if not (os.path.isfile(ref_pdb) and os.path.isfile(xtc)):
            print(f"SKIP {sysname}: missing ref pdb or xtc")
            continue
        out_dir = os.path.join(OUT_ROOT, "snapshots", sysname, "snapshots")
        os.makedirs(out_dir, exist_ok=True)

        ref = mda.Universe(ref_pdb)
        ref_bb = ref.select_atoms(BACKBONE)
        ref_coord = ref_bb.positions

        u = mda.Universe(ref_pdb, xtc)
        n = u.trajectory.n_frames
        frames = list(range(FRAME_START, n, FRAME_STEP))
        bb = u.select_atoms(BACKBONE)

        fe = u.select_atoms(FE_SEL)
        if len(fe) == 0:
            print(f"WARN {sysname}: no FE1 atom found, using WT Fe for offset")
            fe_pos = WT_FE
        else:
            u.trajectory[0]
            fe_pos = fe.positions[0]

        for i in frames:
            u.trajectory[i]
            align.alignto(u, ref, select=BACKBONE, weights="mass")
            fname = os.path.join(out_dir, f"snapshot_{i+1:04d}.pdb")
            u.atoms.write(fname)

        shift = fe_pos - WT_FE
        cavity = list(WT_CAVITY + shift)
        start_points[sysname] = cavity
        print(f"OK {sysname}: {len(frames)} snapshots; cavity start "
              f"{np_round(cavity)}")

    with open(os.path.join(OUT_ROOT, "starting_points.tsv"), "w") as f:
        f.write("system\tx\ty\tz\n")
        for s in SYSTEMS:
            if s in start_points:
                c = start_points[s]
                f.write(f"{s}\t{c[0]:.3f}\t{c[1]:.3f}\t{c[2]:.3f}\n")


def np_round(v):
    import numpy as np
    return [round(float(x), 3) for x in v]


if __name__ == "__main__":
    main()
