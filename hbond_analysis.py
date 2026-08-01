"""
Hydrogen bond frequency analysis using MDAnalysis (replaces gmx hbond, which
could not identify any donors/acceptors on this topology -- see README for
the troubleshooting history).

Donors/acceptors are auto-guessed from atom names (not GROMACS element
metadata, which this system's .tpr does not appear to carry in a way the
newer selection-based gmx hbond tool can use). Includes protein + heme
(CM1/HM1/FE1 residues), since donor/acceptor guessing is applied to the whole
loaded structure, not just standard amino acids.

Usage: python3 hbond_analysis.py <SYSTEM_DIR>
  e.g. python3 hbond_analysis.py WT

Requires md_protein_ref.gro and md_noWAT.xtc to already exist in SYSTEM_DIR
(run fix_topology.sh first).

Outputs (written into SYSTEM_DIR):
  hbond_count_<SYSTEM>.csv   -- number of H-bonds per frame
  hbond_pairs_<SYSTEM>.csv   -- unique donor-acceptor pairs, occurrence count,
                                and frequency (fraction of frames present)
"""
import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

if len(sys.argv) != 2:
    print("Usage: python3 hbond_analysis.py <SYSTEM_DIR>")
    sys.exit(1)

system_dir = sys.argv[1]
system_name = os.path.basename(system_dir.rstrip("/"))
gro = os.path.join(system_dir, "md_protein_ref.gro")
xtc = os.path.join(system_dir, "md_noWAT.xtc")

if not os.path.exists(gro) or not os.path.exists(xtc):
    print(f"ERROR: missing {gro} or {xtc}. Run fix_topology.sh first.")
    sys.exit(1)

print(f"Loading {system_name} ({gro}, {xtc}) ...")
u = mda.Universe(gro, xtc)
print(f"  {len(u.atoms)} atoms, {len(u.trajectory)} frames")

# The .gro topology carries no partial charges, so MDAnalysis's default
# acceptor-guessing (which filters by charge < -0.5) cannot run. Guess
# elements and bonds from atom names/distances instead, and specify
# acceptors explicitly by element (N, O) rather than by charge. Donors are
# still guessed normally: guess_hydrogens() finds H atoms by name, and
# guess_donors() pairs them to their bonded heavy atom via the guessed bonds
# -- that path does not require charges.
if not hasattr(u.atoms, "elements"):
    u.guess_TopologyAttrs(to_guess=["elements"])
if not hasattr(u.atoms, "bonds") or len(u.atoms.bonds) == 0:
    # Default MDAnalysis vdW radius table has no entry for Fe (the heme
    # cofactor's iron). Supply one manually so bond-guessing doesn't fail;
    # 1.5 Angstrom is a generous vdW-scale radius for coordinated heme Fe.
    u.atoms.guess_bonds(vdwradii={"FE": 1.5})

# This topology carries no partial charges, and MDAnalysis's guess_hydrogens/
# guess_donors/guess_acceptors all default to charge-based heuristics. Skip
# guessing entirely and specify all three selections explicitly by element:
# hydrogens are H atoms; donors/acceptors are the N/O atoms they can be
# bonded to (pairing is resolved via the guessed bonds, not by these
# selections alone).
hbonds = HydrogenBondAnalysis(
    universe=u,
    hydrogens_sel="element H",
    donors_sel="element O or element N",
    acceptors_sel="element O or element N",
    d_a_cutoff=3.5,        # donor-acceptor distance cutoff (Angstrom)
    d_h_a_angle_cutoff=150,  # donor-H-acceptor angle cutoff (degrees)
    update_selections=False,
)

print(f"  donors_sel: {hbonds.donors_sel!r}")
print(f"  hydrogens_sel: {hbonds.hydrogens_sel!r}")
print(f"  acceptors_sel: {hbonds.acceptors_sel!r}")
print(f"  n_bonds guessed: {len(u.atoms.bonds)}")

hbonds.run(verbose=True)

results = hbonds.results.hbonds  # columns: frame, donor_ix, hydrogen_ix, acceptor_ix, distance, angle
n_frames = len(u.trajectory)

if results.shape[0] == 0:
    print(f"WARNING: 0 hydrogen bonds detected in {system_name}. Check donor/acceptor selections above.")
    sys.exit(1)

# --- per-frame count ---
frames = results[:, 0].astype(int)
counts_per_frame = np.bincount(frames, minlength=n_frames)
count_path = f"{system_dir}/hbond_count_{system_name}.csv"
np.savetxt(count_path, np.column_stack([np.arange(n_frames), counts_per_frame]),
           header="frame,hbond_count", delimiter=",", comments="")
print(f"Saved {count_path}")

# --- unique donor-acceptor pair frequencies ---
donor_ix = results[:, 1].astype(int)
acceptor_ix = results[:, 3].astype(int)
pairs = list(zip(donor_ix, acceptor_ix))
unique_pairs, counts = np.unique(pairs, axis=0, return_counts=True)

pair_rows = []
for (d_ix, a_ix), c in zip(unique_pairs, counts):
    d_atom = u.atoms[d_ix]
    a_atom = u.atoms[a_ix]
    freq = c / n_frames
    pair_rows.append((
        d_atom.resname, d_atom.resid, d_atom.name,
        a_atom.resname, a_atom.resid, a_atom.name,
        c, freq,
    ))

pair_rows.sort(key=lambda r: -r[-1])

pairs_path = f"{system_dir}/hbond_pairs_{system_name}.csv"
with open(pairs_path, "w") as f:
    f.write("donor_resname,donor_resid,donor_atom,acceptor_resname,acceptor_resid,acceptor_atom,count,frequency\n")
    for row in pair_rows:
        f.write(",".join(str(x) for x in row) + "\n")
print(f"Saved {pairs_path}")

print(f"\n{system_name}: {len(pair_rows)} unique donor-acceptor pairs, "
      f"mean {counts_per_frame.mean():.1f} H-bonds/frame")
