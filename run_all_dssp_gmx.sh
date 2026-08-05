#!/bin/bash
# Replaces the cpptraj-based run_all_dssp.sh, which produced unusable output:
# every residue across all 22 systems showed ~0% helix/sheet content, which
# is impossible for CYP2B6 (a heavily alpha-helical P450). Root cause: this
# is a GROMACS simulation (md.tpr / cyp2b6_GMX.top, no AMBER prmtop exists),
# and md_protein_ref.pdb's bonding was inferred by cpptraj from atom
# distances rather than read from a real topology -- that broke backbone
# peptide-bond connectivity between some residues, so cpptraj's DSSP H-bond
# ladder detection couldn't form proper helix/sheet patterns.
#
# Fix: use GROMACS' own gmx dssp (needs GROMACS 2023+ with the built-in DSSP
# module, or a system dssp/mkdssp binary that gmx dssp can call out to).
# This also exactly matches Shaylyn Govender's thesis method (Chapter 4,
# "Analysis through Clustering and DSSP" used gmx dssp), so results are
# directly comparable to her prior findings, not just internally consistent.
#
# Uses the same "Backbone" group (4) as RMSD/RMSF, since gmx dssp operates
# on the whole protein backbone directly (no group selection prompt needed
# for the assignment itself -- gmx dssp takes -s/-f and computes over the
# full system, we restrict with an index group only if needed).
#
# Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.

set -e

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

for sys in $SYSTEMS; do
    dir="$sys"
    if [ ! -f "$dir/md.tpr" ] || [ ! -f "$dir/md_noWAT.xtc" ]; then
        echo "SKIP: $dir missing md.tpr or md_noWAT.xtc"
        continue
    fi
    if [ -f "$dir/dssp_${sys}.dat" ]; then
        echo "SKIP: $dir already has gmx dssp output"
        continue
    fi

    echo "=== $sys: gmx dssp ==="
    # -tu ns just for readable timestamps; -o writes per-residue, per-frame
    # SS as a compact string matrix (dssp_<sys>.dat).
    echo "1" | gmx dssp -s "$dir/md.tpr" -f "$dir/md_noWAT.xtc" \
        -o "$dir/dssp_${sys}.dat" -tu ns
done

echo "Done."
