#!/bin/bash
# Generates a protein+heme-only topology (md_protein.tpr, md_protein_ref.gro)
# for every system, matching the atom count of the already-stripped
# md_noWAT.xtc trajectories. Needed because md.tpr is the full system
# (protein+water+ions) while md_noWAT.xtc is protein+heme only (7587 atoms) --
# MDAnalysis requires the topology and trajectory atom counts to match exactly.
#
# Group 26 (Protein_CM1_HM1_FE1) in the original md.tpr = protein + heme
# cofactor (CM1/HM1/FE1), which is the group md_noWAT.xtc was actually
# generated from. After convert-tpr subsets to that group, its own internal
# group numbering resets, and "System" (group 0) in the new subset tpr
# becomes the full protein+heme set -- confirmed working on WT (7587 == 7587,
# heme residues present).
#
# Usage: ./fix_topology.sh
# Run from ~/Desktop/Research/Research_Projects/CYP2B6

set -e

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

for sys in $SYSTEMS; do
    dir="$sys"
    if [ ! -d "$dir" ]; then
        echo "SKIP: $dir not found"
        continue
    fi
    if [ ! -f "$dir/md.tpr" ] || [ ! -f "$dir/md_noWAT.xtc" ]; then
        echo "SKIP: $dir missing md.tpr or md_noWAT.xtc"
        continue
    fi
    if [ -f "$dir/md_protein_ref.gro" ]; then
        echo "SKIP: $dir already has md_protein_ref.gro"
        continue
    fi

    echo "=== Processing $dir ==="
    (
        cd "$dir"
        echo 26 | gmx convert-tpr -s md.tpr -n index.ndx -o md_protein.tpr
        echo 0 | gmx trjconv -s md_protein.tpr -f md_noWAT.xtc -o md_protein_ref.gro -dump 0
    )
    echo "=== Done $dir ==="
done

echo "All systems processed. Verify atom counts with check_topology.py before running H-bond analysis."
