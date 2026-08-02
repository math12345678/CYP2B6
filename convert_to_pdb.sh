#!/bin/bash
# Converts each system's md_protein_ref.gro (protein+heme, produced by
# fix_topology.sh for the H-bond pipeline) to a PDB file, since MD-TASK's
# calc_network.py requires a PDB topology rather than .gro/.tpr.
#
# Usage: ./convert_to_pdb.sh
# Run from ~/Desktop/Research/Research_Projects/CYP2B6

set -e

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

for sys in $SYSTEMS; do
    dir="$sys"
    if [ ! -f "$dir/md_protein_ref.gro" ]; then
        echo "SKIP: $dir missing md_protein_ref.gro (run fix_topology.sh first)"
        continue
    fi
    if [ -f "$dir/md_protein_ref.pdb" ]; then
        echo "SKIP: $dir already has md_protein_ref.pdb"
        continue
    fi
    echo "=== Converting $dir ==="
    gmx editconf -f "$dir/md_protein_ref.gro" -o "$dir/md_protein_ref.pdb"
done

echo "Done."
