#!/bin/bash
# Requires AmberTools (cpptraj). If not installed:
#   conda create -n ambertools -c conda-forge ambertools -y
#   conda activate ambertools
#
# Conformational clustering via AmberTools cpptraj, matching the approach used
# in Shaylyn Govender's predecessor MSc thesis on this exact system (Chapter 4:
# "Analysis through Clustering and DSSP"): hierarchical agglomerative
# clustering using RMSD as the distance metric, on the RUBi Jango server.
#
# Deviation from her exact protocol, documented rather than silently assumed:
# her thesis also used the heme Fe and "L-helix" as an additional distance
# metric, but doesn't give the L-helix residue range in a form that could be
# extracted from the text. This run uses RMSD on the Backbone atom mask only
# (matching this project's established Backbone convention for RMSD/RMSF),
# not the heme Fe + L-helix combination. Cluster *counts* and *fractions*
# should still be broadly comparable to her Table 4.1, but not identical.
#
# cpptraj cannot parse GROMACS .tpr topologies directly (confirmed: "Could
# not determine format of topology" on a real test run). Uses
# md_protein_ref.pdb instead -- the protein+heme-only PDB already generated
# for the DRN analysis (via gmx editconf, see convert_to_pdb.sh), which
# matches md_noWAT.xtc's atom composition exactly.
#
# Frame stride: hierarchical agglomerative clustering needs a full pairwise
# RMSD matrix (O(n^2) memory/time), which is intractable at the full 30001
# frames/system. Using every 10th frame (~3000 frames, same ~1 ns-scale
# subsampling reasoning as the DRN analysis's --step 100) keeps this
# tractable while still covering the whole 300 ns trajectory.
#
# Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6, inside whichever
# environment has cpptraj on PATH (AmberTools, not the cyp2b6/mdmtaskweb
# conda envs).

set -e

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

for sys in $SYSTEMS; do
    dir="$sys"
    if [ ! -f "$dir/md_protein_ref.pdb" ] || [ ! -f "$dir/md_noWAT.xtc" ]; then
        echo "SKIP: $dir missing md_protein_ref.pdb or md_noWAT.xtc"
        continue
    fi
    if [ -f "$dir/cluster_${sys}_summary.dat" ]; then
        echo "SKIP: $dir already clustered"
        continue
    fi

    echo "=== $sys: clustering ==="
    cat > "$dir/cluster_${sys}.in" << EOF
parm $dir/md_protein_ref.pdb
trajin $dir/md_noWAT.xtc 1 last 10
rms myrmsd @C,CA,N,O first out $dir/cluster_${sys}_rmsd.dat
cluster hieragglo clusters 3 \
    rms @C,CA,N,O \
    out $dir/cluster_${sys}_assignments.dat \
    summary $dir/cluster_${sys}_summary.dat \
    info $dir/cluster_${sys}_info.dat \
    repout $dir/cluster_${sys}_rep repfmt pdb
go
EOF
    cpptraj -i "$dir/cluster_${sys}.in"
done

echo "Done."
