#!/bin/bash
# Batch driver for Rg (gmx gyrate) and SASA (gmx sasa) across all 22 systems,
# per the handover doc's "Recommended Analyses" list (Rg, SASA under "Global
# stability", alongside RMSD/RMSF which are already done).
#
# Uses the same command pattern as the existing RMSD/RMSF commands (README.md
# "Commands run" section): -s md.tpr -n index.ndx -f md_noWAT.xtc.
#
# Group selection: RMSD/RMSF used group 4 (Backbone) for consistency across
# systems. Rg and SASA are computed on group 1 (Protein) instead -- Backbone-
# only would ignore side chains, which matters for both overall compactness
# (Rg) and solvent-accessible surface (SASA); using the full protein is the
# standard convention for both. This is a deliberate, documented deviation
# from the Backbone choice used for RMSD/RMSF, not an inconsistency.
#
# Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6, inside whichever
# environment has GROMACS on PATH (not the cyp2b6/mdmtaskweb conda envs used
# for the Python analyses).

set -e

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

for sys in $SYSTEMS; do
    dir="$sys"
    if [ ! -f "$dir/md.tpr" ] || [ ! -f "$dir/md_noWAT.xtc" ]; then
        echo "SKIP: $dir missing md.tpr or md_noWAT.xtc"
        continue
    fi
    if [ -f "$dir/gyrate_${sys}.xvg" ] && [ -f "$dir/sasa_${sys}.xvg" ]; then
        echo "SKIP: $dir already has gyrate/sasa output"
        continue
    fi

    echo "=== $dir: Rg ==="
    echo "1" | gmx gyrate -s "$dir/md.tpr" -n "$dir/index.ndx" -f "$dir/md_noWAT.xtc" \
        -o "$dir/gyrate_${sys}.xvg"

    echo "=== $dir: SASA ==="
    echo "1" | gmx sasa -s "$dir/md.tpr" -n "$dir/index.ndx" -f "$dir/md_noWAT.xtc" \
        -o "$dir/sasa_${sys}.xvg" -or "$dir/sasa_res_${sys}.xvg"
done

echo "Done."
