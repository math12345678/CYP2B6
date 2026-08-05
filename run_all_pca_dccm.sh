#!/bin/bash
# PCA (essential dynamics) and DCCM (dynamic cross-correlation matrix),
# following the handover document's "Recommended Analyses" list (Part 7).
#
# Both are derived from the same covariance-matrix machinery (gmx covar),
# but computed on different atom selections for a documented reason:
#   - PCA uses the Backbone group (index group 4), matching this project's
#     established RMSD/RMSF convention -- essential dynamics of the whole
#     chain backbone.
#   - DCCM uses the C-alpha group (index group 3) specifically, since one
#     CA atom per residue makes the resulting covariance matrix directly
#     usable as a residue-residue matrix without further reduction (the
#     same "pick the group that matches the question" reasoning already
#     documented for Rg/SASA using the Protein group instead of Backbone).
#
# gmx covar performs a least-squares fit using the first prompted group,
# then computes the covariance matrix for the second prompted group -- same
# group is used for both here in each case.
#
# Outputs:
#   PCA:  eigenval_<SYS>.xvg (eigenvalues/variance spectrum, small, tracked)
#         eigenvec_<SYS>.trr (eigenvectors, binary, large, gitignored)
#         proj_<SYS>.xvg (trajectory projected onto PC1/PC2, small, tracked)
#   DCCM: covar_<SYS>.dat (-ascii full 3N x 3N covariance matrix for the
#         C-alpha selection, large [~15-20MB/system], gitignored -- reduced
#         to a small per-residue-pair correlation summary by
#         pca_dccm_summary.py)
#
# Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.

set -e

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

for sys in $SYSTEMS; do
    dir="$sys"
    if [ ! -f "$dir/md.tpr" ] || [ ! -f "$dir/md_noWAT.xtc" ] || [ ! -f "$dir/index.ndx" ]; then
        echo "SKIP: $dir missing md.tpr, md_noWAT.xtc, or index.ndx"
        continue
    fi

    if [ ! -f "$dir/eigenval_${sys}.xvg" ]; then
        echo "=== $sys: PCA (gmx covar + anaeig, Backbone) ==="
        printf "4\n4\n" | gmx covar -s "$dir/md.tpr" -n "$dir/index.ndx" -f "$dir/md_noWAT.xtc" \
            -o "$dir/eigenval_${sys}.xvg" -v "$dir/eigenvec_${sys}.trr" -tu ns
        # anaeig prompts twice: once for the fit group used in covar, once
        # for the group matching the eigenvectors -- both Backbone here.
        printf "4\n4\n" | gmx anaeig -s "$dir/md.tpr" -n "$dir/index.ndx" -f "$dir/md_noWAT.xtc" \
            -v "$dir/eigenvec_${sys}.trr" -eig "$dir/eigenval_${sys}.xvg" \
            -first 1 -last 2 -proj "$dir/proj_${sys}.xvg" -tu ns
    else
        echo "SKIP: $dir already has PCA output"
    fi

    if [ ! -f "$dir/covar_${sys}.dat" ]; then
        echo "=== $sys: DCCM (gmx covar -ascii, C-alpha) ==="
        # Fit on Backbone (group 4, matching this project's RMSD/RMSF
        # convention) to remove rigid-body translation/rotation first --
        # DCCM must be computed on internal fluctuations only, not raw
        # coordinates, or overall tumbling would swamp the real residue-
        # residue correlation signal. Analysis group is C-alpha (group 3).
        # -v given an explicit name so this run's eigenvectors don't
        # collide with (or get overwritten by) the next system's, or with
        # the PCA run's own eigenvec_<sys>.trr above.
        printf "4\n3\n" | gmx covar -s "$dir/md.tpr" -n "$dir/index.ndx" -f "$dir/md_noWAT.xtc" \
            -o "$dir/eigenval_ca_${sys}.xvg" -ascii "$dir/covar_${sys}.dat" \
            -v "$dir/eigenvec_ca_${sys}.trr" -tu ns
    else
        echo "SKIP: $dir already has DCCM output"
    fi
done

echo "Done."
