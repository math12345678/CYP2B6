#!/bin/bash
# Part 8 binding-pocket (active-site) analysis -- batch driver for all
# systems, following the reference paper's MDpocket protocol (Rehema et al.,
# JMB 2025, CYP3A4): two rounds per system (whole-protein pocket
# exploration + selected active-site pocket characterization) on 601-frame
# subsampled trajectories (0.5 ns cadence, matching the paper's sampling
# density).
#
# Per-system work (see run_one_mdpocket.sh):
#   cpptraj subsample -> fpocket reference pockets -> active-site pocket
#   selection (nearest to heme Fe) -> MDpocket characterization (per-frame
#   active-site volume) -> MDpocket exploration (pocket frequency map).
#
# Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.
#
# Parallelism: MDpocket is single-threaded and the machine has 14 cores, so
# systems are processed 8 at a time by default (~12 min per system, so the
# full 24-system batch takes ~40 min). Set P=1 for strictly sequential
# execution.

set -e
P=${P:-8}

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

echo "$SYSTEMS" | tr ' ' '\n' | xargs -P "$P" -I{} bash run_one_mdpocket.sh {}

echo "All systems done."
