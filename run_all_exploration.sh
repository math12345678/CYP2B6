#!/bin/bash
# Part 8 -- Round 1 pocket-exploration fixup driver.
#
# The original run_all_mdpocket.sh launched exploration in parallel (P=8)
# with a 1500 s per-attempt timeout; under that load MDpocket exploration
# over 601 snapshots never finished before the timeout, so every system's
# mdpout_dens_grid.dx was left 0 bytes. This driver reruns only the systems
# that still have an empty/missing grid, using the lldb wrapper (ASLR off,
# which dodges the mdpocket startup crash) with a much longer timeout and a
# lower parallelism so each run actually completes.
#
# IMPORTANT: exploration writes fixed output names (mdpout_*), so never run
# two exploration jobs in the SAME system directory. This driver runs systems
# in parallel only across DIFFERENT directories (P=3 to keep each run fast).
#
# Usage (from the repo root): bash run_all_exploration.sh

set -u
P=${P:-3}

SYSTEMS="WT WT_2 G99E G99E_2 K139E K139E_2 M46V M46V_2 I328T I328T_2 I391N I391N_2 K262R K262R_2 R140Q R140Q_2 R487C R487C_2 P428T P428T_2 S259R S259R_2 T306S-R378K T306S-R378K_2"

NEED=""
for S in $SYSTEMS; do
    if [ ! -s "$S/mdpout_dens_grid.dx" ]; then
        NEED="$NEED $S"
    fi
done

echo "Systems needing exploration:$NEED"
echo "$NEED" | tr ' ' '\n' | sed '/^$/d' | xargs -P "$P" -I{} bash run_exploration_lldb.sh {}

echo "Exploration driver done."
