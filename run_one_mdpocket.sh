#!/bin/bash
# Per-system runner for the Part 8 binding-pocket (active-site) analysis,
# following the MDpocket two-round protocol of the reference paper
# (Rehema et al., JMB 2025, CYP3A4):
#
#   Round 1 -- pocket exploration: MDpocket over the full trajectory,
#              producing pocket-density/frequency grids that show where
#              transient pockets open and close over the simulation.
#   Round 2 -- pocket characterization: MDpocket restricted to the
#              active-site pocket (selected by select_active_pocket.py as
#              the fpocket pocket closest to the heme iron), producing the
#              per-frame active-site pocket volume time series
#              (mdpock_<SYS>_descriptors.txt, "pock_volume" column) -- the
#              analog of the reference paper's Fig 4D active-site volume.
#
# Trajectory handling:
#   - The full 30001-frame water-stripped trajectory (md_noWAT.xtc, 10 ps
#     cadence) is subsampled to every 50th frame with cpptraj -> 601 frames
#     at 0.5 ns cadence. This matches the paper's sampling density (1001
#     snapshots every 0.5 ns over a 500 ns trajectory) while keeping each
#     system's MDpocket run at a few minutes instead of hours.
#   - The protein-only reference (md_protein_ref.pdb) is used as the
#     structure/topology for both cpptraj and MDpocket; it contains the
#     heme (CM1/HM1/FE1, resid 408), so the selected active-site pocket is
#     properly flanked by the cofactor.
#
# Usage: run_one_mdpocket.sh <SYS>
# (dispatched in parallel by run_all_mdpocket.sh; run from repo root)
set -e

SYS="$1"
CPPTRAJ=/opt/anaconda3/envs/ambertools/bin/cpptraj
FPOCKET=/opt/anaconda3/envs/cyp2b6/bin/fpocket
MDPOCKET=/opt/anaconda3/envs/cyp2b6/bin/mdpocket
PY=/opt/anaconda3/envs/cyp2b6/bin/python

# Pin OpenMP to 1 thread: fpocket/mdpocket link the OpenMP runtime, and 8
# concurrent instances each spawning a thread team oversubscribed the 14-core
# machine and crashed some runs (Trace/BPT trap, ~snapshot 41). Per-frame
# pocket volume is computed single-threaded anyway, so this costs nothing.
export OMP_NUM_THREADS=1

if [ ! -f "$SYS/md_noWAT.xtc" ] || [ ! -f "$SYS/md_protein_ref.pdb" ]; then
    echo "SKIP $SYS: missing md_noWAT.xtc or md_protein_ref.pdb"
    exit 0
fi

cd "$SYS" || exit 1

# 1. Subsampled trajectory (601 frames, 0.5 ns cadence)
if [ ! -f md_pocket.xtc ]; then
    cat > pocket_traj.in <<EOF
trajin md_noWAT.xtc 1 30001 50
trajout md_pocket.xtc
run
EOF
    echo "== $SYS: cpptraj subsample (stride 50)"
    $CPPTRAJ -p md_protein_ref.pdb -i pocket_traj.in > pocket_traj.log 2>&1
fi

# 2. fpocket pocket detection on the reference structure (fast)
if [ ! -d md_protein_ref_out ]; then
    echo "== $SYS: fpocket"
    $FPOCKET -f md_protein_ref.pdb > fpocket.log 2>&1
fi

# 3. Select the active-site pocket (closest to heme Fe)
if [ ! -f "selected_pocket_${SYS}.pdb" ]; then
    echo "== $SYS: select active-site pocket"
    # Pass the absolute path: the script resolves inputs relative to the
    # repo root, and needs the system name for the output filename.
    $PY ../select_active_pocket.py "$(pwd)"
fi

# 4. Round 2 -- characterization: per-frame active-site volume.
#    Completed descriptors file = header + 601 snapshot rows. Retry if the
#    run crashed or was truncated: the conda-forge mdpocket binary has a
#    flaky startup bug (crashes in get_mdpocket_args/my_malloc with
#    Trace/BPT trap before any analysis, ~30-40% of launches under parallel
#    load -- confirmed identical across every macOS crash report). A run
#    that survives startup always completes.
attempt=0
while [ "$(wc -l < "mdpock_${SYS}_descriptors.txt" 2>/dev/null || echo 0)" -lt 602 ]; do
    attempt=$((attempt + 1))
    if [ "$attempt" -gt 5 ]; then
        echo "FAILED $SYS: mdpocket characterization incomplete after 5 attempts"
        exit 1
    fi
    rm -f "mdpock_${SYS}_descriptors.txt" "mdpock_${SYS}_mdpocket.pdb" "mdpock_${SYS}_mdpocket_atoms.pdb"
    echo "== $SYS: mdpocket characterization (selected pocket), attempt $attempt"
    $MDPOCKET --trajectory_file md_pocket.xtc --trajectory_format xtc \
        --selected_pocket "selected_pocket_${SYS}.pdb" \
        -f md_protein_ref.pdb -o "mdpock_${SYS}" > "mdpock_${SYS}.log" 2>&1 || true
done

# 5. Round 1 -- exploration: whole-protein pocket frequency map (retry on
#    the same startup flakiness; grid written at the end of a successful run)
#    NOTE: check for a NON-EMPTY grid (-s). A 0-byte file is a failed run;
#    MDpocket writes the grid only once the full trajectory is processed.
if [ ! -s mdpout_dens_grid.dx ]; then
    attempt=0
    while [ ! -s mdpout_dens_grid.dx ]; do
        attempt=$((attempt + 1))
        if [ "$attempt" -gt 5 ]; then
            echo "WARN $SYS: mdpocket exploration incomplete after 5 attempts"
            break
        fi
        rm -f mdpout_dens_grid.dx mdpout_freq_grid.dx mdpout_freq_iso_0_5.pdb \
              mdpout_all_atom_pdensities.pdb mdpout_dens_iso_8.pdb
        echo "== $SYS: mdpocket exploration (frequency map), attempt $attempt"
        $MDPOCKET --trajectory_file md_pocket.xtc --trajectory_format xtc \
            -f md_protein_ref.pdb > mdpout.log 2>&1 || true
    done
fi

echo "DONE $SYS"
