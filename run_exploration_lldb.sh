#!/bin/bash
# Exploration-pass fixup for run_all_mdpocket.sh.
#
# Why this exists: the conda-forge mdpocket binary crashes at startup in
# get_mdpocket_args/my_malloc (Trace/BPT trap, EXC_BREAKPOINT inside
# libsystem_malloc's mfm_alloc) whenever the process is launched with a
# random ASLR layout -- the exploration invocation (no --selected_pocket)
# fails ~100% of direct launches. Running the identical command under lldb
# (which disables ASLR by default) is much more reliable: the remaining
# crashes are flaky (rpdb_read null-FILE* on a failed fopen, or mfm_alloc),
# and are absorbed by the retry loop here.
#
# IMPORTANT: never run two mdpocket exploration processes in the SAME
# system directory -- they clobber each other's output files and the
# failures look like rpdb_read crashes. This wrapper is per-system serial;
# the batch driver runs systems in parallel only across DIFFERENT dirs.
#
# Usage (from the repo root): run_exploration_lldb.sh <SYS>
set -u
SYS="$1"
cd "$SYS" || { echo "FATAL $SYS: cannot cd into $SYS"; exit 1; }

# Timeout per attempt: exploration over 601 snapshots takes ~25-60 min under
# parallel load; the original 1500 s (25 min) killed every run mid-grid.
TIMEOUT_SEC="${TIMEOUT_SEC:-5400}"
ATTEMPTS="${ATTEMPTS:-4}"

grid_ok() { [ -s mdpout_dens_grid.dx ]; }

rm -f mdpout_dens_grid.dx mdpout_freq_grid.dx mdpout_freq_iso_0_5.pdb \
      mdpout_all_atom_pdensities.pdb mdpout_dens_iso_8.pdb mdpout.log

for attempt in $(seq 1 "$ATTEMPTS"); do
    echo "== $SYS: mdpocket exploration via lldb (ASLR off), attempt $attempt (timeout ${TIMEOUT_SEC}s)"
    timeout "$TIMEOUT_SEC" lldb -b -o run -o quit -- /opt/anaconda3/envs/cyp2b6/bin/mdpocket \
        --trajectory_file md_pocket.xtc --trajectory_format xtc \
        -f md_protein_ref.pdb > mdpout.log 2>&1
    rc=$?
    # NOTE: do NOT pkill here -- other parallel systems run the same binary.
    if grid_ok; then
        echo "DONE $SYS (attempt $attempt, lldb rc=$rc)"
        exit 0
    fi
    echo "  attempt $attempt failed (lldb rc=$rc, grid $(stat -f%z mdpout_dens_grid.dx 2>/dev/null || echo missing) bytes)"
    tail -3 mdpout.log
done

echo "WARN $SYS: exploration still failed under lldb after $ATTEMPTS attempts"
exit 0
