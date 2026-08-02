"""
Batch driver for MD-TASK's calc_network.py: computes dynamic residue network
(DRN) centrality metrics (betweenness, closeness, eigenvector centrality --
the three specified in the project timeline) for every system, in priority
order (WT first as the universal baseline, then P428T/I328T -- the two
alleles with convergent RMSD/RMSF/H-bond evidence -- then T306S-R378K/S259R
-- the two still-inconclusive uncertain-function alleles -- then the rest).

--step 100 samples every ~1 ns across the 300 ns trajectory (300 frames),
chosen from a quick timing test (~37s for 30 frames / 3 metrics on P428T,
i.e. roughly 6 min/system at this step size -- ~2.3 hours for all 22
systems).

Safe to re-run after an interruption: skips any system whose output CSV
already exists.

Run from ~/Desktop/Research/Research_Projects/CYP2B6, inside the mdmtaskweb
conda env, after convert_to_pdb.sh.
"""
import os
import subprocess
import sys

PRIORITY_ORDER = [
    "WT", "WT_2",
    "P428T", "P428T_2", "I328T", "I328T_2",
    "T306S-R378K", "T306S-R378K_2", "S259R", "S259R_2",
    "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
    "I391N", "I391N_2", "K262R", "K262R_2", "R140Q", "R140Q_2",
    "R487C", "R487C_2",
]

STEP = "100"
# calc_network.py writes all outputs (per-frame .dat files and the final
# _mean.csv/_mean_*.cif) into the current working directory, using only the
# trajectory's bare filename -- it ignores any directory prefix passed in
# the path. So each system MUST be run with that system's own folder as the
# working directory, or every system's outputs collide in one shared folder
# (this happened on the first real attempt, at the CYP2B6 top level).
CALC_NETWORK_ABS = os.path.abspath("MD-TASK/src/calc_network.py")

results = []
for sys_name in PRIORITY_ORDER:
    sys_dir = os.path.abspath(sys_name)
    pdb = "md_protein_ref.pdb"
    xtc = "md_noWAT.xtc"
    out_csv = os.path.join(sys_dir, "md_noWAT_mean.csv")

    if not os.path.exists(os.path.join(sys_dir, pdb)) or not os.path.exists(os.path.join(sys_dir, xtc)):
        print(f"SKIP {sys_name}: missing {pdb} or {xtc} (run convert_to_pdb.sh first)")
        results.append((sys_name, "skipped-no-input"))
        continue
    if os.path.exists(out_csv):
        print(f"SKIP {sys_name}: already has md_noWAT_mean.csv")
        results.append((sys_name, "already-done"))
        continue

    print(f"=== Running {sys_name} ===")
    ret = subprocess.run([
        sys.executable, CALC_NETWORK_ABS,
        "--topology", pdb,
        "--step", STEP,
        "--calc-BC", "--calc-CC", "--calc-EC",
        xtc,
    ], cwd=sys_dir)
    if ret.returncode != 0:
        print(f"FAILED {sys_name} (exit code {ret.returncode})")
        results.append((sys_name, "failed"))
    else:
        results.append((sys_name, "done"))

print("\n=== Summary ===")
for name, status in results:
    print(f"{name:15s} {status}")
