"""
Batch driver for hbond_analysis.py: runs it across every system that has a
ready md_protein_ref.gro (produced by fix_topology.sh), skipping any that are
missing it or that already have output CSVs (so this is safe to re-run after
an interruption). Prints a running summary and a final table.

Run from ~/Desktop/Research/Research_Projects/CYP2B6, inside the cyp2b6
conda env. Expect roughly 3 minutes per system (~20 systems -> ~1 hour).
"""
import os
import subprocess
import sys

SYSTEMS = [
    "WT", "WT_2", "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
    "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2", "R140Q",
    "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2", "S259R", "S259R_2",
    "T306S-R378K", "T306S-R378K_2",
]

results = []
for sys_name in SYSTEMS:
    gro = f"{sys_name}/md_protein_ref.gro"
    pairs_csv = f"{sys_name}/hbond_pairs_{sys_name}.csv"

    if not os.path.exists(gro):
        print(f"SKIP {sys_name}: no md_protein_ref.gro (run fix_topology.sh first)")
        results.append((sys_name, "skipped-no-topology"))
        continue
    if os.path.exists(pairs_csv):
        print(f"SKIP {sys_name}: already has {pairs_csv}")
        results.append((sys_name, "already-done"))
        continue

    print(f"=== Running {sys_name} ===")
    ret = subprocess.run([sys.executable, "hbond_analysis.py", sys_name])
    if ret.returncode != 0:
        print(f"FAILED {sys_name} (exit code {ret.returncode})")
        results.append((sys_name, "failed"))
    else:
        results.append((sys_name, "done"))

print("\n=== Summary ===")
for name, status in results:
    print(f"{name:15s} {status}")
