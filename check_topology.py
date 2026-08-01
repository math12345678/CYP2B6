"""
Verifies, for every system, that md_protein_ref.gro and md_noWAT.xtc have
matching atom counts before trusting them for H-bond analysis. fix_topology.sh
assumes group 26 (Protein_CM1_HM1_FE1) is at the same index in every system's
md.tpr -- this is very likely true since all systems were built through the
same pipeline, but this script checks the actual result rather than assuming
the blind loop worked identically everywhere.

Run from ~/Desktop/Research/Research_Projects/CYP2B6 after fix_topology.sh.
"""
import MDAnalysis as mda

SYSTEMS = [
    "WT", "WT_2", "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
    "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2", "R140Q",
    "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2", "S259R", "S259R_2",
    "T306S-R378K", "T306S-R378K_2",
]

def gro_atom_count(path):
    with open(path) as f:
        f.readline()  # title
        return int(f.readline().strip())

print(f"{'System':15s} {'gro atoms':>10s} {'status':>10s}")
all_ok = True
for sys in SYSTEMS:
    gro = f"{sys}/md_protein_ref.gro"
    xtc = f"{sys}/md_noWAT.xtc"
    try:
        n_gro = gro_atom_count(gro)
    except FileNotFoundError:
        print(f"{sys:15s} {'MISSING':>10s} {'FAIL':>10s}")
        all_ok = False
        continue

    try:
        u = mda.Universe(gro, xtc)
        status = "OK"
    except ValueError as e:
        status = "MISMATCH"
        all_ok = False
    print(f"{sys:15s} {n_gro:10d} {status:>10s}")

print("\nALL OK" if all_ok else "\nSOME SYSTEMS FAILED -- fix before running H-bond analysis")
