#!/usr/bin/env python
"""
Part 8 binding-pocket analysis -- summary aggregation for all systems.

For each system this computes, from files already on disk:
  - Active-site pocket volume statistics from MDpocket characterization
    (mdpock_<SYS>_descriptors.txt, written by run_all_mdpocket.sh):
      * pock_volume: per-snapshot active-site volume (A^3); 0.00 means the
        selected active-site pocket is closed/absent in that snapshot --
        so the "open fraction" (fraction of snapshots with volume > 0) is a
        direct measure of how often the active site is accessible.
      * mean volume over all snapshots, mean volume over open snapshots
        only, and median open volume.
  - Heme H-bond contact count from the existing Part 2 output
    (hbond_pairs_<SYS>.csv): number of trajectory frames (of 30001) in
    which any protein residue forms an H-bond with a heme atom
    (resnames CM1/HM1/FE1, GROMACS resid 408). Summed contact count and
    number of unique heme-partner residues. This is the analog of the
    reference paper's heme-contact heatmap (Rehema et al. Fig 4E/5).
  - Heme COM <-> protein COM distance over the 601 subsampled frames
    (mean +/- std), quantifying heme positional stability within the
    protein (the analog of the paper's heme-position tracking).
  - Active-site RMSF: mean RMSF over the residues lining the WT
    active-site pocket (from selected_pocket_WT.pdb), taken from the
    existing Part 1 per-residue RMSF output. Uses the WT pocket lining as
    the fixed residue set so all systems are compared on identical atoms.

Writes:
  - pocket_<SYS>.csv            per-snapshot volume time series
  - pocket_summary_all.csv      one row per system, all metrics

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.
"""
import glob
import os

import numpy as np

SYSTEMS = ["WT", "WT_2", "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
           "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2",
           "R140Q", "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2",
           "S259R", "S259R_2", "T306S-R378K", "T306S-R378K_2"]

HEME_RESNAMES = {"CM1", "HM1", "FE1"}


def parse_descriptors(path):
    """snapshot, volume, asa, nb_AS arrays from mdpock descriptors.
    Header: snapshot pock_volume pock_asa pock_pol_asa pock_apol_asa
    pock_asa22 pock_pol_asa22 pock_apol_asa22 nb_AS ... -> volume=1,
    asa=2, nb_AS=8."""
    snap, vol, asa, nb = [], [], [], []
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            cols = line.split()
            snap.append(int(cols[0]))
            vol.append(float(cols[1]))
            asa.append(float(cols[2]))
            nb.append(int(cols[8]))
    return np.array(snap), np.array(vol), np.array(asa), np.array(nb)


def heme_hbonds(csv_path):
    """Sum of H-bond counts involving heme, and unique partner residues."""
    total, partners = 0, set()
    with open(csv_path) as fh:
        fh.readline()  # header
        for line in fh:
            c = line.split(",")
            if len(c) < 8:
                continue
            d_res, a_res, count = c[0].strip(), c[3].strip(), int(c[6])
            if d_res in HEME_RESNAMES:
                total += count
                partners.add(c[1].strip())
            if a_res in HEME_RESNAMES:
                total += count
                partners.add(c[4].strip())
    return total, len(partners)


def active_site_rmsf(sysdir, residues):
    """Mean per-residue RMSF over the WT active-site lining residues."""
    rmsf_files = glob.glob(os.path.join(sysdir, "rmsf_*.xvg"))
    if not rmsf_files:
        return np.nan
    arr = np.loadtxt(rmsf_files[0], comments=["#", "@"])
    mask = np.isin(arr[:, 0].astype(int), residues)
    if not mask.any():
        return np.nan
    return arr[mask, 1].mean()


def heme_drift(sysdir):
    """Mean (+/- std) heme COM to protein COM distance over subsampled frames."""
    import MDAnalysis as mda

    top = os.path.join(sysdir, "md_protein_ref.pdb")
    traj = os.path.join(sysdir, "md_pocket.xtc")
    if not (os.path.exists(top) and os.path.exists(traj)):
        return np.nan, np.nan
    u = mda.Universe(top, traj)
    heme = u.select_atoms("resname CM1 HM1 FE1")
    protein = u.select_atoms("protein")
    dists = np.empty(u.trajectory.n_frames)
    for i, ts in enumerate(u.trajectory):
        dists[i] = np.linalg.norm(heme.center_of_mass() - protein.center_of_mass())
    return dists.mean(), dists.std()


def main():
    # WT active-site lining residues: fixed reference set for all systems.
    wt_residues = sorted({int(l[22:26]) for l in open("WT/selected_pocket_WT.pdb")
                          if l.startswith(("ATOM", "HETATM"))})

    rows = []
    for sys in SYSTEMS:
        desc = os.path.join(sys, f"mdpock_{sys}_descriptors.txt")
        if not os.path.exists(desc):
            print(f"WARN {sys}: missing descriptors -- skipping")
            continue
        snap, vol, asa, nb = parse_descriptors(desc)
        n = len(vol)
        open_mask = vol > 0.0
        n_open = int(open_mask.sum())

        # Heme contacts (Part 2 output, full 30001-frame counts)
        hb_path = os.path.join(sys, f"hbond_pairs_{sys}.csv")
        if os.path.exists(hb_path):
            hb_sum, hb_partners = heme_hbonds(hb_path)
        else:
            hb_sum, hb_partners = np.nan, np.nan

        drift_mean, drift_std = heme_drift(sys)
        asrmsf = active_site_rmsf(sys, wt_residues)

        row = {
            "system": sys,
            "n_frames": n,
            "n_open": n_open,
            "open_frac": n_open / n if n else np.nan,
            "mean_vol_all": vol.mean() if n else np.nan,
            "mean_vol_open": vol[open_mask].mean() if n_open else np.nan,
            "median_vol_open": np.median(vol[open_mask]) if n_open else np.nan,
            "mean_asa_open": asa[open_mask].mean() if n_open else np.nan,
            "heme_hbond_sum": hb_sum,
            "heme_hbond_partners": hb_partners,
            "heme_drift_mean": drift_mean,
            "heme_drift_std": drift_std,
            "active_site_rmsf": asrmsf,
        }
        rows.append(row)

        # per-snapshot time series for plotting
        with open(os.path.join(sys, f"pocket_{sys}.csv"), "w") as fh:
            fh.write("snapshot,volume_A3\n")
            for s, v in zip(snap, vol):
                fh.write(f"{s},{v:.3f}\n")

        print(f"{sys:12s} n_open={n_open:4d}/{n:3d}  "
              f"mean_vol_all={row['mean_vol_all']:6.2f}  "
              f"mean_vol_open={row['mean_vol_open']:6.2f}  "
              f"hb_sum={hb_sum:6.0f}  drift={drift_mean:5.3f}+/-{drift_std:.3f}")

    cols = ["system", "n_frames", "n_open", "open_frac", "mean_vol_all",
            "mean_vol_open", "median_vol_open", "mean_asa_open",
            "heme_hbond_sum", "heme_hbond_partners", "heme_drift_mean",
            "heme_drift_std", "active_site_rmsf"]
    with open("pocket_summary_all.csv", "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")
    print(f"\nWrote pocket_summary_all.csv ({len(rows)} systems)")


if __name__ == "__main__":
    main()
