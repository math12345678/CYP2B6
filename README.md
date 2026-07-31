# CYP2B6 Allele MD Analysis

Reproducibility log for RMSD/RMSF analysis of 11 CYP2B6 alleles vs. wild-type,
following the CYP3A4 paper's (Mwaniki et al.) pipeline. All commands below were
run exactly as shown, from each system's extracted trajectory folder.

## Setup

- GROMACS 2026.1 (trajectories originally produced with GROMACS 2018.6)
- Structure: PDB 5UFG-based homology model, MODELLER 10.4, first modeled residue
  ARG at true/UniProt position 29 (residues 1-28 not resolved/modeled)
- Force field: AmberFF14SB, TIP3P water
- Production runs: 300 ns each, 2 independent replicates per system (WT + 10 mutant alleles)
- Analysis environment: conda (mdanalysis, mdtraj, ambertools/cpptraj), Python 3
  (numpy, matplotlib, scipy)

## Important: residue numbering offset

GROMACS renumbers residues sequentially starting at 1 in the simulation topology
(`md.gro`/`.tpr`), rather than preserving the original PDB/UniProt numbering,
which starts at 29 (since residues 1-28 are unresolved in the structure).

**Confirmed offset: GROMACS residue = true/UniProt residue - 28** (verified as a
uniform shift by comparing residue identities at multiple positions between the
raw PDB and the `.gro` topology, including confirming the G99E mutation site
literally becomes `GLH` (protonated Glu) at GROMACS residue 71, i.e. 99 - 28).

**All mutation-site coordinates used in figure legends throughout this repo are
GROMACS numbering unless explicitly labeled "true"/"UniProt".** Always convert
true allele nomenclature (e.g. G99E) by subtracting 28 before marking it on a
GROMACS-derived RMSF plot.

## Commands run (identical pattern for every system)

For each system (WT, WT_2, G99E, G99E_2, K139E, K139E_2, M46V, M46V_2, I328T,
I328T_2, I391N, I391N_2, K262R, K262R_2, R140Q, R140Q_2, R487C, R487C_2, S259R,
S259R_2, T306S-R378K, T306S-R378K_2), run from inside the extracted system
folder:

```bash
gmx rms -s md.tpr -n index.ndx -f md_noWAT.xtc -o rmsd_<SYSTEM>.xvg -tu ns
gmx rmsf -s md.tpr -n index.ndx -f md_noWAT.xtc -o rmsf_<SYSTEM>.xvg -res
```

Both commands prompt for a group selection twice (least-squares fit group, then
calculation group). **Group 4 (Backbone) was selected in every case**, for
consistency across all systems.

## Known data issues encountered and resolved

- WT.tar/WT_2.tar were initially missing from the shared Drive folder; confirmed
  via independent `find`/`tar -tf` verification before flagging to Shaylyn, who
  re-uploaded them.
- P428T_2.tar was corrupted (`tar -tf` "Truncated input file"); flagged and
  re-uploaded by Shaylyn.
- Several large trajectory extractions timed out when run directly against the
  Google Drive–synced folder mid-download; fixed by copying/extracting `.tar`
  archives to a local (non-Drive-synced) folder instead.
- Disk space was exhausted partway through the panel (~100+ GB of raw
  `md.xtc`/`md_noPBC.xtc` trajectories accumulating); resolved by deleting the
  large raw trajectory files (keeping only `md_noWAT.xtc` and the `.xvg`
  results) for each system once its analysis was complete.

## Analysis outputs

- `plot_<ALLELE>.py` — per-allele WT-vs-mutant line plots (RMSD over time,
  RMSF per residue), generated for all 11 alleles.
- `plot_kde_rmsd_all.py` — KDE density plots for RMSD across all 11 alleles vs.
  WT, per Prof. Bishop's request (July 29 meeting) — line plots show
  equilibration, density plots better reveal whether a mutant is sampling a
  distinct/broader conformational ensemble than WT.
- `plot_rmsf_heatmap.py` — single heatmap of RMSF delta (mutant avg - WT avg)
  across all 11 alleles and all residues, per Prof. Bishop's request.

## Panel-wide RMSD/RMSF summary (all 11 alleles, complete)

- **P428T** (uncertain-function allele): the strongest overall signal in the
  panel. Both replicates show sustained elevated RMSD vs. both WT replicates
  for essentially the entire 300 ns (cleaner, more complete separation than
  I328T, which partially converged back to WT range in one replicate). Its KDE
  density is the most clearly shifted/broadened of any allele. The RMSF
  heatmap additionally shows P428T's single largest per-allele deviation
  (residues ~400-415, GROMACS numbering) sitting close to its own true
  mutation site (residue 428 / 400 GROMACS) — the combination of a clean
  global destabilization plus a local flexibility signal near the mutation
  site makes this the most complete "mutation has a measurable effect" case
  found so far, and it resolves one of the three previously uncertain-function
  alleles.
- **I328T**: also a clear global stability effect — sustained elevated RMSD vs.
  both WT replicates for most of the 300 ns.
- **S259R** (uncertain-function allele): both replicates show amplified local
  RMSF right at their own (shared, inherently flexible) mutation site.
- **G99E** and **K139E** both implicate the same loop (true residues ~136-140):
  G99E allosterically (~65 residues away), K139E directly (mutation site sits
  inside the loop).
- The 108-112 (GROMACS)/136-140 (true) loop shows high WT-replicate variability
  across nearly every mutant tested (a different WT replicate is "the high one"
  in most comparisons) — this baseline noise means any allele-specific claim
  about this loop needs the reference-triplicate + 3 SD significance threshold
  (per the CYP3A4 paper's framework), not pairwise-curve reading.
- **K262R and R140Q** show a real, quantified local *rigidification* at their
  own mutation sites (delta vs. WT: K262R -0.067 nm — the largest mutation-site
  RMSF change in the entire panel; R140Q -0.051 nm) — opposite direction from
  S259R's flexibility increase. This was caught during a numeric audit of
  mutation-site RMSF deltas after an earlier informal visual read had
  mischaracterized both as "flat, no effect."
- Remaining alleles (M46V, I391N, R487C, T306S-R378K) show no clear local or
  global effect, or single-replicate-only signals that need more data before
  treating as reproducible.

## KDE and heatmap readout (all 11 alleles, complete)

- **KDE**: P428T and I328T show the clearest broadened/right-shifted densities
  relative to WT — the two strongest "distinct conformational ensemble" cases
  in the panel, with P428T's separation the more complete of the two. M46V,
  K139E, I391N, R140Q, and S259R all show a mutant peak shifted toward lower
  RMSD than WT rep1, often tighter/narrower. G99E, K262R, R487C, and
  T306S-R378K show substantial overlap with WT — no clean density separation.
- **Heatmap**: G99E shows the single strongest deviation block (dark red,
  ~residues 230-260 GROMACS numbering — more flexible than WT on average),
  closely followed by P428T's own strongest deviation near residues ~400-415,
  close to its true mutation site. K262R has its own distinct hotspot near
  residue ~163 (red) plus a rigidifying band near ~230 (blue); R140Q shows a
  similar blue band at ~230. The 108-112 hotspot loop column is muted/mixed
  across alleles rather than one consistent color, supporting the read that
  this loop's signal is WT replicate noise rather than an allele-specific
  effect.

## Next steps (per July 29 meeting with Prof. Bishop and Shaylyn)

1. Finish RMSD (KDE) and RMSF (heatmap) figures, improve labeling (bigger axis
   labels/titles/legends, one figure legend per comparison stating which
   mutant vs. WT).
2. Hydrogen bond analysis (next method in the pipeline).
3. DRN analysis via MDM-TASK-web.
