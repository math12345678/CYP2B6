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
- `significance_check.py` — quantitative robustness check for both RMSD and
  RMSF mutation-site deltas, using WT rep1-vs-rep2 disagreement as an
  empirical noise floor and requiring both mutant replicates to agree in
  direction. Added after a visual-only read produced an incorrect conclusion
  (see summary below) — run this before trusting any directional claim in
  this analysis.

## Panel-wide RMSD/RMSF summary (all 11 alleles, statistically checked)

A quantitative robustness check (`significance_check.py`) was applied after
an initial visual-comparison pass produced at least one wrong conclusion
(see below): a mutant effect only counts as "robust" if both of its own
replicates agree in direction and the replicate-averaged delta exceeds the
disagreement between the two WT replicates themselves (the best available
noise floor with only 2 WT replicates). Everything below reflects that
checked version, not the original visual read.

- **P428T** (uncertain-function allele): the strongest signal in the panel,
  and the only allele with *both* a robust global RMSD elevation (+0.048 nm,
  both replicates) *and* a robust local RMSF elevation at its own mutation
  site (+0.058 nm, both replicates). Resolves one of the three
  uncertain-function alleles.
- **I328T**: robust global RMSD elevation (+0.035 nm) plus a smaller but still
  robust local RMSF increase at its own site (+0.026 nm) — the local signal
  was missed in the first pass (called "flat") until the windowed check.
- **K262R**: robust local RMSF *rigidification* at its own mutation site
  (-0.057 nm, the largest-magnitude local effect in the panel), with no
  corresponding global RMSD effect.
- **M46V**: robust global RMSD compaction (-0.025 nm, both replicates) with no
  corresponding local site effect — easy to miss visually since it had been
  informally lumped with several other alleles that only *looked* similar;
  it's the only one of that group that's actually statistically real.
- **Retracted:** S259R's "amplified local RMSF at its own site" and R140Q's
  "local rigidification at its own site" — both were reported in earlier
  passes of this analysis based on visual/single-point comparison, and both
  do not survive the robustness check (S259R's replicates disagree in
  direction once peak-shift is accounted for; R140Q's own site sits inside a
  region with WT-replicate noise far larger than any mutation could plausibly
  contribute). Both alleles are inconclusive for RMSF, not directional.
- **G99E** and **K139E** both implicate the same loop (true residues ~136-140):
  G99E allosterically (~65 residues away), K139E directly (mutation site sits
  inside the loop) — but this loop's own WT-replicate noise floor is large
  enough that no allele-specific claim about it is made without a proper
  reference-triplicate framework.
- Remaining alleles (I391N, R487C, T306S-R378K) show no robust local or
  global effect by this check.

## KDE and heatmap readout (all 11 alleles, complete)

- **KDE**: P428T and I328T show the clearest broadened/right-shifted densities
  relative to WT — the two strongest "distinct conformational ensemble" cases
  in the panel, with P428T's separation the more complete of the two. M46V,
  K139E, I391N, R140Q, and S259R all *visually* show a mutant peak shifted
  toward lower RMSD than WT rep1, but the quantitative check
  (`significance_check.py`) shows only M46V's shift is actually
  replicate-consistent and beyond the WT-replicate noise floor — the other
  four in that group do not pass the robustness check. G99E, K262R, R487C,
  and T306S-R378K show substantial overlap with WT — no clean density
  separation.
- **Heatmap**: G99E shows the single strongest deviation block (dark red,
  ~residues 230-260 GROMACS numbering — more flexible than WT on average),
  closely followed by P428T's own strongest deviation near residues ~400-415,
  close to its true mutation site. K262R has its own distinct hotspot near
  residue ~163 (red) plus a rigidifying band near ~230 (blue); R140Q shows a
  similar blue band at ~230. The 108-112 hotspot loop column is muted/mixed
  across alleles rather than one consistent color, supporting the read that
  this loop's signal is WT replicate noise rather than an allele-specific
  effect.

## Hydrogen bond analysis

GROMACS's own H-bond tools do not work on this dataset: `gmx hbond-legacy`
segfaults during its grid search, and the newer `gmx hbond` reports 0
donors/acceptors regardless of selection syntax (these 2018.6-generated
topologies appear to lack element metadata the newer tool needs). Switched
to MDAnalysis's `HydrogenBondAnalysis` instead, which identifies donors,
hydrogens, and acceptors by atom name/element rather than partial charge.

Two topology fixes were required first (see `fix_topology.sh`,
`check_topology.py`):
- The raw `.tpr` describes the full system (protein+water+ions), but the
  pre-stripped `md_noWAT.xtc` is protein+heme only. A matching protein+heme
  subset topology (`md_protein.tpr` / `md_protein_ref.gro`) is generated per
  system via `gmx convert-tpr` (group `Protein_CM1_HM1_FE1`) + `gmx trjconv`.
- MDAnalysis's bond-guessing has no default van der Waals radius for the
  heme iron (Fe); supplied manually in `hbond_analysis.py`.

Run order: `fix_topology.sh` → `check_topology.py` → `run_all_hbonds.py`
(batch driver over all 22 systems, calls `hbond_analysis.py` per system) →
`hbond_significance_check.py` (same WT-replicate-noise-floor robustness
framework as `significance_check.py`, applied to a global metric — mean
H-bonds/frame — and a local metric — summed pair frequency in a +/-3 residue
window around each allele's mutation site).

**Results:** no allele shows a robust *global* H-bond count change (WT
noise floor = 11.14 bonds/frame, larger than nearly every allele's effect).
At the *local* (mutation-site) level, four results are robust: **P428T**
(+0.498), **I328T** (+0.542), **G99E** (+0.361, a new finding not seen in
RMSD/RMSF), and **T306S-R378K's R378 site** (-0.571, its first robust
finding by any metric). P428T and I328T's H-bond results corroborate their
existing robust RMSD/RMSF findings; K262R's RMSF rigidification does not
have a matching robust H-bond effect. Full writeup in
`METHODS_RESULTS_DISCUSSION.md`.

## Next steps (per July 29 meeting with Prof. Bishop and Shaylyn)

1. ~~Finish RMSD (KDE) and RMSF (heatmap) figures, improve labeling~~ — done.
2. ~~Hydrogen bond analysis~~ — done, see above.
3. DRN analysis via MDM-TASK-web, prioritizing P428T and I328T (convergent
   multi-metric evidence), then T306S-R378K (one narrow finding at R378),
   then S259R (no robust finding yet by any metric).
