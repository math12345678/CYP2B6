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

## Dynamic residue network (DRN) analysis

Per-residue betweenness (BC), closeness (CC), and eigenvector (EC) centrality
computed with MD-TASK (`RUBi-ZA/MD-TASK`, `mdm-task-web` branch — RUBi's own
tool), installed locally via `git clone -b mdm-task-web ...` +
`conda env create -f environment.yml` (env `mdmtaskweb`, separate from the
`cyp2b6` env). `calc_network.py` requires a PDB topology (`convert_to_pdb.sh`
converts each system's `md_protein_ref.gro` via `gmx editconf`), and writes
outputs to whatever the current working directory is regardless of any path
prefix passed to it — `run_all_drn.py` runs each system with `cwd` set to
that system's own folder to avoid output collisions. Centrality computed
every 100 frames (~1 ns steps, 300 frames/system) across all 22 systems.

**Heme node caveat:** the CB/CA-based network-node reduction does not fully
exclude the heme cofactor. Confirmed via cross-reference against each
system's `.cif` output (which carries the true residue number as
`auth_seq_id`): CSV row N = GROMACS resid N, but resid 408 is heme's CM1
component (it has an atom named "CB" in this force field), not a real amino
acid — 462 real residues + 1 heme node = 463 total rows. `drn_significance_check.py`
explicitly excludes resid 408 from both the global and local metrics.

Run order: `convert_to_pdb.sh` → `run_all_drn.py` (batch driver, all 22
systems) → `drn_significance_check.py` (same WT-replicate-noise-floor
framework, applied to a global metric — mean centrality across all real
residues — and a local metric — window-max value at each allele's own
mutation site — for each of BC/CC/EC).

**Results:** essentially no allele shows a robust *global* effect in any of
the three metrics. At the *local* level: **P428T and I328T — the two
alleles with the strongest prior convergent RMSD/RMSF/H-bond evidence — show
no robust DRN finding at their own site in any of BC/CC/EC** (P428T has the
largest-magnitude local deltas in the panel, but its two replicates disagree
in direction on all three). Several alleles with no prior robust finding
pick up exactly one robust DRN result at their own site: K139E, R140Q, S259R
(local BC), I391N (local CC), R487C (local EC). M46V picks up two (local BC
and CC), and G99E strengthens with a robust local BC increase alongside its
existing H-bond finding. Only `*_mean.csv`/`*_mean_*.cif` summary files are
committed (per-frame `.dat` output and the cloned `MD-TASK/` tool are
gitignored). Full writeup in `METHODS_RESULTS_DISCUSSION.md`.

## Radius of gyration (Rg) and SASA analysis

Per the handover document's "Recommended Analyses" list (RMSD/RMSF, Rg,
SASA, and H-bonds under "Global stability"/"local structural effects" --
several of which, PCA/DCCM/clustering/binding-pocket analysis, are still
outstanding, see Next steps). Computed with `gmx gyrate` and `gmx sasa`
(`run_all_rg_sasa.sh`, run from the terminal directly since GROMACS isn't
available in the sandbox used for the Python-based analyses). Unlike
RMSD/RMSF (Backbone group), Rg/SASA use the Protein group (whole protein,
all atoms) -- side chains matter for both compactness and solvent-exposed
surface, so Backbone-only would be the wrong choice here specifically.
`rg_sasa_significance_check.py` applies the same WT-replicate-noise-floor
framework: Rg has only a global metric (whole-molecule scalar); SASA has
both global (mean total area) and local (per-residue, window-max at each
allele's own site) metrics.

**Results:** robust Rg compaction in K139E, M46V, R140Q, P428T, S259R.
Robust global SASA decrease in M46V, I391N, R140Q, S259R, T306S-R378K.
Robust local SASA effects: M46V and I391N (decrease), S259R and
T306S-R378K's R378 site (increase). Two findings stand out: **M46V** now has
four independent convergent robust results (RMSD compaction, DRN local
BC/CC, Rg compaction, SASA decrease) -- the cleanest distributed/allosteric
case in the panel besides P428T/I328T. **S259R**, which had zero robust
findings through Parts 1-2, now has four (DRN local BC, Rg compaction,
global SASA decrease, and a local SASA *increase* at its own site -- the
opposite direction from the whole-protein decrease, a real and specific
pattern worth a closer look, not noise). Full writeup in
`METHODS_RESULTS_DISCUSSION.md`.

## Conformational clustering

Following Shaylyn Govender's predecessor MSc thesis on this same system
(Chapter 4): AmberTools `cpptraj` hierarchical agglomerative clustering
(average-linkage, 3 clusters), best-fit RMSD on backbone atoms
(`@C,CA,N,O`) as the distance metric, every 10th frame (~3000 frames/system,
same subsampling reasoning as DRN's `--step 100`). `cpptraj` cannot parse
GROMACS `.tpr` directly (confirmed by a real failed run), so uses
`md_protein_ref.pdb` (already generated for DRN) as topology instead.
Deviation from Shaylyn's exact protocol, documented: her thesis also used
heme Fe + "L-helix" as an additional distance metric, but the L-helix
residue range wasn't recoverable from the available thesis text, so this
run uses Backbone RMSD alone. `cluster_summary.py` applies the same
WT-replicate-noise-floor framework to each system's dominant (largest)
cluster's frame fraction, the same quantity in Shaylyn's Table 4.1.

**Results:** the two WT replicates disagree sharply (dominant fraction 0.454
vs. 0.808), giving a noise floor (0.354) larger than every allele's
replicate-averaged delta. **No allele passes the robustness check for this
metric** -- reported as a genuine null result, not evidence the underlying
MD data lacks signal (RMSD/RMSF/H-bonds/DRN/Rg/SASA all found real effects
in several of these alleles). The largest raw deltas (M46V 0.315, R140Q
0.314, I328T 0.244) do overlap with alleles Shaylyn's thesis flagged
independently (I328T, K262R, P428T, R140Q), which is a suggestive
coincidence worth keeping in mind, but none of them clear this project's own
noise floor and per-allele replicate agreement is inconsistent (P428T and
S259R even have opposite-signed deltas between their own two replicates).
Full writeup in `METHODS_RESULTS_DISCUSSION.md`.

## Secondary structure (DSSP)

Following Shaylyn Govender's predecessor thesis (Chapter 4, second half).
First attempt used AmberTools `cpptraj`'s `secstruct` command, but this
produced unusable output: near-zero helix/sheet content across every residue
in every system, implausible for a heavily alpha-helical P450. Traced to
broken backbone bonding inference in `md_protein_ref.pdb` (this is a GROMACS
simulation with no AMBER prmtop, so `cpptraj` had to guess bonds from atom
distances, and got peptide-bond connectivity wrong somewhere). Switched to
GROMACS's own `gmx dssp` (`run_all_dssp_gmx.sh`, group 1/"Protein"), which
fixes the bonding issue and also exactly matches Shaylyn's original method.
Output is one line/frame, one DSSP letter code/residue (H/G/I helix, E/B
strand, T/S/P turn-ish, ~ coil, = heme gap). `dssp_summary.py` reduces this
to a per-residue "ordered secondary structure" fraction (H+G+I+E+B) and
writes a small per-system CSV (`dssp_<SYSTEM>_orderedss.csv`) rather than
tracking the full ~14MB/system per-frame matrix; applies the same
WT-replicate-noise-floor framework (global: mean across all 463 residues;
local: window-max at each allele's own mutation site).

**Results:** no allele passes the *global* robustness check (WT noise floor
0.0208, larger than every allele's delta) — expected, a single point
mutation shouldn't move the whole-protein fold balance. At the *local*
(own-site) level, three alleles pass: **M46V** (decrease, -0.180 — its
first-ever effect localized to the mutation site itself, on top of four
existing distributed/global findings), **I391N** (increase, +0.464, the
largest local DSSP effect in the panel — a third independent method now
agreeing I391N has a real local effect), and **K262R** (increase, +0.016,
smaller but clears its own tight noise floor — directionally consistent with
its existing robust RMSF rigidification at the same site). T306S-R378K's
T306 site also passes (decrease, -0.046); R378 does not. Several alleles show
large raw local deltas (G99E, R140Q, P428T) but fail because their two
replicates move in opposite directions at that window — the same failure
pattern already seen for P428T at the DRN layer. Full writeup in
`METHODS_RESULTS_DISCUSSION.md`.

## Next steps (per July 29 meeting with Prof. Bishop and Shaylyn, and the
## handover document's Recommended Analyses)

1. ~~Finish RMSD (KDE) and RMSF (heatmap) figures, improve labeling~~ — done.
2. ~~Hydrogen bond analysis~~ — done, see above.
3. ~~DRN analysis via MDM-TASK-web~~ — done, see above.
4. ~~Rg and SASA~~ — done, see above.
5. ~~Conformational clustering~~ — done, see above (genuine null result for
   this metric with only 2 WT replicates).
6. ~~DSSP secondary structure analysis~~ — done, see above (after fixing a
   broken cpptraj-based first attempt). Still outstanding from the handover
   document's Recommended Analyses: PCA, DCCM, binding pocket/substrate
   access channel analysis.
7. Consider whether the single-metric DRN-only findings (K139E, R140Q,
   R487C) warrant a deeper look at network structure before being written up
   as standalone findings.
8. S259R's local-vs-global SASA divergence (own site more exposed while the
   whole protein is more compact) deserves a targeted follow-up look.
9. Decide on final reporting format/figures for all six completed parts for
   the supervision meeting write-up.
