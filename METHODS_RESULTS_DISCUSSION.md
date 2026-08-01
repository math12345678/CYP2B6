# CYP2B6 Allele Molecular Dynamics Analysis — Methods, Results, and Discussion

**Author:** Smyan Reddy
**Supervision:** Shaylyn Govender, Prof. Özlem Tastan Bishop (RUBi, Rhodes University)
**Date:** July 2026

## Methods

Eleven CYP2B6 mutant alleles (G99E, K139E, M46V, I328T, I391N, K262R, R140Q,
R487C, P428T, S259R, T306S-R378K), plus wild-type, were each simulated for 300 ns in two
independent replicates, using a homology model built from PDB 5UFG (MODELLER
10.4), the AmberFF14SB force field, and TIP3P water, in GROMACS. Trajectories
were provided pre-run and post-processed (PBC removal, water stripping); this
analysis covers everything from that point forward.

For every system, backbone RMSD and per-residue backbone RMSF were computed
against each trajectory's own reference structure using `gmx rms` and
`gmx rmsf -res`, selecting the Backbone atom group in both cases. RMSD was
computed over the full 300 ns without excluding an equilibration window,
following a conservative default pending confirmation from Shaylyn on
convention. Per Prof. Bishop's guidance, RMSD was additionally visualized as
kernel density estimates (KDE) rather than only time-series line plots, since
density plots more directly reveal whether a mutant is sampling a distinct
conformational ensemble relative to WT (a sharp, well-overlapping WT density
across both replicates versus a broader or shifted mutant density). RMSF
results across the full panel were consolidated into a single heatmap,
plotting the delta between each mutant's replicate-averaged RMSF and the
WT replicate average, following the same delta-based significance framework
used in the CYP3A4 reference paper (Mwaniki et al.).

One methodological issue was identified and corrected during this analysis:
GROMACS renumbers protein residues sequentially starting at 1 in the
simulation topology, while the source structure (and true UniProt/PharmVar
allele nomenclature) begins at residue 29, since residues 1-28 are not
resolved in the model. This was confirmed as a uniform +28 offset (true
residue = GROMACS residue + 28) by cross-referencing residue identities
between the source PDB and the `.gro` topology at multiple positions,
including direct confirmation that the G99E mutation appears as GLH
(protonated glutamate) at GROMACS residue 71 (99 - 28). All mutation-site
markers in this analysis use this corrected mapping.

**Statistical robustness check.** An initial pass through this analysis
characterized mutants as "elevated," "flat," or "rigidified" relative to WT
based on visual comparison of overlaid line plots. This was found to be
unreliable on its own (see Results) and was replaced with a quantitative
check (`significance_check.py`): since only two WT replicates are available
here (rather than the reference triplicate used for the CYP3A4 paper's 3-SD
threshold), the disagreement between the two WT replicates themselves was
used as an empirical noise floor. A mutant effect is only reported as
"robust" if (a) both of its own replicates agree in the direction of the
change relative to WT, and (b) the replicate-averaged delta exceeds this
noise floor. For RMSF mutation-site comparisons specifically, a +/-3 residue
window maximum was used instead of the single-point value at the exact
mutation-site residue, since RMSF peaks can shift by a residue or two between
systems and a single-point comparison proved sensitive to that shift (this
was discovered when it produced a result contradicted by the windowed
version for the same allele). All findings below use this robustness check;
anything that does not pass it is explicitly labeled inconclusive rather than
described directionally.

**Hydrogen bond analysis.** GROMACS's own hydrogen-bond tools could not be
used on this dataset: the deprecated `gmx hbond-legacy` segfaulted during its
grid-search step, and the newer `gmx hbond` (GROMACS 2024+) reported zero
donors and zero acceptors regardless of selection syntax, apparently because
these 2018.6-generated topologies lack the per-atom element metadata the
newer tool's donor/acceptor detection depends on. H-bond analysis was
instead performed with MDAnalysis's `HydrogenBondAnalysis`, which identifies
donors/acceptors/hydrogens by atom name and element rather than by partial
charge or embedded element records (both of which this topology format also
lacks in a form MDAnalysis could read directly). Two further topology issues
were resolved before this would run: (1) the raw `.tpr` files describe the
full system (protein + water + ions), while the pre-stripped
`md_noWAT.xtc` trajectories contain only protein + the heme cofactor
(CM1/HM1/FE1 residues) — a protein+heme-only subset topology was generated
per system via `gmx convert-tpr` (selecting group `Protein_CM1_HM1_FE1`)
and `gmx trjconv`, confirmed to match the trajectory's atom count before
use; (2) MDAnalysis's bond-guessing (needed since the `.gro` topology
carries no bond records) has no default van der Waals radius for the heme
iron (Fe), which was supplied manually. Donor/hydrogen/acceptor selections
were specified explicitly as `element H` (hydrogens) and `element O or
element N` (donors and acceptors), since MDAnalysis's automatic guessing
defaults to a partial-charge heuristic this topology cannot supply. H-bonds
were computed over the full 300 ns trajectory (donor-acceptor distance
cutoff 3.5 Å, donor-H-acceptor angle cutoff 150°) for all 22 systems (11
alleles + WT, 2 replicates each), producing a per-frame H-bond count and a
per-pair occurrence frequency for each system. The same replicate-noise-floor
robustness framework used for RMSD/RMSF (`significance_check.py`) was applied
to these results via `hbond_significance_check.py`: a global metric (mean
H-bonds/frame, compared the same way as RMSD) and a local metric (the sum of
donor-acceptor pair frequencies where either residue falls within a +/-3
residue window of the allele's own mutation site, the H-bond analogue of the
windowed RMSF check).

## Results

**Global stability (RMSD).** Using the robustness check described above
(WT rep1-vs-rep2 noise floor = 0.0189 nm), only three of the eleven mutant
alleles show a replicate-consistent RMSD deviation from WT that exceeds this
noise floor: **P428T** (delta = +0.048 nm, both replicates elevated),
**I328T** (delta = +0.035 nm, both replicates elevated), and **M46V**
(delta = -0.025 nm, both replicates decreased/more rigid). P428T's separation
is the most complete of the two elevated cases — both replicates stay
elevated for essentially the full 300 ns, whereas I328T's rep1 partially
converges back toward WT range around 100-170 ns. M46V's robust global
compaction was not obvious from visual inspection alone (it had been grouped
informally with several other alleles that also appeared to trend toward
lower RMSD) but is the only one of that group whose replicates actually agree
in direction and clear the noise floor. All other alleles (G99E, K139E,
I391N, K262R, R140Q, R487C, S259R, T306S-R378K) do not show a
replicate-consistent RMSD effect that exceeds the noise floor — several
looked directionally suggestive on the KDE grid, but the two mutant
replicates either disagreed in direction (R140Q, T306S-R378K) or the
averaged effect size did not clear what two WT replicates alone would differ
by.

**Local flexibility (RMSF) at each allele's own mutation site.** Applying the
same robustness check (+/-3 residue window, both replicates must agree in
direction and exceed the WT rep1-vs-rep2 noise floor at that window) to every
allele's own mutation site: **P428T** shows a robust local flexibility
increase (delta = +0.058 nm), **K262R** shows a robust local flexibility
*decrease*/rigidification (delta = -0.057 nm, the largest-magnitude local
effect in the panel), and **I328T** shows a smaller but still robust local
increase (delta = +0.026 nm). No other allele's own mutation site passes both
criteria. Notably, **S259R** does not: an earlier informal read of the line
plots described "both replicates show amplified RMSF at S259's own site," but
the windowed, replicate-by-replicate check shows the two S259R replicates
actually disagree in direction (rep1 delta -0.029, rep2 delta +0.026) once
peak-shift is accounted for. That claim is retracted here; S259R's own
mutation site is inconclusive, not amplified. **R140Q** is similarly
inconclusive at its own site: it sits adjacent to a separately-recurring
flexibility peak (true residues ~136-140, discussed below) with very high
WT-replicate disagreement there (noise floor of 0.135-0.220 nm depending on
window), which swamps any signal a single point mutation could contribute.
An earlier version of this document reported R140Q as showing "a real
decrease" at its own site based on a single-point (non-windowed) comparison;
that comparison did not account for peak-shift or the loop's inherent noise
and is also retracted.

A recurring flexibility peak at true residues ~136-140 (a loop connecting two
alpha-helices, per DSSP secondary structure assignment) appears in the WT
reference itself, with high replicate-to-replicate disagreement across nearly
every mutant tested — consistent with this being inherent simulation
variability rather than an allele-specific effect. Two alleles' own mutation
sites sit at or near this loop (K139E directly inside it, G99E ~65 residues
away allosterically), but given the loop's own noise floor is far larger than
typical mutation-site deltas seen elsewhere in this panel, no allele-specific
claim about this loop is made here without the reference-triplicate framework
the CYP3A4 paper uses.

**Hydrogen bonding.** At the global level (mean H-bonds/frame across the
whole protein), the WT rep1-vs-rep2 noise floor is 11.14 bonds/frame —
large relative to most alleles' apparent effects (typical deltas of
1-7 bonds/frame), and no allele's global H-bond count passes the robustness
check: every allele's two replicates either disagree in direction outright,
or the averaged effect does not clear the noise floor. This mirrors the
RMSD result in one sense (most alleles show no robust global effect) but is
a stricter finding: even M46V and S259R, whose replicate-averaged deltas are
among the largest (+6.84 and +4.59 bonds/frame respectively), fail because
their own two replicates move in the same direction only loosely and do not
clear 11.14.

At the local (mutation-site window) level, four results are robust: **P428T**
(delta = +0.498, both replicates agree), **I328T** (delta = +0.542, both
replicates agree), **G99E** (delta = +0.361, both replicates agree), and
**T306S-R378K's R378 site specifically** (delta = -0.571, both replicates
agree; its paired T306 site does not pass, and the allele's global H-bond
count also does not pass). No other allele's own mutation site — including
K262R, despite its robust RMSF rigidification — shows a robust local H-bond
effect.

Two of these four (P428T, I328T) corroborate their existing robust RMSD/RMSF
findings with an independent method: both show elevated local flexibility at
their own site (RMSF) and an increased local H-bond count there, consistent
with a picture of local backbone mobility accompanied by additional
transient hydrogen bonding (rather than, say, loss of a stabilizing contact).
G99E is a new finding not previously flagged by RMSD/RMSF: it did not show a
robust effect by either of those metrics, but does show a robust local
H-bond increase at its own site (true residue 99, GROMACS residue 71).
K262R's RMSF rigidification does not have a corresponding robust local
H-bond signal, so that finding remains supported by RMSF alone rather than
by convergent evidence across methods.

## Discussion

Applying a real robustness check rather than visual comparison substantially
changed which findings can be trusted, and in one case (S259R) reversed a
conclusion reported earlier. The findings that survive replicate-consistency
and noise-floor testing are:

1. **P428T** remains the strongest case in the panel, now supported by three
   independent, mutually-corroborating robust findings: elevated global RMSD,
   elevated local RMSF at its own site, and elevated local H-bond count at
   its own site. This is one of the three uncertain-function alleles and
   should be the top priority for DRN follow-up.
2. **I328T** is the second-strongest case, also now with convergent evidence
   across three metrics: robust global RMSD elevation, a smaller but robust
   local RMSF increase (previously reported as "flat" before the windowed
   check), and a robust local H-bond increase at its own site. Worth
   checking whether this correlates with a DRN centrality shift near the
   active site.
3. **G99E** is a new finding that only emerges from the H-bond analysis: it
   showed no robust global or local effect by RMSD/RMSF, but does show a
   robust local H-bond increase at its own site. Since this is a single-
   method finding so far, it should be treated as provisional pending DRN
   data rather than placed on the same footing as P428T/I328T.
4. **K262R** shows a robust local rigidification at its own mutation site by
   RMSF (the largest-magnitude local RMSF effect in the panel) without a
   corresponding global RMSD effect or a robust local H-bond effect — a
   mutation that stiffens an otherwise mobile site by some mechanism other
   than a straightforward gain of hydrogen bonding at that site. Worth
   checking against DRN centrality changes, since the rigidification is not
   explained by the H-bond data collected so far.
5. **M46V** shows a robust global RMSD compaction/rigidification without a
   corresponding local RMSF signal or a robust local H-bond effect at its own
   mutation site — a candidate for an allosteric or distributed stabilizing
   effect rather than one localized to the mutation position, worth checking
   against DRN centrality changes.
6. **T306S-R378K**: while its global RMSD and global H-bond count both fail
   the robustness check, its R378 site specifically shows a robust local
   H-bond *decrease* (delta = -0.571) not seen at its paired T306 site. This
   is the first robust finding of any kind for this previously fully
   inconclusive allele, and narrows the open question to the R378 mutation
   specifically rather than the double mutant as a whole.

**Retracted from earlier versions of this analysis:** S259R's local RMSF
"amplification" claim (the two replicates actually disagree in direction once
peak-shift is properly accounted for) and R140Q's local RMSF "rigidification"
claim (based on an unreliable single-point comparison inside a very noisy
region). Both alleles' own mutation sites should currently be treated as
inconclusive for RMSF, not directional.

Neither of the two remaining point mutations in T306S-R378K shows a robust
local RMSF effect, and its global RMSD and global H-bond count both show
replicate-disagreeing or noise-floor-bound excursions rather than one
reproducible signal. However, its R378 site does show a robust local H-bond
decrease (see above) — the first robust signal for this allele by any
metric, though a narrower one (single site, single metric) than P428T or
I328T's convergent multi-metric evidence. S259R remains the allele with no
robust finding by any metric collected so far (RMSD, RMSF, or H-bond).

## Next steps

1. DRN analysis via MDM-TASK-web (betweenness, closeness, eigenvector
   centrality), prioritizing the three uncertain-function alleles in this
   order given the RMSD/RMSF/H-bond evidence so far: P428T and I328T first
   (both have convergent multi-metric support), then T306S-R378K (one
   narrow robust finding, at R378 only), then S259R (no robust finding yet
   by any metric).
2. G99E's new H-bond-only finding (robust local increase at its own site)
   should be checked against DRN data before being treated as more than
   provisional, since it is not yet corroborated by a second method.
3. Once DRN data exists, request or approximate a proper reference-triplicate
   3-SD threshold (the CYP3A4 paper's framework) rather than continuing to
   rely on a 2-replicate noise floor, which is a workable but weaker
   substitute.
4. Re-run this same robustness check (`significance_check.py` /
   `hbond_significance_check.py`) on DRN metrics once available, rather than
   reverting to visual comparison.
