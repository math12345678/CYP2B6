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

## Discussion

Applying a real robustness check rather than visual comparison substantially
changed which findings can be trusted, and in one case (S259R) reversed a
conclusion reported earlier. The findings that survive replicate-consistency
and noise-floor testing are:

1. **P428T** remains the strongest case in the panel, and is now on firmer
   footing than before: it is the only allele with *both* a robust global
   RMSD elevation *and* a robust local RMSF elevation at its own mutation
   site, with both replicates agreeing in direction on both metrics. This is
   one of the three uncertain-function alleles and should be the top
   priority for H-bond and DRN follow-up.
2. **I328T** is the second-strongest case — robust global RMSD elevation
   plus a smaller but still robust local RMSF increase at its own site
   (previously reported as "flat" before the windowed check). Worth checking
   whether this correlates with loss of a stabilizing hydrogen-bond network
   or a DRN centrality shift near the active site.
3. **K262R** shows a robust local rigidification at its own mutation site
   (the largest-magnitude local effect in the panel) without a corresponding
   global RMSD effect — a mutation that stiffens an otherwise mobile site
   rather than destabilizing the protein overall. Worth checking via
   hydrogen bonding whether this corresponds to a newly formed or
   strengthened stabilizing contact.
4. **M46V** shows a robust global RMSD compaction/rigidification without a
   corresponding local signal at its own mutation site — a candidate for an
   allosteric or distributed stabilizing effect rather than one localized to
   the mutation position, worth checking against DRN centrality changes.

**Retracted from earlier versions of this analysis:** S259R's local RMSF
"amplification" claim (the two replicates actually disagree in direction once
peak-shift is properly accounted for) and R140Q's local RMSF "rigidification"
claim (based on an unreliable single-point comparison inside a very noisy
region). Both alleles' own mutation sites should currently be treated as
inconclusive for RMSF, not directional.

No conclusions are drawn about the third uncertain allele, T306S-R378K:
neither of its two point mutations shows a robust local RMSF effect, and its
global RMSD shows replicate-disagreeing excursions rather than one
reproducible signal. It remains the least resolved of the three uncertain
alleles pending hydrogen-bond and DRN data.

## Next steps

1. Hydrogen bond frequency analysis, all 11 alleles vs. WT.
2. DRN analysis via MDM-TASK-web (betweenness, closeness, eigenvector
   centrality), prioritizing the three uncertain-function alleles in this
   order given the RMSD/RMSF evidence: P428T first (strongest, multi-metric
   signal), then S259R and T306S-R378K (both currently inconclusive and in
   most need of an independent method to resolve).
3. Once H-bond and DRN data exist, request or approximate a proper
   reference-triplicate 3-SD threshold (the CYP3A4 paper's framework) rather
   than continuing to rely on a 2-replicate noise floor, which is a workable
   but weaker substitute.
4. Re-run this same robustness check (`significance_check.py`) on H-bond and
   DRN metrics once available, rather than reverting to visual comparison.
