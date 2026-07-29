# CYP2B6 Allele Molecular Dynamics Analysis — Methods, Results, and Discussion

**Author:** Smyan Reddy
**Supervision:** Shaylyn Govender, Prof. Özlem Tastan Bishop (RUBi, Rhodes University)
**Date:** July 2026

## Methods

Eleven CYP2B6 alleles (G99E, K139E, M46V, I328T, I391N, K262R, R140Q, R487C,
S259R, T306S-R378K, and wild-type) were each simulated for 300 ns in two
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

## Results

**Global stability (RMSD).** Of the ten mutant alleles, I328T shows the
clearest and most consistent destabilization: both replicates maintain an
elevated RMSD plateau relative to both WT replicates for the entire 300 ns,
and its KDE density is visibly broader and right-shifted relative to WT with
minimal overlap — the strongest evidence in the panel for a mutation that
shifts the protein toward a distinct, less stable conformational ensemble.
By contrast, M46V, K139E, I391N, R140Q, and S259R each show a KDE density
shifted toward *lower* RMSD than WT (tighter, more rigid sampling), though in
each case this pattern held for only one of the two mutant replicates against
one of the two WT replicates, not a clean four-way separation. G99E, K262R,
R487C, and T306S-R378K show substantial density overlap with WT, i.e. no
strong evidence of an overall stability effect from RMSD alone.

**Local flexibility (RMSF).** A recurring flexibility peak at true residues
~136-140 (a loop connecting two alpha-helices, per DSSP secondary structure
assignment) appeared in the WT reference itself, but with inconsistent
replicate agreement: across the ten mutants analyzed, a different WT
replicate showed the higher peak at this loop in nine distinct combinations,
indicating this is largely inherent simulation-to-simulation variability
rather than an allele-specific effect. Two alleles nonetheless implicate this
same loop directly: K139E's mutation site sits inside it, and G99E's mutation
site (~65 residues away, at a beta-strand/loop junction) correlates with
flexibility there allosterically.

The strongest local-effect candidate is S259R: both replicates show RMSF
amplified above WT specifically at S259's own mutation site (a loop that is
already among the more flexible regions of the protein at baseline), making
it the only allele where a local flexibility increase is at least directionally
consistent across both replicates and localizes to the mutation's own
position rather than a distal site. The panel-wide RMSF heatmap additionally
highlights G99E as showing the single largest flexibility increase of any
allele (concentrated around residues 230-260, GROMACS numbering) and K262R
as showing a distinct, isolated flexibility increase near residue ~163 not
seen in any other allele.

All other alleles (M46V, I391N, R140Q, R487C, T306S-R378K) show flat RMSF at
their own mutation sites, with only single-replicate, unconfirmed local peaks
elsewhere in the sequence.

## Discussion

The RMSD/RMSF results provide two candidate mechanistic hypotheses worth
prioritizing for the next stage of analysis (hydrogen bonding and DRN):

1. **I328T** as a case of global destabilization — worth checking whether
   this correlates with loss of a key stabilizing hydrogen-bond network or a
   shift in dynamic residue network centrality near the active site.
2. **S259R** as a case of local, mutation-site flexibility — one of three
   alleles with currently uncertain functional classification, and the only
   one in this analysis with a plausible local structural mechanism
   (increased flexibility at an already-mobile loop containing the mutation
   itself).

A caution that should carry into every subsequent analysis: the recurring
~136-140 loop signal is present in WT itself and varies substantially between
WT replicates, meaning any claim that a specific allele perturbs this loop
needs to be tested against the reference-replicate spread (e.g., a 3-standard-
deviation threshold over WT replicates, as used in the CYP3A4 paper) rather
than a pairwise visual comparison. This same caution likely applies to any
region where WT replicate agreement is not first confirmed.

No conclusions are drawn here about the third uncertain allele,
T306S-R378K, beyond noting that neither of its two point mutations shows a
local RMSF effect, and its RMSD shows two non-overlapping single-replicate
excursions rather than one reproducible signal — it remains the least
resolved of the three uncertain alleles pending hydrogen-bond and DRN data.

## Next steps

1. Hydrogen bond frequency analysis, all 11 alleles vs. WT.
2. DRN analysis via MDM-TASK-web (betweenness, closeness, eigenvector
   centrality), prioritizing the three uncertain-function alleles (S259R,
   P428T, T306S-R378K).
3. Re-evaluate the ~136-140 loop and the I328T/S259R hypotheses above once
   H-bond and DRN data are available, applying the reference-triplicate
   significance threshold rather than pairwise comparison.
