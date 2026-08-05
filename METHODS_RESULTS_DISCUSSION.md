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

**Dynamic residue network (DRN) analysis.** Betweenness (BC), closeness (CC),
and eigenvector (EC) centrality were computed per residue with `calc_network.py`
from MD-TASK (RUBi-ZA/MD-TASK, `mdm-task-web` branch — the RUBi group's own
tool, installed locally rather than via its web-upload interface since the
trajectories here are multi-gigabyte). The tool requires a PDB topology
rather than `.gro`/`.tpr`; each system's `md_protein_ref.gro` was converted
via `gmx editconf` (`convert_to_pdb.sh`). Network nodes are one representative
atom per residue (CB, or CA for glycine), connected by an edge when the two
representative atoms fall within the tool's default 6.7 Å cutoff. Metrics
were computed every 100 frames (~1 ns intervals across the 300 ns trajectory,
300 frames total) rather than every frame, based on a timing test showing the
full-resolution run would take prohibitively long across all 22 systems;
`run_all_drn.py` batched this across all systems in priority order.

A methodological issue was found and must be accounted for in any residue-
level DRN comparison: the heme cofactor is not fully excluded by the CB/CA-
based network reduction. Cross-referencing each system's output CSV against
its paired `.cif` structure file (which carries the true residue number as
`auth_seq_id`) confirmed that CSV row N corresponds exactly to GROMACS resid
N — but row/resid 408 is not a real amino acid; it is the heme cofactor's
CM1 component, which happens to have an atom named "CB" in this force field
and therefore survives the reduction step as a spurious 463rd network node
(462 real residues + 1 heme node). Heme's other components (HM1 hydrogens,
FE1 iron) do not have a CA/CB-named atom and are correctly excluded. This
mirrors, but is distinct from, the heme-adjacency issue found during the
H-bond analysis (there, heme resids sat *near* certain mutation-site windows;
here, heme is an actual row inside the per-residue centrality dataset).
`drn_significance_check.py` explicitly excludes resid 408 from both the
global (whole-protein mean) and local (mutation-site window) metrics, rather
than relying on no mutation site's window happening to reach it (none
currently do — the closest is P428T at 397-403 — but this is guarded
defensively regardless, following the same precaution taken for the heme-
adjacency issue in `hbond_significance_check.py`). The same
replicate-noise-floor robustness framework used for RMSD/RMSF/H-bonds was
applied to each of BC, CC, and EC: a global metric (mean centrality across
all real residues) and a local metric (window-max value in a +/-3 residue
window at each allele's own mutation site).

**Radius of gyration (Rg) and solvent-accessible surface area (SASA).**
Following the handover document's recommended analyses (RMSD/RMSF, Rg, SASA,
H-bonds, and several analyses not yet performed -- PCA, DCCM, conformational
clustering, binding-pocket/channel analysis -- listed in Next steps), Rg
(`gmx gyrate`) and SASA (`gmx sasa`) were computed for all 22 systems.
Unlike RMSD/RMSF, which used the Backbone index group for consistency, Rg
and SASA were computed on the Protein group (whole protein, all atoms):
Backbone-only would omit side chains, which matters for both overall
compactness and solvent-exposed surface, so the full protein is the
standard choice for these two metrics specifically -- a deliberate,
documented deviation from the Backbone convention used elsewhere, not an
inconsistency. Rg is a whole-molecule scalar with no per-residue breakdown,
so only a global metric applies (mean Rg over the trajectory). SASA has
both a global metric (mean total SASA) and a local, per-residue metric
(`gmx sasa -or`, window-max at +/-3 residues around each allele's own
mutation site, the SASA analogue of the windowed RMSF/H-bond/DRN checks);
the per-residue output has exactly 462 rows (one per real amino acid, no
heme contamination of the kind found in the DRN dataset, since the
"Protein" selection group excludes the heme cofactor entirely). The same
replicate-noise-floor robustness framework was applied to both metrics via
`rg_sasa_significance_check.py`.

**Conformational clustering.** Following Shaylyn Govender's predecessor MSc
thesis on this same system (Chapter 4), conformational clustering was
performed with AmberTools `cpptraj`: hierarchical agglomerative clustering
(average-linkage), cut to 3 clusters, using best-fit RMSD on backbone atoms
(`@C,CA,N,O`) as the distance metric. `cpptraj` cannot parse GROMACS `.tpr`
topologies directly (confirmed by a real failed run: "Could not determine
format of topology"), so `md_protein_ref.pdb` (the protein+heme-only PDB
already generated for the DRN analysis) was used as the topology instead,
matching `md_noWAT.xtc`'s atom composition exactly. Two deviations from her
exact protocol are documented rather than silently assumed: (1) her thesis
also used the heme Fe and "L-helix" as an additional distance metric, but
the L-helix residue range isn't available in a form that could be reliably
extracted from the thesis text, so this run uses Backbone RMSD alone,
matching this project's established Backbone convention for RMSD/RMSF; (2)
hierarchical clustering requires a full pairwise distance matrix (O(n^2)),
which is intractable at the full 30001 frames/system, so every 10th frame
was used (~3000 frames/system, the same subsampling reasoning as the DRN
analysis's `--step 100`). `cluster_summary.py` parses each system's
`cluster_<SYSTEM>_summary.dat` and applies the same replicate-noise-floor
robustness framework used throughout, using each system's dominant
(largest) cluster's frame fraction as the comparison metric -- a high value
means the trajectory is dominated by one conformational state; a lower
value split across multiple substantial clusters means more than one
distinct conformation was sampled, the same quantity reported in Shaylyn's
Table 4.1.

**Secondary structure (DSSP).** Following the second half of Shaylyn
Govender's predecessor thesis Chapter 4 ("Analysis through Clustering and
DSSP"), per-residue secondary structure was assigned via the Kabsch-Sander
algorithm. An initial attempt used AmberTools `cpptraj`'s built-in
`secstruct` command (already set up for the clustering step above), but this
produced unusable results: every one of the 463 real residues across all 22
systems showed near-zero Alpha/3-10/Extended/Bridge/Pi content, with almost
everything assigned to the generic "Bend" category — implausible for a
heavily alpha-helical P450 like CYP2B6. The root cause was traced to
`md_protein_ref.pdb`'s bonding, which `cpptraj` infers from atom distances
rather than reading from a real topology (this system has no AMBER prmtop —
it was simulated in GROMACS, `md.tpr`/`cyp2b6_GMX.top`); this likely broke or
missed backbone peptide-bond connectivity between some residues, preventing
`cpptraj`'s H-bond-ladder-based helix/sheet detection from working. This was
caught by a sanity check (a real protein's DSSP output should show
substantial helix/sheet content, not none), not assumed to be a valid null
result. The fix was to switch to GROMACS's own `gmx dssp` (`run_all_dssp_gmx.sh`,
group 1/"Protein"), which both resolves the bonding issue (reads the correct
GROMACS topology directly) and exactly matches Shaylyn's original method,
making these results directly comparable to her thesis. `gmx dssp -o`
produces one line per frame, one character per residue (standard DSSP
one-letter codes: H/G/I = alpha/3-10/pi helix, E/B = strand/bridge, T/S/P =
turn/bend/polyproline, ~ = coil, = = unassignable gap — heme's CM1 at resid
408 falls in this last category and is excluded). `dssp_summary.py` computes,
per residue, the fraction of the 30001-frame trajectory spent in any defined
helix or strand/bridge state ("ordered secondary structure" = H+G+I+E+B),
writing a small per-system per-residue CSV (`dssp_<SYSTEM>_orderedss.csv`)
rather than tracking the full ~14MB/system per-frame matrix in git. The same
replicate-noise-floor robustness framework used throughout was applied: a
global metric (mean ordered-SS fraction across all 463 residues) and a local
metric (window-max |delta| in a +/-3 residue window at each allele's own
mutation site).

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

**DRN centrality (BC, CC, EC).** At the global (whole-protein mean) level,
essentially no allele shows a robust effect in any of the three centrality
metrics — expected, since a single point mutation rarely shifts a
whole-protein average computed over 462 residues. The few that nominally
pass (S259R for BC; G99E, I328T, K262R, R487C, P428T for EC) sit close to
their respective noise floors and should be treated cautiously rather than
as strong global findings.

At the local (mutation-site window) level, robust effects are: **BC** —
G99E, K139E, M46V, R140Q, S259R; **CC** — M46V, I391N; **EC** — R487C (and
the R378 site of T306S-R378K). Notably, **P428T and I328T — the two alleles
prioritized for DRN analysis on the strength of their convergent
RMSD/RMSF/H-bond evidence — show no robust local DRN signal in any of the
three metrics.** P428T's local deltas are the largest in the panel for all
three metrics (BC -0.012, CC -0.007, EC -0.015), but its own two replicates
disagree in direction in every case (e.g. BC: rep1 -0.0267 vs rep2 +0.0020),
so despite the large magnitude this does not clear the robustness bar. This
is reported as a genuine negative result rather than downplayed: at the
network-centrality layer specifically, the two highest-priority alleles from
Parts 1-2 do not show a reproducible signal at their own mutation site.

**Radius of gyration and SASA.** At the global level, robust Rg *compaction*
(more compact than WT) is seen in K139E, M46V, R140Q, P428T, and S259R.
Robust global SASA *decrease* is seen in M46V, I391N, R140Q, S259R, and
T306S-R378K (whole allele, not resolved to either individual site at the
global level). At the local (per-residue, own-site window) level, robust
SASA effects are: M46V (decrease), I391N (decrease), S259R (increase, in the
opposite direction from its own global decrease), and T306S-R378K's R378
site (increase).

Two patterns stand out. First, **M46V** now has robust findings across four
independent measures -- global RMSD compaction, local DRN BC and CC
increases, Rg compaction, and both global and local SASA decrease -- the
most convergent evidence of distributed, non-local stabilization of any
allele in the panel, despite showing no RMSF or H-bond signal at its own
site. Second, **S259R**, the allele with zero robust findings through Parts
1-2, now has three robust findings across Rg/SASA (compaction, global SASA
decrease, local SASA increase) in addition to its earlier DRN-only local BC
finding -- and its own-site SASA increase moving opposite to its whole-
protein SASA decrease is a real, specific pattern (a locally more exposed
mutation site inside an overall more compact, less solvent-exposed protein)
worth a closer look rather than a data artifact to wave away.

P428T's robust Rg compaction is worth flagging alongside its already-
established elevated global RMSD: these are not contradictory (RMSD
measures deviation from a reference conformation; Rg measures compactness
of whatever conformation is currently sampled), but it does mean P428T's
picture is more complex than "more flexible than WT" alone -- it deviates
more from the reference structure while adopting an overall more compact
shape.

**Conformational clustering.** The two WT replicates disagree sharply on
dominant-cluster-fraction (0.454 vs. 0.808 -- WT_2 is dominated by a single
conformational state nearly twice as strongly as WT), giving a noise floor
(0.354) larger than every single allele's replicate-averaged delta from WT.
**No allele passes the robustness check for this metric.** M46V (0.315),
R140Q (0.314), and I328T (0.244) have the largest raw deltas, and all three
happen to overlap with alleles Shaylyn's thesis flagged as showing the most
cluster/structural deviation (I328T, K262R, P428T, R140Q) -- but none of
these deltas clears this project's own noise floor, and per-allele
replicate agreement is inconsistent (e.g. P428T and S259R have *opposite-
signed* deltas between their own two replicates). This is reported as a
genuine null result for this specific metric rather than cherry-picked into
a positive finding: with only two WT replicates, a single hierarchical
clustering run's dominant-cluster fraction is evidently too sensitive to
trajectory-specific sampling noise to support any robust per-allele
conclusion here. This does not mean the underlying MD data lacks signal --
RMSD/RMSF/H-bonds/DRN/Rg/SASA all found robust effects in several of these
same alleles -- it means this particular clustering metric, run this way,
isn't discriminating enough given the available replicate count.

**Secondary structure (DSSP).** At the global level (mean ordered-SS fraction
across all 463 residues), the WT rep1-vs-rep2 noise floor is 0.0208 — larger
than every single allele's replicate-averaged delta (all under 0.013).
**No allele passes the global robustness check.** This is consistent with
expectation: a single point mutation is unlikely to shift the whole-protein
fold-element balance measurably, and no earlier analysis (RMSD, RMSF,
H-bonds, DRN, Rg/SASA, clustering) found a robust *global* structural
collapse or expansion either.

At the local (own mutation-site window) level, three alleles pass the
robustness check with large, clearly-signed effects: **M46V** (delta =
-0.180, both replicates decreased — a robust local loss of ordered
secondary structure at its own site), **I391N** (delta = +0.464, both
replicates increased — a robust local *gain* of ordered structure), and
**K262R** (delta = +0.016, both replicates increased, smaller magnitude but
still clears its own tight noise floor of 0.004). T306S-R378K's T306 site
also passes (delta = -0.046, both replicates decreased), while its R378 site
does not. No other allele's own site passes; several show large raw local
deltas (G99E +0.240, R140Q +0.234, P428T +0.349) but fail because their two
replicates move in opposite directions at that window (e.g. P428T: rep1
+0.746, rep2 -0.047) — the same failure pattern already seen for P428T at
the DRN layer.

M46V's local ordered-SS *decrease* is a new, independent finding at its own
mutation site — notable because M46V previously showed no robust RMSF or
H-bond effect at that same site, only robust *global* RMSD/Rg/SASA
compaction and robust *local* DRN (BC, CC) increases elsewhere in the
protein. This adds a fifth independent measure to M46V's case, and is the
first of the five to show a robust effect *at the mutation site itself*
rather than only in a distributed/global sense — a modest local unfolding
event co-occurring with (not necessarily causing) the protein-wide
compaction. I391N's local ordered-SS *increase* is likewise a new finding at
its own site, adding to its existing robust global/local SASA decreases —
three independent methods now agree I391N has some real, if subtle, local
structural effect. K262R's local ordered-SS increase adds a fourth angle on
a site that already showed the panel's largest robust local RMSF
rigidification, though no robust H-bond or DRN effect there; a residue
locally stiffening (RMSF) while gaining defined secondary structure (DSSP)
is at least directionally consistent, even though the two metrics were not
designed to move together.

## Discussion

Applying a real robustness check rather than visual comparison substantially
changed which findings can be trusted, and in one case (S259R) reversed a
conclusion reported earlier. The findings that survive replicate-consistency
and noise-floor testing are:

1. **P428T** was the strongest case across RMSD/RMSF/H-bonds (elevated global
   RMSD, elevated local RMSF and H-bond count at its own site), but this
   convergence does **not** extend to DRN centrality: none of BC, CC, or EC
   pass the robustness check at its own site, and its two replicates
   disagree in direction on all three despite showing the largest-magnitude
   deltas in the panel. Read together, this now looks like a site with real,
   reproducible local backbone flexibility and transient H-bonding, but
   without a corresponding reproducible shift in network centrality/topology
   at that residue — the mutation may be affecting local dynamics without
   rewiring the broader residue-interaction network in a consistent way.
2. **I328T** similarly showed convergent RMSD/RMSF/H-bond evidence but no
   robust DRN finding at its own site in any of the three centrality
   metrics. As with P428T, the RMSD/RMSF/H-bond case stands on its own but
   is not corroborated by network centrality.
3. **G99E** is a new finding that first emerged from the H-bond analysis (a
   robust local H-bond increase at its own site, without a robust RMSD/RMSF
   effect), and now also shows a robust local BC increase at its own site —
   the first allele with a robust finding across two independent methods
   without also having a robust RMSD/RMSF signal. This strengthens G99E as
   worth continued attention, on different grounds than P428T/I328T.
4. **K262R** shows a robust local rigidification at its own mutation site by
   RMSF (the largest-magnitude local RMSF effect in the panel) without a
   corresponding global RMSD effect or a robust local H-bond effect, and
   without a robust DRN finding at its own site either — a mutation whose
   local stiffening is not explained by, or reflected in, either the
   H-bond or centrality data collected so far. DSSP now adds a robust local
   ordered-secondary-structure *increase* at the same site (smaller
   magnitude, but real), directionally consistent with the RMSF
   rigidification even though the two methods aren't formally linked.
5. **M46V** shows a robust global RMSD compaction/rigidification without a
   corresponding local RMSF or H-bond effect at its own mutation site, and
   robust local BC *and* CC increases at its own site (the only allele with
   robust local hits in two of the three DRN metrics). Rg/SASA now add
   robust global Rg compaction and robust global *and* local SASA decrease —
   four independent measures now converge on the same picture: a
   distributed, non-local stabilizing/compacting effect, with no
   corresponding RMSF or H-bond signal at the mutation site itself. DSSP now
   adds a fifth measure and, notably, the first one that IS local to the
   mutation site itself: a robust local ordered-secondary-structure
   *decrease*, meaning M46V shows some real local unfolding right where the
   protein-wide compaction is also happening — a more complete picture than
   "purely allosteric/distributed" alone.
6. **T306S-R378K**: its R378 site (not T306) continues to be the only part
   of this allele with robust findings: a local H-bond decrease (delta =
   -0.571), a robust global and local EC increase at the R378 site
   specifically, and now a robust local SASA increase at R378 too. T306
   shows no robust finding in any metric across any of the four analyses.
   The allele's global SASA (whole double mutant, not resolved by site) also
   shows a robust decrease. This continues to narrow the case for this
   allele's functional effect specifically to the R378 substitution.
7. **K139E, R140Q, R487C** — no robust RMSD/RMSF or H-bond finding; each
   shows exactly one robust DRN result at its own site (K139E and R140Q:
   local BC; R487C: local EC). R140Q additionally now shows a robust global
   Rg compaction. These remain single-network-metric findings (plus, for
   R140Q, one global Rg result) and should be treated as provisional leads.
8. **I391N** — no robust RMSD/RMSF or H-bond finding, one robust DRN result
   (local CC), and robust global and local SASA decreases — and now a robust
   local ordered-secondary-structure *increase* at its own site by DSSP,
   the largest-magnitude local DSSP effect in the panel. Three independent
   methods now agree I391N has a real, subtle local structural/compacting
   effect, a step up from a single-metric provisional finding.
9. **S259R** — the allele with zero robust findings through Parts 1-2, then
   one DRN-only finding (local BC), now has three more: robust Rg
   compaction, robust global SASA decrease, and a robust *local* SASA
   *increase* at its own site — moving in the opposite direction from the
   whole-protein SASA decrease. That specific pattern (a more solvent-
   exposed mutation site sitting inside an overall more compact, less
   exposed protein) is a real, reproducible signature worth investigating
   further (e.g. whether the mutation site sits at a loop or surface patch
   that locally unfolds/opens while the rest of the protein compacts),
   rather than being dismissed as noise.

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
I328T's convergent multi-metric evidence. S259R was the allele with no
robust finding by any metric through Parts 1-2 (RMSD, RMSF, or H-bond); DRN
now gives it its first robust finding (local BC increase at its own site),
though this is a single-method result and should be treated with the same
caution as the other single-metric DRN-only findings above.

## Next steps

1. ~~DRN analysis via MDM-TASK-web~~ — done (Part 3, above). Headline result:
   P428T and I328T's strong RMSD/RMSF/H-bond convergence does not extend to
   DRN centrality (no robust finding at their own sites in BC/CC/EC), while
   several alleles with no prior robust finding (K139E, R140Q, S259R, I391N,
   R487C) each pick up exactly one robust DRN-only result, and M46V and G99E
   strengthen with a second independent robust method.
2. Request or approximate a proper reference-triplicate 3-SD threshold (the
   CYP3A4 paper's framework) rather than continuing to rely on a 2-replicate
   noise floor, which remains a workable but weaker substitute — this now
   applies across all three analyses (RMSD/RMSF, H-bonds, DRN).
3. Consider whether the single-method DRN-only findings (K139E, R140Q,
   R487C) merit a follow-up look at the actual network structure (e.g.
   which specific edges/neighbors drive the local centrality shift) rather
   than resting on the summary centrality value alone, before writing these
   up as findings in their own right.
4. ~~Rg and SASA~~ — done (Part 4, above). Headline results: M46V now has
   four independent convergent robust findings (RMSD, DRN-BC/CC, Rg, SASA),
   the strongest distributed/allosteric case in the panel after P428T/I328T;
   S259R goes from zero robust findings to four (DRN-BC plus Rg compaction,
   global SASA decrease, and an opposite-direction local SASA increase at
   its own site); I391N gains a second independent method (SASA) alongside
   its DRN-only finding; P428T's new Rg compaction adds nuance (not
   contradiction) to its established elevated-RMSD finding.
5. ~~Rg and SASA~~ / ~~conformational clustering~~ — done (Parts 4-5, above).
   Clustering headline: with only 2 WT replicates, dominant-cluster-fraction
   is too noisy a metric to call any allele robust (noise floor 0.354 larger
   than every allele's delta), despite the largest raw deltas (M46V, R140Q,
   I328T) overlapping with alleles Shaylyn's thesis flagged independently.
   Reported as a genuine null result for this metric, not evidence the
   underlying data lacks signal (it doesn't -- see the other five analyses).
6. ~~DSSP secondary structure analysis~~ — done (Part 6, above), after
   discovering and fixing a broken first attempt (cpptraj's `secstruct`
   produced near-zero helix/sheet content everywhere due to bad backbone
   connectivity inference from a distance-only PDB; switched to `gmx dssp`,
   which also matches Shaylyn's exact method). Headline: no allele shows a
   robust *global* effect (expected), but three show robust *local* effects
   at their own mutation site: M46V (ordered-SS decrease, its first-ever
   effect localized to the mutation site itself, adding to its four existing
   distributed/global findings), I391N (ordered-SS increase, a third
   independent method now agreeing on a real local effect), and K262R
   (ordered-SS increase, directionally consistent with its existing robust
   RMSF rigidification at the same site). Remaining analyses from the
   handover document's "Recommended Analyses" list, still not started: PCA,
   DCCM (dynamic cross-correlation matrices), and binding pocket/substrate
   access channel analysis.
7. If more WT replicates become available at any point, re-run
   `cluster_summary.py`'s robustness check -- a 2-replicate noise floor is
   an especially weak substitute specifically for the clustering metric
   (see above), more so than for the scalar/per-residue metrics elsewhere
   in this project.
8. S259R's local-vs-global SASA divergence (site more exposed while the
   whole protein is more compact) is worth a closer, targeted look --
   possibly via DSSP or a residue-level RMSF/SASA overlay at that specific
   site -- rather than folding it into a generic future-analysis queue item.
9. Request or approximate a proper reference-triplicate 3-SD threshold (the
   CYP3A4 paper's framework) rather than continuing to rely on a 2-replicate
   noise floor, which remains a workable but weaker substitute — this now
   applies across all five analyses completed so far.
10. Decide on final reporting format/figures for all five parts (RMSD/RMSF,
    H-bonds, DRN, Rg/SASA, clustering) for the supervision meeting write-up.
