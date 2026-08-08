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

**Principal component analysis (PCA) and dynamic cross-correlation matrix
(DCCM).** Following the handover document's remaining "Recommended
Analyses," both were computed from GROMACS's `gmx covar` covariance-matrix
machinery, but on different atom selections for a documented reason: PCA
used the Backbone group (matching this project's RMSD/RMSF convention) via
`gmx covar` + `gmx anaeig` (least-squares fit and covariance/eigenvector
calculation both on Backbone), giving the eigenvalue spectrum (variance
captured by each collective motion, descending) and a PC1/PC2 trajectory
projection. DCCM used the C-alpha group specifically (462 atoms = 462 real
residues; heme has no atom literally named "CA," so it is automatically and
completely excluded from this group without any manual exclusion step, a
cleaner situation than some earlier per-residue datasets in this project
that did include a heme placeholder row), fit on Backbone first to remove
rigid-body translation/rotation (`-ascii` full covariance matrix output).
`pca_dccm_summary.py` reduces the raw (3x462)^2 ascii covariance matrix to
a standard 462x462 residue-residue correlation matrix (`C_ij =
trace(Cov_ij) / sqrt(trace(Cov_ii)*trace(Cov_jj))`) entirely in memory --
the full matrix is not written to disk, only derived summary statistics,
consistent with this project's practice of tracking small summaries rather
than large raw outputs. Two PCA metrics (PC1 eigenvalue, PC1 variance
fraction) and two DCCM metrics (global mean |correlation| across all
residue pairs; local mean |correlation| between each allele's own mutation
site and every residue more than 3 positions away, excluding trivially-
correlated backbone neighbors) were each run through the same replicate-
noise-floor robustness framework used throughout.

A residue-numbering subtlety required care here: the C-alpha group (462
atoms) is numbered contiguously 1-462 with no gap at the heme's position,
unlike the 463/464-row datasets elsewhere in this project (DRN, DSSP) whose
numbering runs 1-463 through a heme placeholder at resid 408. Any mutation
site above resid 408 needs 1 subtracted before indexing into the C-alpha
array; of this project's sites, only R487C (resid 459) is affected
(T306S-R378K's sites, 278 and 350, are both below 408). This adjustment is
applied via a small helper function and documented in-script. Checking this
also surfaced a genuine bug in the earlier Rg/SASA local-metric script
(`rg_sasa_significance_check.py`), which was not making this adjustment;
fixed alongside this work (see Next steps) -- re-running confirmed no
conclusions actually changed as a result, but the underlying indexing logic
was wrong and has been corrected.

**Binding-pocket (active-site) analysis.** Following the reference CYP3A4
paper's MDpocket protocol (Rehema et al., JMB 2025) -- the final item on the
handover document's "Recommended Analyses" list. Two rounds per system, per
that paper: Round 1, whole-protein *pocket exploration*, where MDpocket over
the full trajectory produces pocket density/frequency grids
(`mdpout_dens_grid.dx`/`mdpout_freq_grid.dx`, viewable in VMD) showing where
transient pockets open and close; Round 2, *pocket characterization*, where
MDpocket is restricted to the active-site pocket to produce the per-frame
pocket volume time series (the analog of the paper's Fig 4D). The full
30001-frame water-stripped trajectory was subsampled to every 50th frame (601
frames at 0.5 ns cadence, matching the paper's 1001-snapshot-per-500-ns
density) with cpptraj, using `md_protein_ref.pdb` (protein+heme) as topology.
`fpocket` then detected pockets on the reference structure, and
`select_active_pocket.py` selected the pocket whose centroid is closest to the
heme Fe (validated against the canonical CYP2B6 active-site lining residues,
e.g. the heme-flanking residues previously identified). Round 2 ran MDpocket
with `--selected_pocket`, writing `mdpock_<SYS>_descriptors.txt` (per-frame
`pock_volume` in A^3; 0.00 means the active-site pocket is closed/absent in
that snapshot, so the open fraction directly measures active-site
accessibility). Round 1 ran MDpocket without the selected pocket to write the
whole-protein grids.

A significant tooling issue had to be solved before exploration would run at
all: the conda-forge `mdpocket` binary crashes at startup (Trace/BPT trap,
EXC_BREAKPOINT inside `libsystem_malloc`'s `mfm_alloc`) whenever the process
is launched with a random ASLR layout. The characterization invocation (with
`--selected_pocket`) survived ~60-70% of launches and a retry loop absorbed
the rest, but the exploration invocation (no selected pocket) failed ~100% of
direct launches. Running the identical command under `lldb` (which disables
ASLR by default) is much more reliable; the remaining crashes (an `rpdb_read`
null-FILE* on a failed fopen, or the same `mfm_alloc` trap) are absorbed by
the per-system retry loop in `run_exploration_lldb.sh`. Two other
practicalities: exploration writes fixed output names (`mdpout_*`), so two
exploration jobs must never run in the same system directory (the batch driver
runs systems in parallel only across different dirs), and the grids are only
written at the *end* of a successful run -- so a 0-byte grid file is a failed
run, and the scripts check for a non-empty file (`-s`), not mere existence
(`-f`). Under parallel load (8 systems at once) the exploration step was being
killed by a 1500-s timeout before it could finish; the fixup driver
(`run_all_exploration.sh`) reruns only the systems with empty/missing grids at
lower parallelism (P=3) and a 5400-s timeout.

`pocket_summary.py` aggregates, per system, from files already on disk:
active-site volume statistics (mean over all snapshots where closed snapshots
count as 0 -- combining "how big" with "how often open"; mean/median over open
snapshots only), the open fraction, protein-heme H-bond contact count (from
the Part 2 `hbond_pairs_<SYS>.csv`, summed over all 30001 frames and the
number of unique heme-partner residues), heme COM-to-protein COM distance over
the 601 subsampled frames ("heme drift", quantifying heme positional stability
-- the analog of the paper's heme-position tracking), and the mean RMSF over
the WT active-site lining residues (from `selected_pocket_WT.pdb`, used as a
fixed reference residue set so all systems are compared on identical atoms).
The same WT-replicate-noise-floor robustness framework used throughout was
applied via `pocket_significance_check.py` to six metrics -- mean_vol_all,
mean_vol_open, open_frac, heme_hbond_sum, heme_drift_mean, active_site_rmsf --
with metric-specific minimum effect sizes (2-3 A^3 for volumes, 0.03 for open
fraction, 3 frames for H-bond counts, 0.01 nm for drift/RMSF) to prevent the
tiny-noise-floor false-positive class documented earlier.
`pocket_plots.py` generates the figure set
(`pocket_volume_all_alleles.png`, `pocket_volume_timeseries_WT.png`,
`pocket_heme_contacts_all_alleles.png`, `active_site_rmsf_all_alleles.png`).

**Substrate-access-channel analysis (CAVER 3.0).** The reference CYP3A4
paper explicitly includes substrate channel dynamics among the analyses used
to characterize allele effects, so the substrate-access-channel geometry was
added as a second round of the binding-pocket item, using **CAVER 3.0**
(Chovancova et al., PLoS Comput Biol 2012, 8:e1002708) -- the established
tool for CYP substrate-channel analysis (CYP2B6: IterTunnel, J Cheminform
2014; CYP2D6: PLoS One 2014; CYP3A4: PLoS One 2024). All 24 systems
(WT/WT_2 plus 11 alleles x 2 replicates) were analyzed with no failures.
Snapshots: 121 structures per system prepared by `prep_caver_snapshots.py`
from the same 601-frame water-stripped trajectory used by MDpocket, taking
every 5th frame (2.5 ns cadence). The starting point was each system's
buried active-site cavity centroid (the WT centroid [48.27, 43.87, 42.05]
plus the per-system heme-Fe offset from WT; `starting_points.tsv`,
`config_caver_template.txt`). A starting point on the bare heme Fe was tried
first and yields zero tunnels (CAVER warns the point lies in a buried
cavity); the cavity-centroid start resolves this and matches published CAVER
CYP practice. CAVER 3.0 needs a working JVM; it runs under the `cyp2b6`
conda env with `openjdk=8` (Zulu macos-aarch64 build), driven by
`caver_analysis/run_all_caver.sh` (4 parallel systems). Per snapshot,
`analysis/tunnel_characteristics.csv` gives the tunnel count, widest
bottleneck radius and length; a channel is "open" when its widest bottleneck
radius >= 1.3 A (the van der Waals radius of water). `analyze_caver_results.py`
aggregates over the full expected 121-snapshot set (CAVER omits rows for
snapshots with zero tunnels, so the denominator is the fixed frame range
1..601 step 5, not the observed rows) and applies the same
WT-replicate-noise-floor robustness framework, writing `caver_summary_all.csv`;
`caver_plots.py` makes `caver_open_frac_all_alleles.png` (open fraction, mean
widest bottleneck over all snapshots, mean bottleneck over open snapshots).


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

**DRN network-structure deep-dive (K139E, R140Q, R487C).** The three
single-metric DRN-only findings above (README Next-steps item 9) were
re-examined at the per-residue level to distinguish a *coherent neighborhood
effect* from a *single-residue artifact* — the same question the S259R SASA
follow-up had to answer, and with the same failure mode available. Per-residue
centrality deltas were built from the per-timepoint `.dat` traces (301
timepoints/system, the same data the committed `*_mean.csv` files average),
using the WT-replicate disagreement as the per-residue noise floor (same
convention as `s259r_sasa_followup.py`). Each mutation-site ±3 window was
scored two ways: window-max (the original significance-check metric) and
window-mean (a spike-invariant test of whether the neighborhood moves
together rather than one residue moving alone). Per-residue movers were
ranked by |Δ| above the noise floor and annotated against literature-defined
CYP2B6 regions (`drn_network_deepdive.py`, figure
`drn_network_deepdive.png`).

**K139E (CYP2B6\*8, local BC) — a coherent neighborhood effect.** Both window
tests pass: window-max Δ = +0.00849 (noise floor 0.00151) and window-mean Δ
= +0.00118 (noise floor 0.00004), both replicates agreeing in sign. 141/462
residues (30.5%) move robustly, and the |Δ| mass is delocalized (±3 = 2%,
±10 = 5%, ±20 = 8%) — expected for a centrality metric, and not the signature
of a single-residue spike. The strongest movers are BC *decreases* in the
C-terminal L helix (true 442-448: 442, 445, 446 in the top movers) and in
SRS-5 (true 362, 365). True 443 is R443, the CYP2B6 ortholog of CYP2B4 R443,
one of the proximal-face CPR-binding residues mapped by alanine-scanning
(Bridges et al. 1998, JBC 273:17036-17049; R443 orthology verified here by
direct CYP2B4↔CYP2B6 global alignment — the full set R122/R126/R133/K139/K422/
K433/R443 maps 1:1 onto CYP2B6, and 3IBD HELIX 20 = Gly438-Asn456 places
R443 in the L helix). K139 itself is a Bridges et al. CPR contact located in
the C/D loop (3IBD C helix 117-135, D helix 141-160 → loop 136-140), and the
K139E charge reversal (CYP2B6\*8) is known to abolish electron transfer from
CPR while leaving the heme/active site catalytically competent (Zhang et al.
2011, JPET 338:803-809; stopped-flow rate ~30-fold slower than WT, 77% 7-EFC
O-deethylase activity retained with tBHP as oxidant). The DRN reading — a
BC re-wiring that connects the C/D-loop mutation site to the L-helix
CPR-binding surface — is therefore structurally coherent and
literature-consistent, not a spike artifact.

**R140Q (CYP2B6\*14, local BC) — a coherent neighborhood effect.** Both window
tests pass: window-max Δ = +0.00816 (noise 0.00151) and window-mean Δ =
+0.00249 (noise 0.00007). 143/462 (31.0%) robust movers, same delocalized
|Δ| mass profile. Uniquely among the three, the top movers are dominated by
active-site lining residues: V367 and L363 (SRS-5), F297 and T302 (SRS-4),
and I114 (SRS-1) are all members of the Angle & Cox (2023) active-site 5 Å
set, plus F369 (SRS-5). This is consistent with R140Q's independent robust
active-site widening in the pocket analysis (+98 Å³ mean volume, +0.77 open
fraction — the strongest pocket change in the panel): the BC re-wiring
reaching the active-site lining residues and the enlarged, nearly-always-open
active site are mutually corroborating. R140 sits immediately adjacent to
K139 in the same C/D loop, and R140Q (CYP2B6\*14) is a documented variant
with reduced expression/function (Lang et al. 2004, JPET 311:34-43).

**R487C (CYP2B6\*5, local EC) — a single-residue artifact, the same failure
mode as S259R's window-max SASA.** Window-max Δ = −0.00506 (noise 0.00193) is
robust, but window-mean Δ = −0.00205 (noise 0.00133) is NOT: the two R487C
replicates disagree in sign on the window mean. Only 71/462 residues (15.4%)
move robustly, and |Δ| mass is even more delocalized (±3 = 1%, ±10 = 3%,
±20 = 4%). The *local* EC finding at R487C's own site is therefore downgraded
from a neighborhood finding to a spike at the mutated residue itself. The
allele's robust *global* EC decrease stands, and its largest movers cluster
in the same C-terminal L helix (true 442-448: EC decreases at 442/444/445/448,
around the R443 CPR contact) and in SRS-4 — consistent with CYP2B6\*5
(R487C), a terminal β-strand variant (3IBD SHEET C pairs Arg487-Pro490 with
Phe457-Ala460, immediately adjacent to the L helix) with markedly reduced
hepatic expression, in part compensated by higher specific activity (Lang
et al. 2001, Pharmacogenetics 11:399-415; Zanger et al. 2007).

**Verdict.** K139E and R140Q's DRN findings survive the window-mean test and
are written up as genuine findings: both are delocalized (network-wide) but
structurally coherent, and each localizes to a literature-defined functional
surface (proximal-face CPR surface for K139E; active-site 5 Å set for
R140Q). R487C's *local* DRN finding does not survive and is reported as a
single-residue artifact (the same pattern as S259R's window-max SASA); only
its robust global EC finding is retained, with its C-terminal/L-helix
concentration noted as a secondary, non-robust observation.

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

**S259R SASA follow-up (site more exposed inside a more compact protein).**
The apparent contradiction in S259R -- robust global SASA decrease alongside
robust *local* SASA increase at its own site -- was checked residue-by-
residue (`s259r_sasa_followup.py`). The global decrease is confirmed
(d_avg −4.35 nm², rep1 −5.01, rep2 −3.69, WT noise floor 3.54; robust). The
local increase, however, is *strictly a single-residue effect at the
mutation position itself*: per-residue mean SASA at GROMACS resid 231 is
1.79/2.04 nm² in the two S259R replicates vs. 0.76/0.72 nm² in the two WT
replicates -- a ~2.5-fold increase (+1.18 nm²), the single largest
per-residue SASA change in the entire profile. Every residue in its ±3
window (228-230, 232-234) is *more buried* in S259R (deltas −0.08 to
−0.48 nm²). This explains the metric split: the framework's window-max
passes robustly (+0.432, rep deltas +0.393/+0.471, noise floor 0.173), while
the window-mean does not (rep1 +0.039 vs. rep2 −0.058, opposite signs). The
picture is not local loop unfolding -- it is specifically the mutated
sidechain flipping out into solvent while its immediate neighbors and the
whole protein pack tighter. Consistent secondary features: the most-buried
residue in S259R is GROMACS 111 (true 139, the K139E site / the 108-112
hotspot loop, −0.58 nm²), and the region 105-114 as a whole is markedly less
exposed, so S259R simultaneously opens its own site and compacts that
loop. Mildly *more* exposed secondary sites are GROMACS 107 (true 135),
216 (true 244), and 160-161 (true 188-189). This is a real, reproducible,
site-specific signature, not an artifact -- worth folding into the final
write-up as S259R's mechanistic note (a locally-remodeled, more exposed
surface patch at the mutation site inside an overall compacted protein).

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

**PCA and DCCM.** Unlike clustering (Part 5), which found no robust signal
for any allele, PCA and DCCM are the most productive of the "essential
dynamics" analyses run so far. **PC1 eigenvalue** (how much variance the
single dominant collective motion captures) shows a robust decrease
(more rigid dominant mode) in K139E, M46V, I391N, R140Q, and S259R, and a
robust increase in G99E, with T306S-R378K also robust (increase). **PC1
variance fraction** (how dominant that one mode is relative to all others)
shows a robust decrease in M46V, I391N, and R140Q specifically -- a
stronger, more specific finding than the eigenvalue metric alone, since a
lower fraction means motion is genuinely more evenly distributed across
multiple modes, not just globally smaller.

At the **DCCM global** level (mean |correlation| across the whole protein),
robust increases in overall coupling are seen in I328T, K262R, P428T, and
T306S-R378K, and a robust decrease in M46V and S259R. At the **DCCM local**
level (mutation site's long-range coupling to the rest of the protein,
excluding trivial backbone neighbors), only three alleles pass: I328T
(increase), K262R (decrease), and R487C (increase) -- notably, most alleles
with large raw local deltas (G99E, M46V, P428T, T306S-R378K) fail because
their replicates disagree in direction at that specific site, the same
failure pattern already seen repeatedly at the DRN layer.

Several convergences stand out. **M46V** now shows a sixth and seventh
independent robust finding: PC1 eigenvalue decrease, PC1 fraction decrease,
and DCCM global decrease -- all three pointing the same direction (reduced,
more evenly-distributed motion; lower overall residue-residue coupling),
reinforcing its position as the most convergent, clearly rigidifying/
compacting allele in the panel, though still without a robust *local* DCCM
signal at its own site (consistent with its DSSP finding being the only
one localized to the mutation site itself). **I328T** and **K262R** both
now show robust findings at BOTH the DCCM global and DCCM local level with
the same sign in each case (I328T: increase/increase; K262R: decrease/
decrease) -- the first time either of these two alleles shows a *local*,
mutation-site-specific network effect with this much internal consistency;
K262R in particular now has independent, direction-consistent evidence
across RMSF (local rigidification), DSSP (local ordered-SS increase), and
DCCM (local coupling decrease), a genuinely convergent multi-method local
story despite showing no DRN or H-bond signal at that same site. **R487C**
picks up its second robust local finding (DCCM, alongside its existing DRN
EC finding) -- both indicate increased local network engagement at its own
site, though see the residue-numbering caveat above regarding this
specific allele's site index. **R140Q** and **I391N** both add a robust
PC1-fraction decrease to their existing single- or dual-metric profiles.

**Binding-pocket (active-site) analysis.** The active-site volume and open
fraction are the most productive of the pocket metrics. At the **mean volume
over all snapshots** level (closed snapshots = 0; WT rep1-vs-rep2 noise floor
4.65 A^3), six alleles show robust increases: **R140Q** (delta = +98 A^3, the
strongest pocket-widening in the panel), **P428T** (+88 A^3), **I328T**
(+59 A^3), **R487C** (+31 A^3), **T306S-R378K** (+17 A^3), and **I391N**
(+8 A^3). The same six pass on **open fraction** (noise floor 0.128): R140Q
(+0.77, nearly always open), I328T (+0.48), R487C (+0.40), P428T (+0.38),
T306S-R378K (+0.22), I391N (+0.13). These are large, replicate-consistent
effects: R140Q's active site is open in ~92% of snapshots (vs. WT's 22%/9%
across its two replicates) at ~100-125 A^3 mean open volume (vs. WT's 40-45).
K139E and M46V show robust increases in **open-pocket volume only** (K139E
delta = +18.6 A^3, M46V +16.6 A^3) -- a larger pocket *when* it is open, but
no replicate-consistent change in how often it opens (their open-fraction
deltas are opposite-signed between replicates). G99E's open-volume increase
is directionally consistent but sits below its own noise floor (delta 3.1 vs.
5.0 A^3 floor) and is not called robust.

**Heme drift** (heme COM-to-protein COM distance, WT mean ~6.6-7.6 nm, noise
floor 1.01 nm) shows robust *decreases* in **K262R** (delta = -1.2 nm) and
**S259R** (delta = -1.5 nm) -- both mutants' heme sits significantly closer to
the protein center than WT on average, i.e. a more stably/centrally positioned
cofactor. No other allele passes (most agree in direction but fall inside the
noise floor). **Active-site RMSF** (mean over the WT pocket-lining residues,
noise floor 0.013 nm) adds only **I328T** (robust increase, +0.019 nm). The
**protein-heme H-bond contact** metric was not discriminative in this dataset:
the two WT replicates differ by ~42,000 frames (122,573 vs. 164,822 of 30,001)
-- a noise floor that dwarfs every allele's delta -- so no allele claim is
made on it; it is reported as a non-result rather than interpreted.

**Substrate-access-channel analysis (CAVER 3.0).** Two CAVER metrics survive
the WT-replicate-noise-floor framework. (1) **Open fraction** -- the fraction
of the 121 snapshots in which at least one water-passable (bottleneck >= 1.3
A) channel exists (WT noise floor 0.017). Robust *increases* are seen in
**I328T** (+0.45), **R487C** (+0.16), **K139E** (+0.15), **T306S-R378K**
(+0.13) and **M46V** (+0.10); **S259R** shows a robust *decrease* (-0.05).
Caveat: I328T's open fraction is sign-consistent between replicates but
magnitude-driven by a single replicate (0.917 vs. 0.091); its direction
(greater channel accessibility) is robust, its magnitude is not. (2) **Open
snapshot bottleneck radius** -- the mean widest bottleneck over open
snapshots only (WT noise floor just 0.003 A). Eight of eleven alleles show
robustly *wider* open channels: **K262R** (+0.113 A), **R140Q** (+0.120 A),
**T306S-R378K** (+0.067 A), **R487C** (+0.053 A), **K139E** (+0.038 A),
**M46V** (+0.035 A), **G99E** (+0.020 A) and **I391N** (+0.019 A). The mean
widest-bottleneck-over-all-snapshots metric (closed snapshots = 0) and the
>= 1.7 A open threshold do not support claims (WT noise floor too large /
zero). The open-fraction increases converge with the MDpocket active-site
widening above: I328T, R487C and T306S-R378K widen the pocket *and* open the
substrate-access channel more often, and S259R's channel closes more often,
consistent with its heme sitting more centrally and its RMSD/RMSF localization
at a single residue (below). K139E is notable in this context: a CPR-face
(proximal) DRN effect plus a distal-face channel that opens more often and
wider -- consistent with the allosteric-coupling interpretation of that
allele from the DRN deep-dive.

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
   rigidification even though the two methods aren't formally linked. DCCM
   now adds robust findings at BOTH the global (increase) and local
   (decrease) level, the local one directly at K262R's own mutation site --
   this is now the most internally consistent multi-method *local* case in
   the panel (RMSF rigidification, DSSP ordered-SS increase, DCCM local
   coupling decrease, all pointing toward a genuinely stiffer, more
   locally-isolated residue), despite no DRN or H-bond signal there.
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
   "purely allosteric/distributed" alone. PCA and DCCM now add a sixth and
   seventh: robust PC1 eigenvalue decrease, robust PC1 variance fraction
   decrease, and robust DCCM global decrease — all three global/distributed
   measures again, all pointing the same direction (less, more evenly-spread
   motion; lower overall coupling), with still no robust *local* DCCM
   signal at its own site. M46V is now the clearest, most convergent
   allosteric-type case in the panel, seven independent measures deep, with
   exactly one of those seven (DSSP) localized to the mutation site itself.
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
   Rg compaction. These were treated as provisional single-network-metric
   leads until the DRN deep-dive (see Results) resolved them: **K139E and
   R140Q's BC findings survive the window-mean test** and are written up as
   genuine — delocalized but structurally coherent, localizing to the
   proximal-face/L-helix CPR-binding surface (containing the mapped CPR
   contact R443) for K139E and to the active-site 5 Å lining set for R140Q
   (which also has the strongest pocket-widening in the panel). **R487C's
   local EC finding does not survive** (single-residue spike, same failure
   mode as S259R's window-max SASA); only its robust global EC decrease is
   retained.
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
     whole-protein SASA decrease. The targeted residue-level follow-up
     (above) resolves the pattern: the local increase is strictly the single
     mutated residue (GROMACS 231) flipping out to solvent (~2.5x WT
     per-residue exposure, +1.18 nm²), with its ±3 neighbors all slightly
     more buried and the whole protein more compact — a locally remodeled,
     more exposed surface patch at the mutation site inside an overall
     compacted protein, not local loop unfolding and not an artifact.
 10. **Binding pocket.** The pocket analysis adds a mechanistically
    meaningful layer that the earlier whole-protein metrics could not see:
    a robust, large active-site *widening/opening* in six alleles (R140Q,
    P428T, I328T, R487C, T306S-R378K, I391N). This is notable for **P428T
    and I328T** in particular — both were already the strongest cases for
    elevated global RMSD, local RMSF, and local H-bonds (Parts 1-2), and the
    pocket result corroborates that picture with an *independent, direct
    measure of active-site geometry*: a more open, larger active-site
    pocket is exactly what a more flexible loop/beta-sheet-adjacent
    mutation region around the substrate access route would be expected to
    produce. **R140Q**, previously only a single-DRN-metric lead, now shows
    the largest pocket widening in the panel (nearly-always-open active
    site at ~2.5x WT open volume) — its strongest evidence to date, and
    directionally consistent with a more accessible active site. **K262R
    and S259R** add the complementary result on the heme side: both show
    robustly decreased heme drift (heme pulled toward the protein center),
    consistent with (though not proof of) a more stably anchored cofactor —
    interesting for K262R in particular, whose robust local RMSF
    rigidification, DSSP ordered-SS increase, and DCCM local coupling
    decrease already suggested a stiffer, more locally-isolated residue.
    The heme H-bond metric is genuinely uninformative here (WT noise floor
    exceeds every allele's delta) and is reported as such. The CAVER 3.0
    substrate-access-channel round then adds a *distal* geometry measure to
    the pocket picture: the channel open fraction increases in I328T,
    R487C, K139E, T306S-R378K and M46V and decreases in S259R, and open
    channels are wider in 8/11 alleles -- so the mutations that widen and
    open the active-site pocket (I328T, R487C, T306S-R378K) also open its
    access channel more often. The whole-protein exploration frequency
    grids (Round 1) also completed for all 24 systems (non-empty grids,
    reported qualitatively as in the reference paper).

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
3. ~~Consider whether the single-method DRN-only findings (K139E, R140Q,
   R487C) merit a follow-up look at the actual network structure (e.g.
   which specific edges/neighbors drive the local centrality shift) rather
   than resting on the summary centrality value alone, before writing these
   up as findings in their own right~~ — done (Part 3 follow-up, see the DRN
   deep-dive in Results): `drn_network_deepdive.py`/`drn_network_deepdive.png`
   resolves K139E and R140Q into genuine, coherent BC neighborhood effects
   (K139E localized to the proximal-face/L-helix CPR surface, R140Q to the
   active-site 5 Å set) and R487C's local EC into a single-residue artifact
   (only its global EC is retained).
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
   RMSF rigidification at the same site).
7. ~~PCA and DCCM~~ — done (Part 7, above). Headline: unlike clustering,
   these are highly productive -- M46V picks up two more global convergent
   findings (sixth/seventh independent measure); K262R becomes the most
   internally-consistent *local* multi-method case in the panel (RMSF,
   DSSP, and DCCM local all agree at its own site); I328T gets its first
   robust DRN-layer-style local finding via DCCM (global+local both
   increase); R487C gets a second robust local finding.
8. ~~Binding-pocket (active-site) analysis~~ — done (Part 8, above).
    Headline: robust active-site volume/open-fraction *increases* in R140Q
    (the largest, nearly-always-open at ~2.5x WT open volume), P428T,
    I328T, R487C, T306S-R378K, I391N -- directly corroborating the Parts
    1-2 flexible-site story for P428T/I328T with an independent
    active-site-geometry measure; robust heme-drift *decreases* in K262R
    and S259R (heme pulled toward the protein center); substrate-access
    channel (CAVER 3.0) open-fraction *increases* in I328T, R487C, K139E,
    T306S-R378K, M46V and a decrease in S259R, with wider open channels in
    8/11 alleles; the Round 1 whole-protein exploration frequency grids
    also completed for all 24 systems (non-empty grids, viewable in VMD;
    reported qualitatively as in the reference paper).
8. ~~Double-check the R487C local-metric residue indexing in
   rg_sasa_significance_check.py~~ — done and fixed. `load_sasa_res` was
   silently dropping `gmx sasa -or`'s own resid column and assuming
   contiguous 1-based array indexing; the file's real resid column skips
   408 (heme, confirmed directly: `sasa_res_WT.xvg` jumps "407 ... 409"
   with no row for 408), so any site above 408 (only R487C, resid 459, in
   this project) was being read from the wrong array position -- a real
   off-by-one bug, not just a theoretical concern. Fixed by looking up by
   actual resid via a dict instead of by array position. Re-ran the full
   local SASA check after the fix: every allele's ROBUST/non-ROBUST
   conclusion is unchanged, including R487C's (still correctly
   non-robust) -- so no retraction is needed here, but the underlying code
   was wrong and is now fixed for future runs/re-analysis.
9. If more WT replicates become available at any point, re-run
   `cluster_summary.py`'s robustness check -- a 2-replicate noise floor is
   an especially weak substitute specifically for the clustering metric
   (see above), more so than for the scalar/per-residue metrics elsewhere
   in this project.
 10. ~~S259R's local-vs-global SASA divergence (site more exposed while the
     whole protein is more compact)~~ — done and resolved (see the S259R
     follow-up subsection in Results): `s259r_sasa_followup.py` shows the
     local increase is strictly the single mutated residue (GROMACS 231,
     ~2.5x WT exposure, +1.18 nm²) flipping out to solvent while its ±3
     neighbors and the whole protein are all more buried; the framework's
     window-max metric passes robustly while window-mean does not, because
     the effect is a single-residue outlier, not local loop unfolding.
     Figure: `s259r_sasa_followup.png`; the dashboard's S259R panel now
     plots the real per-residue SASA delta profile.
11. Request or approximate a proper reference-triplicate 3-SD threshold (the
    CYP3A4 paper's framework) rather than continuing to rely on a 2-replicate
    noise floor, which remains a workable but weaker substitute — this now
    applies across all eight analyses completed so far.
12. Decide on final reporting format/figures for all eight completed parts
    (RMSD/RMSF, H-bonds, DRN, Rg/SASA, clustering, DSSP, PCA/DCCM,
    binding-pocket) for the supervision meeting write-up.
