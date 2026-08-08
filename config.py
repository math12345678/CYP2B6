#!/usr/bin/env python3
"""
Centralized configuration module for CYP2B6 allele analysis pipeline.

This module consolidates all allele definitions, file paths, residue mappings,
and analysis parameters used across the 15+ analysis scripts in the
CYP2B6 MD analysis pipeline, eliminating duplication and standardizing
configuration.

Usage:
    from config import (
        ALLELES,
        ALLELES_SIMPLE,
        ALLELES_MUTATION_SITES,
        PCA_ALLELES,
        CLUSTERING_ALLELES,
        SYSTEMS_ALL,
        WT_REPLICATES,
        RESID_EXCLUDE,
        BASE_DIRECTORY,
        MIN_EFFECT_NM
    )

Run from ~/Desktop/Research/Research_Projects/RU-CYP2B6.
"""

# Base directory for the project (can be overridden)
BASE_DIRECTORY = "."

# All 11 CYP2B6 alleles being analyzed (mutant + wild-type)
ALLELES = [
    "G99E",
    "K139E", 
    "M46V",
    "I328T",
    "I391N",
    "K262R",
    "R140Q",
    "R487C",
    "P428T",
    "S259R",
    "T306S-R378K"
]

# Simplified allele list for scripts that only need the allele names
ALLELES_SIMPLE = ALLELES.copy()

# Alleles with mutation site information for different analysis types
# Format: {allele_name: (rep1_dir, rep2_dir, mutation_site_residue, true_mutation_description)}
ALLELES_MUTATION_SITES = {
    # For RMSF mutation-site comparisons (+/-3 residue window max)
    "G99E": ("G99E", "G99E_2", 71, "G99E (GROMACS 71 = true residue 99)"),
    "K139E": ("K139E", "K139E_2", 111, "K139E (GROMACS 111 = true residue 139)"),
    "M46V": ("M46V", "M46V_2", 18, "M46V (GROMACS 18 = true residue 46)"),
    "I328T": ("I328T", "I328T_2", 300, "I328T (GROMACS 300 = true residue 328)"),
    "I391N": ("I391N", "I391N_2", 363, "I391N (GROMACS 363 = true residue 391)"),
    "K262R": ("K262R", "K262R_2", 234, "K262R (GROMACS 234 = true residue 262)"),
    "R140Q": ("R140Q", "R140Q_2", 112, "R140Q (GROMACS 112 = true residue 140)"),
    "R487C": ("R487C", "R487C_2", 459, "R487C (GROMACS 459 = true residue 487)"),
    "P428T": ("P428T", "P428T_2", 400, "P428T (GROMACS 400 = true residue 428)"),
    "S259R": ("S259R", "S259R_2", 231, "S259R (GROMACS 231 = true residue 259)"),
    "T306S-R378K": ("T306S-R378K", "T306S-R378K_2", 278, "T306 (GROMACS 278 = true residue 306)"),
}

# Simple allele-to-site mapping for scripts that only need residue numbers
# GROMACS numbering (true residue = GROMACS residue + 28)
ALLELES_MUTATION_SITES_SIMPLE = {
    name: (dirs[0], dirs[1], dirs[2])
    for name, dirs in ALLELES_MUTATION_SITES.items()
}

# For PCA analysis - same structure but PCA only uses global metrics
PCA_ALLELES = ALLELES_MUTATION_SITES.copy()

# For clustering analysis - simplified, only allele names
CLUSTERING_ALLELES = ["G99E", "K139E", "M46V", "I328T", "I391N", "K262R",
                      "R140Q", "R487C", "P428T", "S259R", "T306S-R378K"]

# All systems for batch processing (WT + 2 replicates each of 11 alleles)
SYSTEMS_ALL = [
    "WT", "WT_2",
    "G99E", "G99E_2", "K139E", "K139E_2", "M46V", "M46V_2",
    "I328T", "I328T_2", "I391N", "I391N_2", "K262R", "K262R_2",
    "R140Q", "R140Q_2", "R487C", "R487C_2", "P428T", "P428T_2",
    "S259R", "S259R_2", "T306S-R378K", "T306S-R378K_2"
]

# WT replicates for statistical robustness checks
WT_REPLICATES = {
    "WT": "WT",
    "WT_2": "WT_2"
}

# Residues to exclude from per-residue analyses (contaminating heme component)
# Row/resid 408 = heme cofactor's CM1 component in combined numbering
RESID_EXCLUDE = {408}

# Minimum effect sizes for robustness checks (from significance_check.py)
MIN_EFFECT_NM = 0.01  # nm, for RMSD/RMSF/DRN/DSSM/H-bond significance

# Analysis-specific constants
DRN_WINDOW_SIZE = 3  # +/-3 residues for local DRN analysis
RMSF_WINDOW_SIZE = 3  # +/-3 residues for mutation-site RMSF analysis
SASA_WINDOW_SIZE = 3  # +/-3 residues for mutation-site SASA analysis
DSSM_WINDOW_SIZE = 3  # +/-3 residues for mutation-site DSSM analysis

# Metric-specific minimum effect sizes (from pocket_significance_check.py)
POCKET_METRIC_MIN_EFFECTS = {
    "mean_vol_all": 2.0,      # A^3
    "mean_vol_open": 3.0,     # A^3
    "open_frac": 0.03,        # fraction
    "heme_hbond_sum": 3.0,    # frames (of 30001)
    "heme_drift_mean": 0.01,  # nm
    "active_site_rmsf": 0.01, # nm
}

# File patterns and naming conventions
FILE_PATTERNS = {
    "rmsd": "{system}/rmsd_{system}.xvg",
    "rmsf": "{system}/rmsf_{system}.xvg",
    "drn_mean": "{system}/md_noWAT_mean.csv",
    "drn_cif": "{system}/md_noWAT_mean.cif",
    "dssp": "{system}/dssp_{system}.dat",
    "dssp_summary": "{system}/dssp_{system}_orderedss.csv",
    "cluster_summary": "{system}/cluster_{system}_summary.dat",
    "eigenval": "{system}/eigenval_{system}.xvg",
    "covar": "{system}/covar_{system}.dat",
    "pocket_descriptors": "{system}/mdpock_{system}_descriptors.txt",
    "pocket_summary": "{system}/pocket_{system}.csv",
    "hbond_pairs": "{system}/hbond_pairs_{system}.csv",
}

# Analysis step ordering (for dependency management)
ANALYSIS_STEPS = [
    "fix_topology",
    "run_all_hbonds",
    "run_all_drn",
    "run_all_rg_sasa",
    "run_all_clustering",
    "run_all_dssp_gmx",
    "run_all_pca_dccm",
    "run_all_mdpocket",
    "run_all_exploration",
]

# Utility functions for file path resolution
def get_file_path(pattern, system=None, **kwargs):
    """
    Get file path for a given pattern and system.
    
    Args:
        pattern: File pattern from FILE_PATTERNS
        system: System name (if pattern requires it)
        **kwargs: Additional format parameters
    
    Returns:
        str: Formatted file path
    """
    if system:
        return pattern.format(system=system, **kwargs)
    return pattern.format(**kwargs)

def get_allele_files(allele, data_type="rmsf"):
    """
    Get file paths for both replicates of an allele for a specific data type.
    
    Args:
        allele: Allele name
        data_type: Type of data ("rmsd", "rmsf", "drn", "pocket_descriptors")
    
    Returns:
        tuple: (rep1_path, rep2_path)
    
    Raises:
        KeyError: If allele not found or data_type not supported
    """
    if allele not in ALLELES_MUTATION_SITES_SIMPLE:
        raise KeyError(f"Unknown allele: {allele}")
    
    rep1_dir, rep2_dir, _ = ALLELES_MUTATION_SITES_SIMPLE[allele]
    
    pattern_map = {
        "rmsd": f"{rep1_dir}/rmsd_{rep1_dir}.xvg",
        "rmsf": f"{rep1_dir}/rmsf_{rep1_dir}.xvg",
        "drn": f"{rep1_dir}/md_noWAT_mean.csv",
        "pocket_descriptors": f"{rep1_dir}/mdpock_{rep1_dir}_descriptors.txt",
    }
    
    if data_type not in pattern_map:
        raise ValueError(f"Unsupported data_type: {data_type}. "
                        f"Supported: {list(pattern_map.keys())}")
    
    rep1_path = pattern_map[data_type]
    rep2_path = rep1_path.replace(f"{rep1_dir}", f"{rep2_dir}")
    
    return rep1_path, rep2_path

def get_mutation_site(allele):
    """
    Get the mutation site (GROMACS residue number) for an allele.
    
    Args:
        allele: Allele name
    
    Returns:
        int: GROMACS residue number where mutation occurs
    
    Raises:
        KeyError: If allele not found
    """
    if allele not in ALLELES_MUTATION_SITES_SIMPLE:
        raise KeyError(f"Unknown allele: {allele}")
    
    return ALLELES_MUTATION_SITES_SIMPLE[allele][2]

# Print configuration summary for debugging
def print_config_summary():
    """Print a summary of the configuration for debugging."""
    print("=" * 70)
    print("CYP2B6 Analysis Pipeline Configuration Summary")
    print("=" * 70)
    print(f"Total alleles: {len(ALLELES)}")
    print(f"Alleles: {', '.join(ALLELES)}")
    print(f"\nAll systems: {len(SYSTEMS_ALL)}")
    print(f"  WT replicates: WT, WT_2")
    for allele in ALLELES:
        rep1, rep2, site = ALLELES_MUTATION_SITES_SIMPLE[allele]
        true_residue = site + 28
        print(f"  {allele:12s}: {rep1:15s}, {rep2:15s} (site {site}, true {true_residue})")
    print(f"\nExcluded residues: {RESID_EXCLUDE}")
    print(f"Minimum effect size (RMSD/RMSF): {MIN_EFFECT_NM} nm")
    print("=" * 70)

if __name__ == "__main__":
    # When run as a script, print configuration summary
    print_config_summary()