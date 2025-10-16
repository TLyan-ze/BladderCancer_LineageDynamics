#!/usr/bin/env python3
"""
Example Configuration for Lineage Tracing Upstream Processing

This file contains example configurations for different sample types
used in the bladder cancer lineage dynamics study.
"""

# Sample configurations based on temp_up analysis
SAMPLE_CONFIGS = {
    "Met-AP": {
        "barcode_interval": (20, 34),
        "cutsite_locations": [112, 166, 220],
        "cutsite_width": 12,
        "chemistry": "10xv3"
    },
    
    "Pre-MGH": {
        "barcode_interval": (21, 35),
        "cutsite_locations": [107, 161, 215],
        "cutsite_width": 12,
        "chemistry": "10xv3"
    },
    
    "Default": {
        "barcode_interval": (21, 35),
        "cutsite_locations": [113, 167, 221],
        "cutsite_width": 12,
        "chemistry": "10xv3"
    }
}

# Quality control parameters
QC_PARAMS = {
    "quality_threshold": 10,
    "min_umi_per_cell": 10,
    "min_avg_reads_per_umi": 2.0,
    "min_reads_per_umi": -1,
    "intbc_prop_thresh": 0.5,
    "intbc_umi_thresh": 10,
    "intbc_dist_thresh": 1,
    "doublet_threshold": 0.35,
    "kinship_thresh": 0.25
}

# UMI processing parameters
UMI_PARAMS = {
    "max_hq_mismatches": 3,
    "max_indels": 2,
    "method": "likelihood",
    "max_umi_distance": 2
}

# Alignment parameters
ALIGNMENT_PARAMS = {
    "gap_open_penalty": 20,
    "gap_extend_penalty": 1
}

# Lineage grouping parameters
LINEAGE_PARAMS = {
    "min_cluster_prop": 0.005,
    "min_intbc_thresh": 0.05,
    "inter_doublet_threshold": 0.35
}


def get_sample_config(sample_name):
    """
    Get configuration for a specific sample
    
    Parameters:
    -----------
    sample_name : str
        Name of the sample
        
    Returns:
    --------
    dict : Configuration dictionary
    """
    # Extract sample type from sample name
    if "Met-AP" in sample_name:
        return SAMPLE_CONFIGS["Met-AP"]
    elif "Pre-MGH" in sample_name:
        return SAMPLE_CONFIGS["Pre-MGH"]
    else:
        return SAMPLE_CONFIGS["Default"]


def get_all_params():
    """
    Get all parameter dictionaries
    
    Returns:
    --------
    dict : Dictionary containing all parameter sets
    """
    return {
        "qc_params": QC_PARAMS,
        "umi_params": UMI_PARAMS,
        "alignment_params": ALIGNMENT_PARAMS,
        "lineage_params": LINEAGE_PARAMS
    }


# Example usage
if __name__ == "__main__":
    # Example for Met-AP sample
    sample_name = "Met-AP_L3_S0006B0006"
    config = get_sample_config(sample_name)
    print(f"Configuration for {sample_name}:")
    print(f"Barcode interval: {config['barcode_interval']}")
    print(f"Cutsite locations: {config['cutsite_locations']}")
    
    # Get all parameters
    all_params = get_all_params()
    print(f"\nQuality control parameters: {all_params['qc_params']}")