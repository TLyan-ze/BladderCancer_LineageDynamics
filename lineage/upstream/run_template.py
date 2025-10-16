#!/usr/bin/env python3
"""
Lineage Tracing Preprocessing Template Script

This script provides a template for running Cassiopeia preprocessing
using Python API instead of configuration files. Based on M-AP sample processing.

Usage:
    python run_template.py --sample_name SAMPLE_NAME --input_dir /path/to/input --output_dir /path/to/output
"""

import pandas as pd
import cassiopeia as cas
import os
import argparse


def run_preprocessing(sample_name, input_dir, output_dir, reference_filepath, 
                     barcode_interval=(20, 34), cutsite_locations=[112, 166, 220],
                     n_threads=32, allow_allele_conflicts=False, verbose=True):
    """
    Run Cassiopeia preprocessing using Python API
    
    Parameters:
    -----------
    sample_name : str
        Sample identifier
    input_dir : str
        Directory containing input FASTQ files
    output_dir : str
        Output directory
    reference_filepath : str
        Path to reference FASTA file
    barcode_interval : tuple
        Barcode interval positions
    cutsite_locations : list
        Cutsite locations
    n_threads : int
        Number of threads
    allow_allele_conflicts : bool
        Whether to allow allele conflicts
    verbose : bool
        Verbose output
    """
    
    # Setup output directory
    cas.pp.setup(output_dir, verbose=verbose)
    
    print(f"Processing sample: {sample_name}")
    print(f"Output directory: {output_dir}")
    print(f"Reference file: {reference_filepath}")
    
    # Load aligned sequences (assuming alignment step is already done)
    align_file = os.path.join(output_dir, f"{sample_name}.align.txt")
    
    if os.path.exists(align_file):
        print(f"Loading aligned sequences from: {align_file}")
        umi_table = pd.read_csv(align_file, sep="\t")
        print(f"Loaded {len(umi_table)} sequences")
    else:
        print(f"Aligned file not found: {align_file}")
        print("Please run alignment step first or provide aligned sequences")
        return None
    
    # Step 1: Call alleles
    print("Step 1: Calling alleles...")
    umi_table = cas.pp.call_alleles(
        umi_table,
        ref_filepath=reference_filepath,
        barcode_interval=barcode_interval,
        cutsite_locations=cutsite_locations,
        cutsite_width=12,
        context=True,
        context_size=5,
    )
    
    # Save intermediate result
    call_alleles_file = os.path.join(output_dir, f"{sample_name}.call_alleles.txt")
    umi_table.to_csv(call_alleles_file, sep="\t", index=False)
    print(f"Saved call_alleles result to: {call_alleles_file}")
    
    # Step 2: Error correct UMIs
    print("Step 2: Error correcting UMIs...")
    umi_table = cas.pp.error_correct_umis(
        umi_table,
        max_umi_distance=2,
        allow_allele_conflicts=allow_allele_conflicts,
        n_threads=n_threads,
    )
    
    # Save intermediate result
    error_correct_file = os.path.join(output_dir, f"{sample_name}.error_correct_umis.txt")
    umi_table.to_csv(error_correct_file, sep="\t", index=False)
    print(f"Saved error_correct_umis result to: {error_correct_file}")
    
    # Step 3: Filter molecule table
    print("Step 3: Filtering molecule table...")
    umi_table = cas.pp.filter_molecule_table(
        umi_table,
        output_directory=output_dir,
        min_umi_per_cell=10,
        min_avg_reads_per_umi=2.0,
        min_reads_per_umi=-1,
        intbc_prop_thresh=0.5,
        intbc_umi_thresh=10,
        intbc_dist_thresh=1,
        doublet_threshold=0.35,
        allow_allele_conflicts=allow_allele_conflicts,
        plot=False,
    )
    
    # Save intermediate result
    filter_file = os.path.join(output_dir, f"{sample_name}.filter_molecule_table.txt")
    umi_table.to_csv(filter_file, sep="\t", index=False)
    print(f"Saved filter_molecule_table result to: {filter_file}")
    
    # Step 4: Call lineage groups
    print("Step 4: Calling lineage groups...")
    allele_table = cas.pp.call_lineage_groups(
        umi_table,
        output_directory=output_dir,
        min_umi_per_cell=10,
        min_avg_reads_per_umi=2.0,
        min_cluster_prop=0.005,
        min_intbc_thresh=0.05,
        inter_doublet_threshold=0.35,
        kinship_thresh=0.25,
        plot=True,
    )
    
    # Save final result
    final_file = os.path.join(output_dir, f"{sample_name}.alleleTable.csv")
    allele_table.to_csv(final_file, sep="\t", index=False)
    print(f"Saved final allele table to: {final_file}")
    
    print("Preprocessing completed successfully!")
    return allele_table


def main():
    """Main function for command line execution"""
    parser = argparse.ArgumentParser(description='Lineage Tracing Preprocessing Template')
    parser.add_argument('--sample_name', required=True, help='Sample name')
    parser.add_argument('--input_dir', required=True, help='Input directory')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--reference_filepath', required=True, help='Reference FASTA file')
    parser.add_argument('--barcode_interval', nargs=2, type=int, default=[20, 34],
                       help='Barcode interval (start end)')
    parser.add_argument('--cutsite_locations', nargs='+', type=int, default=[112, 166, 220],
                       help='Cutsite locations')
    parser.add_argument('--n_threads', type=int, default=32, help='Number of threads')
    parser.add_argument('--allow_allele_conflicts', action='store_true',
                       help='Allow allele conflicts')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    # Convert barcode_interval to tuple
    barcode_interval = tuple(args.barcode_interval)
    
    # Run preprocessing
    allele_table = run_preprocessing(
        sample_name=args.sample_name,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        reference_filepath=args.reference_filepath,
        barcode_interval=barcode_interval,
        cutsite_locations=args.cutsite_locations,
        n_threads=args.n_threads,
        allow_allele_conflicts=args.allow_allele_conflicts,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()


# Example usage for different samples:
"""
# For M-AP sample:
python run_template.py \
    --sample_name M-AP \
    --input_dir /path/to/input \
    --output_dir /path/to/output \
    --reference_filepath /path/to/PCT48.ref.fa \
    --barcode_interval 20 34 \
    --cutsite_locations 112 166 220 \
    --n_threads 32 \
    --verbose

# For Pre-MGH sample:
python run_template.py \
    --sample_name Pre-MGH \
    --input_dir /path/to/input \
    --output_dir /path/to/output \
    --reference_filepath /path/to/PCT48.ref.fa \
    --barcode_interval 21 35 \
    --cutsite_locations 107 161 215 \
    --n_threads 32 \
    --verbose
"""