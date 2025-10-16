#!/usr/bin/env python3
"""
Lineage Tracing Upstream Processing Pipeline

This script implements the upstream processing pipeline for CRISPR/Cas9-based 
lineage tracing data using Cassiopeia. The pipeline processes raw FASTQ files
into character matrices suitable for phylogenetic reconstruction.

Based on the Cassiopeia preprocessing workflow described in:
"Inference of Single-Cell Phylogenies from Lineage Tracing Data"
"""

import pandas as pd
import cassiopeia as cas
import os
import argparse


def setup_preprocessing_pipeline(input_files, sample_name, output_directory, 
                                reference_filepath, n_threads=8, 
                                allow_allele_conflicts=False, verbose=True):
    """
    Complete upstream processing pipeline for lineage tracing data
    
    Parameters:
    -----------
    input_files : list
        List of FASTQ file paths [R1.fastq.gz, R2.fastq.gz]
    sample_name : str
        Sample identifier
    output_directory : str
        Output directory path
    reference_filepath : str
        Path to reference FASTA file
    n_threads : int
        Number of threads for processing
    allow_allele_conflicts : bool
        Whether to allow allele conflicts
    verbose : bool
        Verbose output
    """
    
    # Setup output directory
    cas.pp.setup(output_directory, verbose=verbose)
    
    print(f"Processing sample: {sample_name}")
    print(f"Input files: {input_files}")
    print(f"Output directory: {output_directory}")
    
    # Step 1: Convert FASTQs to unmapped BAM
    print("Step 1: Converting FASTQs to unmapped BAM...")
    bam_fp = cas.pp.convert_fastqs_to_unmapped_bam(
        input_files,
        chemistry='10xv3',
        output_directory=output_directory,
        name=sample_name,
        n_threads=n_threads
    )
    
    # Step 2: Filter BAM by quality
    print("Step 2: Filtering BAM by quality...")
    bam_fp = cas.pp.filter_bam(
        bam_fp,
        output_directory=output_directory,
        quality_threshold=10,
        n_threads=n_threads,
    )
    
    # Step 3: Collapse UMIs (optional - can be skipped)
    print("Step 3: Collapsing UMIs...")
    umi_table = cas.pp.collapse_umis(
        bam_fp,
        output_directory=output_directory,
        max_hq_mismatches=3,
        max_indels=2,
        method='likelihood',
        n_threads=n_threads,
    )
    
    # Step 4: Resolve UMI sequences
    print("Step 4: Resolving UMI sequences...")
    umi_table = cas.pp.resolve_umi_sequence(
        umi_table,
        output_directory=output_directory,
        min_umi_per_cell=10,
        min_avg_reads_per_umi=2.0,
        plot=True,
    )
    
    # Step 5: Align sequences to reference
    print("Step 5: Aligning sequences to reference...")
    umi_table = cas.pp.align_sequences(
        umi_table,
        ref_filepath=reference_filepath,
        gap_open_penalty=20,
        gap_extend_penalty=1,
        n_threads=n_threads,
    )
    
    # Step 6: Call alleles
    print("Step 6: Calling alleles...")
    umi_table = cas.pp.call_alleles(
        umi_table,
        ref_filepath=reference_filepath,
        barcode_interval=(21, 35),  # Adjust based on your design
        cutsite_locations=[113, 167, 221],  # Adjust based on your design
        cutsite_width=12,
        context=True,
        context_size=5,
    )
    
    # Step 7: Error correct UMIs
    print("Step 7: Error correcting UMIs...")
    umi_table = cas.pp.error_correct_umis(
        umi_table,
        max_umi_distance=2,
        allow_allele_conflicts=allow_allele_conflicts,
        n_threads=n_threads,
    )
    
    # Step 8: Filter molecule table
    print("Step 8: Filtering molecule table...")
    umi_table = cas.pp.filter_molecule_table(
        umi_table,
        output_directory=output_directory,
        min_umi_per_cell=10,
        min_avg_reads_per_umi=2.0,
        min_reads_per_umi=-1,
        intbc_prop_thresh=0.5,
        intbc_umi_thresh=10,
        intbc_dist_thresh=1,
        doublet_threshold=0.35,
        allow_allele_conflicts=allow_allele_conflicts,
        plot=True,
    )
    
    # Step 9: Call lineage groups
    print("Step 9: Calling lineage groups...")
    allele_table = cas.pp.call_lineage_groups(
        umi_table,
        output_directory=output_directory,
        min_umi_per_cell=10,
        min_avg_reads_per_umi=2.0,
        min_cluster_prop=0.005,
        min_intbc_thresh=0.05,
        inter_doublet_threshold=0.35,
        kinship_thresh=0.25,
        plot=True,
    )
    
    # Save final allele table
    output_file = os.path.join(output_directory, f"{sample_name}.alleleTable.csv")
    allele_table.to_csv(output_file, sep="\t", index=False)
    print(f"Final allele table saved to: {output_file}")
    
    return allele_table


def main():
    """Main function for command line execution"""
    parser = argparse.ArgumentParser(description='Lineage Tracing Upstream Processing Pipeline')
    parser.add_argument('--input_files', nargs=2, required=True, 
                       help='Input FASTQ files (R1 R2)')
    parser.add_argument('--sample_name', required=True, 
                       help='Sample name')
    parser.add_argument('--output_directory', required=True, 
                       help='Output directory')
    parser.add_argument('--reference_filepath', required=True, 
                       help='Reference FASTA file path')
    parser.add_argument('--n_threads', type=int, default=8, 
                       help='Number of threads')
    parser.add_argument('--allow_allele_conflicts', action='store_true', 
                       help='Allow allele conflicts')
    parser.add_argument('--verbose', action='store_true', 
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Run the pipeline
    allele_table = setup_preprocessing_pipeline(
        input_files=args.input_files,
        sample_name=args.sample_name,
        output_directory=args.output_directory,
        reference_filepath=args.reference_filepath,
        n_threads=args.n_threads,
        allow_allele_conflicts=args.allow_allele_conflicts,
        verbose=args.verbose
    )
    
    print("Pipeline completed successfully!")


if __name__ == "__main__":
    main()