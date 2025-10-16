# Lineage Tracing Upstream Processing

This directory contains the upstream processing pipeline for CRISPR/Cas9-based lineage tracing data using Cassiopeia.

## Overview

The upstream processing pipeline converts raw FASTQ sequencing files into character matrices suitable for phylogenetic reconstruction. This pipeline is based on the Cassiopeia preprocessing workflow.

## Pipeline Steps

The complete upstream processing pipeline includes the following steps:

### 1. Data Conversion
- **convert_fastqs_to_unmapped_bam**: Convert FASTQs into unmapped BAM format while parsing barcode and UMI sequences

### 2. Quality Control
- **filter_bam**: Filter reads with low-quality barcode and/or UMI sequences
- **error_correct_cellbcs_to_whitelist**: Correct sequencing errors using predefined barcode whitelist (optional)

### 3. UMI Processing
- **collapse_umis**: Collapse reads into UMIs by constructing consensus sequences
- **resolve_umi_sequence**: Resolve a single sequence for each UMI by choosing the most likely read

### 4. Sequence Alignment
- **align_sequences**: Align sequences to reference target site using Smith-Waterman algorithm

### 5. Allele Calling
- **call_alleles**: Call alleles with respect to reference target site and report mutations

### 6. Error Correction
- **error_correct_intbcs_to_whitelist**: Correct intBC sequencing errors (if applicable)
- **error_correct_umis**: Error-correct UMIs with identical mutation data and similar sequences

### 7. Quality Filtering
- **filter_molecule_table**: Filter UMIs with conflicting allele information or insufficient reads

### 8. Lineage Grouping
- **call_lineage_groups**: Group cells into lineage clusters based on shared mutations

## Usage

### Command Line Usage

```bash
python lineage_preprocessing_pipeline.py \
    --input_files sample_R1.fastq.gz sample_R2.fastq.gz \
    --sample_name SAMPLE_NAME \
    --output_directory /path/to/output \
    --reference_filepath /path/to/reference.fa \
    --n_threads 8 \
    --verbose
```

## Key Parameters

### Barcode and Cutsite Configuration
- **barcode_interval**: Position of integration barcode in sequence (e.g., (21, 35))
- **cutsite_locations**: Positions of Cas9 cut sites (e.g., [113, 167, 221])
- **cutsite_width**: Width around cut sites to consider for mutations (default: 12)

### Quality Control Thresholds
- **min_umi_per_cell**: Minimum UMIs required per cell (default: 10)
- **min_avg_reads_per_umi**: Minimum average reads per UMI (default: 2.0)
- **quality_threshold**: Minimum quality score for reads (default: 10)

## Output Files

The pipeline generates several intermediate and final output files:

- **{sample_name}.unmapped.bam**: Unmapped BAM file with parsed barcodes
- **{sample_name}.filtered.bam**: Quality-filtered BAM file
- **{sample_name}.collapsed.txt**: UMI-collapsed molecule table
- **{sample_name}.resolved.txt**: Resolved UMI sequences
- **{sample_name}.aligned.txt**: Aligned sequences
- **{sample_name}.alleles.txt**: Called alleles
- **{sample_name}.filtered_molecules.txt**: Filtered molecule table
- **{sample_name}.alleleTable.csv**: Final allele table for downstream analysis

## Dependencies

- cassiopeia
- pandas
- numpy
- matplotlib (for plotting)
- seaborn (for plotting)

## Two Methods for Running Preprocessing

Based on the M-AP sample analysis, there are two main approaches to run the preprocessing pipeline:

### Method 1: Using Configuration File (Recommended)

This method uses Cassiopeia's built-in configuration system with a `.cfg` file:

```bash
# Run using configuration file
cassiopeia-preprocess preprocess_template.cfg
```

**Advantages:**
- Easy to configure and modify parameters
- Built-in parameter validation
- Automatic pipeline execution
- Less prone to errors

**Configuration Template:** `preprocess_template.cfg`

### Method 2: Using Python Script

This method provides more control over the preprocessing steps using Python API:

```bash
# Run using Python script
python run_template.py --sample_name SAMPLE_NAME --input_dir /path/to/input --output_dir /path/to/output --reference_filepath /path/to/reference.fa
```

**Advantages:**
- More flexibility and control
- Can customize intermediate steps
- Better for debugging
- Can integrate with other Python workflows

**Script Template:** `run_template.py`

## Sample-Specific Configurations

### M-AP Sample Configuration
```
barcode_interval = (20, 34)
cutsite_locations = [112, 166, 220]
```

### Pre-MGH Sample Configuration  
```
barcode_interval = (21, 35)
cutsite_locations = [107, 161, 215]
```

## Typical Workflow

1. **Prepare input files**: Ensure FASTQ files are available
2. **Choose method**: Configuration file or Python script
3. **Adjust parameters**: Modify barcode intervals and cutsite locations based on sample type
4. **Run preprocessing**: Execute the chosen method
5. **Check outputs**: Verify intermediate and final results

## Output Structure

Both methods generate the same output files:
- `{sample}.align.txt`: Aligned sequences
- `{sample}.call_alleles.txt`: Called alleles
- `{sample}.error_correct_umis.txt`: Error-corrected UMIs
- `{sample}.filter_molecule_table.txt`: Filtered molecules
- `{sample}.alleleTable.csv`: Final allele table for downstream analysis

