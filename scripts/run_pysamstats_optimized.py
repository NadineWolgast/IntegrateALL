#!/usr/bin/env python3
"""
Optimized pysamstats processing for hotspot mutations.
Processes all mutations in batch for better performance.
"""

import pysam
import pysamstats
import sys
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from functools import partial

# Path to hotspots configuration file (relative to script location)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HOTSPOTS_CONFIG = os.path.join(os.path.dirname(SCRIPT_DIR), "config", "hotspots_config.csv")


def load_hotspots_config(config_file=HOTSPOTS_CONFIG):
    """
    Load hotspots configuration from CSV file.

    Returns:
        list: List of tuples (gene, hotspot, chromosome, start, end, strand)
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Hotspots config file not found: {config_file}")

    df = pd.read_csv(config_file)

    # Convert DataFrame to list of tuples for compatibility with existing code
    mutations_data = []
    for _, row in df.iterrows():
        mutations_data.append((
            row['gene'],
            row['hotspot'],
            str(row['chromosome']),  # Ensure chromosome is string
            int(row['start']),
            int(row['end']),
            row['strand']
        ))

    return mutations_data


def process_single_mutation(mutation_data, bam_file, fasta_file, sample_id, output_dir):
    """Process a single mutation hotspot with strand-specific handling."""
    gene, hotspot, chromosome, start, end, strand = mutation_data

    output_file = os.path.join(output_dir, f"{sample_id}_{gene}_{hotspot}.tsv")

    try:
        # CRITICAL FIX: Convert 1-based coordinates (config) to 0-based (pysamstats/BAM)
        # Config format: start=1-based inclusive, end=1-based exclusive (like VCF ranges)
        # pysamstats expects: start=0-based inclusive, end=0-based exclusive (like Python range)
        # Example NRAS G13 codon (1-based): 114716122, 114716123, 114716124
        #   Config: start=114716122, end=114716125 (1-based, end exclusive)
        #   pysamstats: start=114716121, end=114716124 (0-based, end exclusive)
        start_0based = start - 1  # Convert 1-based inclusive to 0-based inclusive
        end_0based = end - 1      # Convert 1-based exclusive to 0-based exclusive

        # Open files once per mutation (still more efficient than before)
        with pysam.AlignmentFile(bam_file) as mybam, pysam.Fastafile(fasta_file) as myfasta:
            
            # Strand-specific parameters for pysamstats
            # Plus strand genes (like IKZF1) need different handling than minus strand genes
            if strand == "+":
                # Plus strand (IKZF1) - originally used -u flag
                # -u flag typically means "unique mapping" or special read handling
                variations = pysamstats.load_variation(
                    mybam, myfasta,
                    chrom=chromosome,
                    start=start_0based,
                    end=end_0based,
                    truncate=True,
                    max_depth=8000,  # Lower depth threshold for plus strand
                    min_mapq=1,      # Minimum mapping quality (equivalent to -u behavior)
                    min_baseq=13     # Minimum base quality for plus strand
                )
            else:  # strand == "-"
                # Minus strand genes (majority of our genes)
                variations = pysamstats.load_variation(
                    mybam, myfasta,
                    chrom=chromosome,
                    start=start_0based,
                    end=end_0based,
                    truncate=True,
                    max_depth=1000000,  # Higher depth for minus strand
                    min_mapq=0,         # More permissive mapping quality
                    min_baseq=13        # Standard base quality
                )
            
            # Write output efficiently
            with open(output_file, 'w') as f:
                # Write header
                f.write("chrom\tpos\tref\treads_all\treads_pp\tmatches\tmatches_pp\t"
                       "mismatches\tmismatches_pp\tdeletions\tdeletions_pp\t"
                       "insertions\tinsertions_pp\tA\tA_pp\tC\tC_pp\tT\tT_pp\t"
                       "G\tG_pp\tN\tN_pp\n")

                # Write data rows with 1-based position conversion
                for entry in variations:
                    # Convert to list for modification
                    entry = list(entry)
                    # CRITICAL FIX: Convert pos from 0-based to 1-based for VCF/HGVS compatibility
                    # entry[1] is the 'pos' column (0-indexed in array, but column 2 in output)
                    entry[1] = int(entry[1]) + 1
                    # Clean and convert data efficiently
                    clean_entry = [str(val).replace("b'", "").replace("'", "") for val in entry]
                    f.write('\t'.join(clean_entry) + '\n')
        
        strand_info = f" ({strand} strand)" 
        print(f"✅ Processed {gene}:{hotspot}{strand_info} → {output_file}")
        return True
        
    except Exception as e:
        print(f"❌ Error processing {gene}:{hotspot}: {str(e)}")
        return False


def run_pysamstats_batch(bam_file, fasta_file, sample_id, output_dir, threads=4):
    """
    Run pysamstats for all mutations in batch with parallel processing.

    Args:
        bam_file: Path to BAM file
        fasta_file: Path to reference FASTA
        sample_id: Sample identifier
        output_dir: Output directory for results
        threads: Number of threads for parallel processing
    """

    # Load hotspots configuration
    print(f"📋 Loading hotspots configuration from {HOTSPOTS_CONFIG}")
    mutations_data = load_hotspots_config()

    print(f"📊 Processing {len(mutations_data)} hotspot mutations for {sample_id}")
    print(f"Using {threads} threads for parallel processing")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Process mutations in parallel
    process_func = partial(
        process_single_mutation,
        bam_file=bam_file,
        fasta_file=fasta_file,
        sample_id=sample_id,
        output_dir=output_dir
    )

    success_count = 0
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(process_func, mutations_data))
        success_count = sum(results)

    print(f"✅ Successfully processed {success_count}/{len(mutations_data)} mutations")
    
    # Create special IKZF1 output for compatibility (as expected by the rule)
    ikzf1_source = os.path.join(output_dir, f"{sample_id}_IKZF1_N159Y.tsv")
    ikzf1_target = os.path.join(output_dir, f"{sample_id}_IKZF1.csv")
    
    if os.path.exists(ikzf1_source):
        # Convert TSV to CSV for IKZF1 compatibility
        df = pd.read_csv(ikzf1_source, sep='\t')
        df.to_csv(ikzf1_target, index=False)
        print(f"✅ Created IKZF1 compatibility file: {ikzf1_target}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python run_pysamstats_optimized.py <bam_file> <fasta_file> <sample_id> <output_dir>")
        sys.exit(1)
    
    bam_file, fasta_file, sample_id, output_dir = sys.argv[1:5]
    
    # Get threads from environment or default to 4
    threads = int(os.environ.get('THREADS', 4))
    
    run_pysamstats_batch(bam_file, fasta_file, sample_id, output_dir, threads)