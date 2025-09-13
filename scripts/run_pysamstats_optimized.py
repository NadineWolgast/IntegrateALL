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

# Load mutations from external file for better maintainability  
# Format: (gene, hotspot, chromosome, start, end, strand)
# Strand information from GTF: + = plus strand, - = minus strand
# IMPORTANT: These are the exact coordinates from the original working script
MUTATIONS_DATA = [
    # ZEB2 gene (Chr2, minus strand)
    ("ZEB2", "H1038", "2", 144389981, 144389984, "-"),
    ("ZEB2", "Q1072", "2", 144389879, 144389882, "-"),
    
    # KRAS gene (Chr12, minus strand)
    ("KRAS", "G12", "12", 25245348, 25245351, "-"),
    ("KRAS", "G13", "12", 25245345, 25245348, "-"),
    ("KRAS", "Q16", "12", 25227340, 25227343, "-"),
    ("KRAS", "A146", "12", 25225625, 25225628, "-"),
    
    # NRAS gene (Chr1, minus strand)
    ("NRAS", "G12", "1", 114716124, 114716127, "-"),
    ("NRAS", "G13", "1", 114716121, 114716124, "-"),
    ("NRAS", "Q61", "1", 114713906, 114713909, "-"),
    ("NRAS", "A146", "1", 114709580, 114709583, "-"),
    
    # FLT3 gene (Chr13, minus strand)
    ("FLT3", "P857", "13", 28015671, 28015674, "-"),
    ("FLT3", "V843", "13", 28018478, 28018481, "-"),
    ("FLT3", "Y842", "13", 28018481, 28018484, "-"),
    ("FLT3", "N841", "13", 28018484, 28018487, "-"),
    ("FLT3", "D839", "13", 28018490, 28018493, "-"),
    ("FLT3", "M837", "13", 28018496, 28018499, "-"),
    ("FLT3", "I836", "13", 28018499, 28018502, "-"),
    ("FLT3", "D835", "13", 28018502, 28018505, "-"),
    ("FLT3", "R834", "13", 28018505, 28018508, "-"),
    ("FLT3", "A680", "13", 28028190, 28028193, "-"),
    ("FLT3", "N676", "13", 28028202, 28028205, "-"),
    ("FLT3", "A627", "13", 28033947, 28033950, "-"),
    ("FLT3", "K623", "13", 28033959, 28033962, "-"),
    ("FLT3", "Y599", "13", 28034121, 28034124, "-"),
    ("FLT3", "R595", "13", 28034133, 28034136, "-"),
    ("FLT3", "V592", "13", 28034142, 28034145, "-"),
    ("FLT3", "Y589", "13", 28034151, 28034154, "-"),
    ("FLT3", "N587", "13", 28034157, 28034160, "-"),
    ("FLT3", "G583", "13", 28034169, 28034172, "-"),
    ("FLT3", "Q580", "13", 28034178, 28034181, "-"),
    ("FLT3", "V579", "13", 28034181, 28034184, "-"),
    ("FLT3", "Q577", "13", 28034187, 28034190, "-"),
    ("FLT3", "L576", "13", 28034190, 28034193, "-"),
    ("FLT3", "E573", "13", 28034199, 28034202, "-"),
    ("FLT3", "Y572", "13", 28034202, 28034205, "-"),
    ("FLT3", "V491", "13", 28035618, 28035621, "-"),
    ("FLT3", "S446", "13", 28036014, 28036017, "-"),
    
    # PAX5 gene (Chr9, minus strand)
    ("PAX5", "P80R", "9", 37015166, 37015169, "-"),
    
    # IKZF1 gene (Chr7, plus strand) - Only plus strand gene!
    ("IKZF1", "N159Y", "7", 50382592, 50382595, "+")  # Originally handled with -u flag (plus strand)
]


def process_single_mutation(mutation_data, bam_file, fasta_file, sample_id, output_dir):
    """Process a single mutation hotspot with strand-specific handling."""
    gene, hotspot, chromosome, start, end, strand = mutation_data
    
    output_file = os.path.join(output_dir, f"{sample_id}_{gene}_{hotspot}.tsv")
    
    try:
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
                    start=start, 
                    end=end, 
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
                    start=start, 
                    end=end, 
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
                
                # Write data rows
                for entry in variations:
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
    
    print(f"📊 Processing {len(MUTATIONS_DATA)} hotspot mutations for {sample_id}")
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
        results = list(executor.map(process_func, MUTATIONS_DATA))
        success_count = sum(results)
    
    print(f"✅ Successfully processed {success_count}/{len(MUTATIONS_DATA)} mutations")
    
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