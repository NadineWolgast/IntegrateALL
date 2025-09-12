#!/usr/bin/env python3
"""
Optimized script to calculate total mapped reads from BAM file.
Replaces shell-based samtools command with Python pysam for better performance.
"""

import pysam
import sys
import os

def calculate_total_mapped_reads(bam_file, output_file):
    """
    Calculate total mapped reads using pysam (faster than shell samtools).
    Args:
        bam_file: Path to BAM file
        output_file: Path to output text file
    """
    try:
        # Open BAM file with pysam for faster processing
        with pysam.AlignmentFile(bam_file, "rb") as bamfile:
            # Count mapped reads (exclude flag 4 = unmapped)
            mapped_count = sum(1 for read in bamfile if not read.is_unmapped)
        
        # Write result to output file
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            f.write(str(mapped_count) + '\n')
            
        print(f"✅ Calculated {mapped_count:,} mapped reads from {bam_file}")
        
    except Exception as e:
        print(f"❌ Error processing {bam_file}: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_total_mapped_reads.py <bam_file> <output_file>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(bam_file):
        print(f"❌ BAM file not found: {bam_file}", file=sys.stderr)
        sys.exit(1)
    
    calculate_total_mapped_reads(bam_file, output_file)