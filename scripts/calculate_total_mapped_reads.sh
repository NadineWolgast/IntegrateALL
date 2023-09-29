#!/bin/bash
set -euo pipefail

# Replace 'sample_id' with the actual sample ID
sample_id="{wildcards.sample_id}"

# Calculate the total mapped reads using SAMtools
total_mapped_reads=$(samtools view -F 4 -c "STAR_output/${sample_id}/Aligned.sortedByCoord.out.bam")

# Write the total mapped reads to the output file
echo "$total_mapped_reads" > "data/total_mapped_reads/${sample_id}.txt"
