#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create an empty FusionCatcher output file when classification was successful with Arriba only.
"""

import sys
import os

def create_empty_fusioncatcher_file(output_file):
    """
    Create an empty FusionCatcher file with proper header.
    
    Args:
        output_file: Path to create the empty FusionCatcher file
    """
    
    # Ensure output directory exists (Python 2.7 compatible)
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # FusionCatcher header
    header = "Gene_1_symbol(5end_fusion_partner)\tGene_2_symbol(3end_fusion_partner)\tFusion_description\tCounts_of_common_mapping_reads\tSpanning_pairs\tSpanning_unique_reads\tLongest_anchor_found\tFusion_finding_method\tFusion_point_for_gene_1(5end_fusion_partner)\tFusion_point_for_gene_2(3end_fusion_partner)\tGene_1_id(5end_fusion_partner)\tGene_2_id(3end_fusion_partner)\tExon_1_id(5end_fusion_partner)\tExon_2_id(3end_fusion_partner)\tFusion_sequence\tPredicted_effect"
    
    # Write header to file
    with open(output_file, 'w') as f:
        f.write(header + '\n')
    
    print("Created empty FusionCatcher file: {}".format(output_file))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_empty_fusioncatcher.py <output_file>")
        sys.exit(1)
        
    output_file = sys.argv[1]
    create_empty_fusioncatcher_file(output_file)