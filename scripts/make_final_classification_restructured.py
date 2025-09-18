#!/usr/bin/env python3
"""
Restructured IntegrateALL final classification script.
3-Phase approach: Data Collection -> Simple Matching -> Complex Rules
"""

import csv
import os
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Driver fusion gene list - optimized as set for O(1) lookup
DRIVER_FUSION_GENES = {
    "ABL1", "ABL2", "ACIN1", "ADAMTSL5", "AFF1", "AMPH", "ANTXR1", "ARID1B", "ATF7IP",
    "AUTS2", "BCL2", "BCL6", "BCL9", "BCR", "BMP2K", "BRD9", "C7ORF72",
    "CASC15", "CBFA2T2", "CBFA2T3", "CBL", "CD163", "CEBPA", "CEBPB", "CEBPD", "CEBPE",
    "CENPC", "CLTC", "CREBBP", "CRLF2", "CSF1R", "CUX1", "DACH1", "DACH2", "DAZAP1",
    "DBX1", "DCPS", "DDX42", "DGKH", "DMRTA2", "DUX4", "EBF1", "ELMO1", "ELN", "EP300",
    "EPOR", "EPS15", "ERC1", "ESRRA", "ETV6", "EWSR1", "EXOSC2", "EXTL1", "FAM136A",
    "FBRSL1", "FIP1L1", "FKBP15", "FOXJ2", "FOXO3", "FOXP2", "GATA2DA", "GATAD2A",
    "GOPC", "HLF", "HNRNPM", "HNRNPUL1", "HOXA9", "ID4", "IGH@", "IGK", "IKZF1",
    "IL2RB", "IL7-R", "JAK2", "KANK1", "KMT2A", "LAIR1", "LEF1", "LSM14A", "LYN",
    "MBNL1", "MED12", "MEF2D", "MEIS2", "MLLT1", "MLLT10", "MLLT3", "MPRIP", "MYB",
    "MYC", "MYH9", "MYO18B", "NCOA5", "NCOR1", "NIPBL", "NOL4L", "NTRK3", "NUP153",
    "NUP214", "NUTM1", "OFD1", "P2RY8", "PAG1", "PAX5", "PBX1", "PCGF5", "PCM1",
    "PDGFRA", "PDGFRB", "PML", "PPFIBP1", "PTK2B", "PYGO2", "QSOX1", "RANBP2",
    "RCSD1", "RFX3", "RHOXF2B", "RNFT2", "ROS1", "RUNX1", "SFPQ", "SLC12A6", "SLC30A7",
    "SMARCA2", "SNX1", "SNX2", "SNX29", "SRRM1", "SS18", "SSBP2", "STIM2", "STRN3",
    "TAF15", "TAF3", "TBL1XR1", "TCF3", "TCF4", "TERF2", "THADA", "TMEM2", "TMPRSS9",
    "TMTC1", "TNIP1", "TNS3", "TPR", "TYK2", "UBASH3B", "UBTF", "USP2", "USP25",
    "WDR37", "WDR5", "ZC3HAV1", "ZEB2", "ZFAND3", "ZMIZ1", "ZMYND8", "ZNF274",
    "ZNF276", "ZNF340", "ZNF362", "ZNF384", "ZNF521", "ZNF618", "ZPBP"
}

# Subtype-specific classification rules
SUBTYPE_RULES = {
    'PAX5 P80R': {
        'evidence_type': 'snv',
        'required_columns': ['ALLCatchR', 'PAX5_P80R'],
        'required_values': {'PAX5_P80R': True},
        'confidence_policy': 'any',  # any confidence value accepted
        'fusion_policy': 'optional',  # fusions are optional
        'secondary_driver_policy': 'ignore'  # ignore secondary drivers
    },
    'IKZF1 N159Y': {
        'evidence_type': 'snv', 
        'required_columns': ['ALLCatchR', 'IKZF1_N159Y'],
        'required_values': {'IKZF1_N159Y': True},
        'confidence_policy': 'any',
        'fusion_policy': 'optional',
        'secondary_driver_policy': 'ignore'
    },
    'with mutated ZEB2 (p.H1038R)/IGH::CEBPE (provisional entity)': {
        'evidence_type': 'snv',
        'required_columns': ['ALLCatchR', 'Confidence', 'ZEB2_H1038R'],
        'required_values': {'ZEB2_H1038R': True},
        'confidence_policy': 'restricted',
        'allowed_confidence': ['high-confidence', 'candidate'],
        'fusion_policy': 'optional',
        'secondary_driver_policy': 'ignore'
    },
    'Ph-like': {
        'evidence_type': 'fusion',
        'required_columns': ['ALLCatchR', 'fusion'],
        'confidence_policy': 'flexible',  # ignore if empty in Class_test.csv
        'fusion_policy': 'required',
        'min_fusion_reads': 1,
        'secondary_driver_policy': 'manual_curation'  # secondary drivers → manual curation
    },
    'BCL2/MYC': {
        'evidence_type': 'fusion',
        'required_columns': ['ALLCatchR', 'fusion'],
        'confidence_policy': 'flexible',
        'fusion_policy': 'required', 
        'expected_fusions': ['BCL2', 'MYC', 'BCL6', 'IGH@'],
        'min_fusion_reads': 1,
        'secondary_driver_policy': 'manual_curation'
    },
    'Hyperdiploid': {
        'evidence_type': 'karyotype',
        'required_columns': ['ALLCatchR', 'karyotype_classifier'], 
        'required_values': {'karyotype_classifier': 'Hyperdiploid'},
        'confidence_policy': 'any',
        'fusion_policy': 'conditional',  # low reads ignored, high reads → manual
        'fusion_read_threshold': 3,  # >=3 reads → manual curation
        'secondary_driver_policy': 'read_threshold'
    },
    'DUX4': {
        'evidence_type': 'fusion',
        'required_columns': ['ALLCatchR', 'fusion'],
        'confidence_policy': 'flexible',
        'fusion_policy': 'required',
        'expected_fusions': ['DUX4'],
        'min_fusion_reads': 3,  # DUX4 fusions need >2 reads for non-DUX4 subtypes
        'secondary_driver_policy': 'ignore'
    }
    # Add more subtypes as needed...
}


class ClassificationProcessor:
    """Restructured processor with 3-phase approach."""
    
    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.data = {
            'allcatchr': {},
            'karyotype': {},
            'snvs': {},
            'fusions': [],
            'summary': {}
        }
        
    # ========== PHASE 1: DATA COLLECTION ==========
    
    def collect_allcatchr_data(self, allcatchr_file):
        """Collect ALLCatchR predictions."""
        try:
            df = pd.read_csv(allcatchr_file, sep='\t')
            if df.empty:
                self.data['allcatchr'] = {'subtype': '', 'confidence': '', 'ph_pos': ''}
                return
                
            row = df.iloc[0]
            self.data['allcatchr'] = {
                'subtype': row["Prediction"],
                'confidence': row["Confidence"],
                'ph_pos': row.get("BCR_ABL1_maincluster_pred", "not Ph-pos predicted")
            }
            logger.info(f"✅ ALLCatchR: {self.data['allcatchr']['subtype']} ({self.data['allcatchr']['confidence']})")
            
        except Exception as e:
            logger.warning(f"Error loading ALLCatchR data: {e}")
            self.data['allcatchr'] = {'subtype': '', 'confidence': '', 'ph_pos': ''}
    
    def collect_karyotype_data(self, karyotype_file):
        """Collect karyotype classification."""
        try:
            df = pd.read_csv(karyotype_file)
            if df.empty:
                self.data['karyotype'] = {'prediction': '', 'score': '', 'raw_prediction': ''}
                return
                
            raw_prediction = df.iloc[0]["Prediction"]
            raw_score = df.iloc[0]["Score"]
            
            # Store raw data first - filtering happens after ALLCatchR data is collected
            self.data['karyotype'] = {
                'prediction': raw_prediction,
                'score': raw_score,
                'raw_prediction': raw_prediction
            }
            
            logger.info(f"✅ Karyotype (raw): {raw_prediction}")
            
        except Exception as e:
            logger.warning(f"Error loading karyotype data: {e}")
            self.data['karyotype'] = {'prediction': '', 'score': '', 'raw_prediction': ''}
    
    def apply_karyotype_filtering(self):
        """Apply iAMP21 filtering after ALLCatchR data is collected."""
        raw_prediction = self.data['karyotype']['raw_prediction']
        allcatchr_subtype = self.data['allcatchr']['subtype']
        
        # Apply iAMP21 filtering: only use iAMP21 karyotype if ALLCatchR subtype is also iAMP21
        if raw_prediction == 'iAMP21' and allcatchr_subtype != 'iAMP21':
            logger.info(f"🔧 iAMP21 karyotype filtered out - ALLCatchR subtype is '{allcatchr_subtype}', not 'iAMP21'")
            # Use fallback - could be 'other' or check for second highest prediction
            # For now, use 'other' as fallback
            filtered_prediction = 'other'
            self.data['karyotype']['prediction'] = filtered_prediction
            logger.info(f"✅ Karyotype: {filtered_prediction} (filtered from {raw_prediction})")
        else:
            logger.info(f"✅ Karyotype: {raw_prediction} (no filtering needed)")
    
    def collect_snv_data(self, hotspot_dir):
        """Collect hotspot SNV data - both standard and extended."""
        # Standard hotspots for matching
        snv_patterns = {
            "PAX5_P80R": False,
            "IKZF1_N159Y": False, 
            "ZEB2_H1038R": False
        }
        
        # Extended hotspots for reporting
        extended_hotspots = []
        
        try:
            hotspot_path = Path(hotspot_dir)
            if not hotspot_path.exists():
                self.data['snvs'] = snv_patterns
                self.data['extended_hotspots'] = extended_hotspots
                return
                
            for file_path in hotspot_path.iterdir():
                filename = file_path.name
                
                # Standard hotspots for matching
                if filename.startswith("PAX5_P80R"):
                    snv_patterns["PAX5_P80R"] = True
                elif filename.startswith("IKZF1_N159Y"):
                    snv_patterns["IKZF1_N159Y"] = True
                elif filename.startswith("ZEB2_H1038"):
                    snv_patterns["ZEB2_H1038R"] = True
                
                # All hotspots for extended reporting
                if "_" in filename and filename.endswith(('.txt', '.csv', '.tsv')):
                    # Extract mutation info: GENE_MUTATION_xxx.txt -> GENE_MUTATION
                    mutation_parts = filename.split('_')
                    if len(mutation_parts) >= 2:
                        gene = mutation_parts[0]
                        mutation = '_'.join(mutation_parts[1:]).split('.')[0]  # Remove extension
                        extended_hotspots.append(f"{gene}_{mutation}")
            
            self.data['snvs'] = snv_patterns
            self.data['extended_hotspots'] = extended_hotspots
            
            standard_snvs = [k for k, v in snv_patterns.items() if v]
            logger.info(f"✅ Standard SNVs: {standard_snvs if standard_snvs else 'None'}")
            if extended_hotspots:
                logger.info(f"✅ Extended hotspots: {extended_hotspots}")
            
        except Exception as e:
            logger.warning(f"Error checking hotspot files: {e}")
            self.data['snvs'] = snv_patterns
            self.data['extended_hotspots'] = extended_hotspots
    
    def collect_fusion_data(self, fusioncatcher_file, arriba_file, classification_df):
        """Collect driver fusions from both callers, filtered against Class_test.csv."""
        fusions = []
        
        # Load FusionCatcher data
        try:
            fc_df = pd.read_csv(fusioncatcher_file, sep='\t', skiprows=1, header=None)
            if not fc_df.empty and len(fc_df.columns) >= 6:
                logger.info(f"🔍 FusionCatcher: Processing {len(fc_df)} fusion candidates...")
                for _, row in fc_df.iterrows():
                    gene1, gene2, spanning_reads = row[0], row[1], row[5]
                    logger.info(f"   Checking fusion: {gene1}::{gene2}")
                    
                    if self._is_driver_fusion_in_classtest(gene1, gene2, classification_df):
                        logger.info(f"   ✅ Driver fusion accepted: {gene1}::{gene2}")
                        fusion_data = {
                            'gene_1': gene1,
                            'gene_2': gene2,
                            'caller': 'FusionCatcher',
                            'spanning_reads': spanning_reads
                        }
                        
                        # Add breakpoints and sequence if available
                        if len(row) >= 8:
                            fusion_data.update({
                                'breakpoint_1': row[2] if pd.notna(row[2]) else '',
                                'breakpoint_2': row[3] if pd.notna(row[3]) else '', 
                                'sequence_data': row[6] if pd.notna(row[6]) else ''
                            })
                        else:
                            fusion_data.update({
                                'breakpoint_1': '', 'breakpoint_2': '', 'sequence_data': ''
                            })
                        
                        fusions.append(fusion_data)
                    else:
                        logger.info(f"   ❌ Fusion rejected: {gene1}::{gene2}")
        except Exception as e:
            logger.warning(f"Error loading FusionCatcher data: {e}")
        
        # Load Arriba data  
        try:
            arriba_df = pd.read_csv(arriba_file, sep='\t', skiprows=1, header=None)
            if not arriba_df.empty and len(arriba_df.columns) >= 2:
                for _, row in arriba_df.iterrows():
                    gene1, gene2 = row[0], row[1]
                    
                    # Spanning reads: try column 11 first (real Arriba), then column 5 (test data)
                    if len(row) > 11:
                        spanning_reads = row[11]
                    elif len(row) > 5:
                        spanning_reads = row[5]  # For test data
                    else:
                        spanning_reads = 1  # Default fallback
                    
                    if self._is_driver_fusion_in_classtest(gene1, gene2, classification_df):
                        fusion_data = {
                            'gene_1': gene1,
                            'gene_2': gene2,
                            'caller': 'Arriba',
                            'spanning_reads': spanning_reads
                        }
                        
                        # Add breakpoints and sequence if available
                        if len(row) >= 8:
                            fusion_data.update({
                                'breakpoint_1': row[2] if pd.notna(row[2]) else '',
                                'breakpoint_2': row[3] if pd.notna(row[3]) else '',
                                'sequence_data': row[6] if len(row) > 6 and pd.notna(row[6]) else ''
                            })
                        else:
                            fusion_data.update({
                                'breakpoint_1': '', 'breakpoint_2': '', 'sequence_data': ''
                            })
                        
                        fusions.append(fusion_data)
        except Exception as e:
            logger.warning(f"Error loading Arriba data: {e}")
        
        self.data['fusions'] = fusions
        logger.info(f"✅ Driver fusions found: {len(fusions)}")
        for fusion in fusions:
            logger.info(f"   {fusion['gene_1']}::{fusion['gene_2']} ({fusion['caller']}, reads: {fusion['spanning_reads']})")
    
    def _is_driver_fusion_in_classtest(self, gene1, gene2, classification_df):
        """Check if fusion is driver gene pair and exists in Class_test.csv."""
        logger.info(f"      🔍 Checking genes: {gene1}, {gene2}")
        
        # Normalize IGH genes for matching - remove @ for Class_test.csv matching
        gene1_norm = 'IGH' if gene1.startswith('IGH') else gene1
        gene2_norm = 'IGH' if gene2.startswith('IGH') else gene2
        logger.info(f"      🔍 Normalized for matching: {gene1_norm}, {gene2_norm}")
        
        # Check if both genes are driver genes (use original names with @ for driver gene check)
        gene1_driver_check = 'IGH@' if gene1.startswith('IGH') else gene1
        gene2_driver_check = 'IGH@' if gene2.startswith('IGH') else gene2
        logger.info(f"      🔍 Driver gene check: {gene1_driver_check} in DRIVER_GENES={gene1_driver_check in DRIVER_FUSION_GENES}, {gene2_driver_check} in DRIVER_GENES={gene2_driver_check in DRIVER_FUSION_GENES}")
        
        if not (gene1_driver_check in DRIVER_FUSION_GENES and gene2_driver_check in DRIVER_FUSION_GENES):
            logger.info(f"      ❌ Not both driver genes")
            return False
        
        # Check if combination exists in Class_test.csv (both orientations)
        # Use original gene names for Class_test.csv lookup (Class_test.csv has IGH@ not IGH)
        exists1 = ((classification_df['Gene_1_symbol(5end_fusion_partner)'] == gene1) &
                  (classification_df['Gene_2_symbol(3end_fusion_partner)'] == gene2)).any()
        exists2 = ((classification_df['Gene_1_symbol(5end_fusion_partner)'] == gene2) &
                  (classification_df['Gene_2_symbol(3end_fusion_partner)'] == gene1)).any()
        
        logger.info(f"      🔍 Class_test.csv check: {gene1}::{gene2}={exists1}, {gene2}::{gene1}={exists2}")
        
        exists = exists1 or exists2
        logger.info(f"      {'✅' if exists else '❌'} Final result: {exists}")
        return exists
    
    def create_summary_dataframe(self):
        """Create summary DataFrame for matching."""
        # Handle no fusions case
        if not self.data['fusions']:
            self.data['summary'] = pd.DataFrame([{
                'ALLCatchR': self.data['allcatchr']['subtype'],
                'Ph-pos': self.data['allcatchr']['ph_pos'],
                'Confidence': self.data['allcatchr']['confidence'],
                'Gene_1_symbol(5end_fusion_partner)': None,
                'Gene_2_symbol(3end_fusion_partner)': None,
                'Fusioncaller': None,
                'Unique_spanning_reads': None,
                'karyotype_classifier': self.data['karyotype']['prediction'],
                'PAX5_P80R': self.data['snvs']['PAX5_P80R'],
                'IKZF1_N159Y': self.data['snvs']['IKZF1_N159Y'],
                'ZEB2_H1038R': self.data['snvs']['ZEB2_H1038R']
            }])
        else:
            # Create one row per fusion
            rows = []
            for fusion in self.data['fusions']:
                # Use original gene names for matching (Class_test.csv has IGH@ not IGH)
                rows.append({
                    'ALLCatchR': self.data['allcatchr']['subtype'],
                    'Ph-pos': self.data['allcatchr']['ph_pos'],
                    'Confidence': self.data['allcatchr']['confidence'],
                    'Gene_1_symbol(5end_fusion_partner)': fusion['gene_1'],
                    'Gene_2_symbol(3end_fusion_partner)': fusion['gene_2'],
                    'Fusioncaller': fusion['caller'],
                    'Unique_spanning_reads': fusion['spanning_reads'],
                    'karyotype_classifier': self.data['karyotype']['prediction'],
                    'PAX5_P80R': self.data['snvs']['PAX5_P80R'],
                    'IKZF1_N159Y': self.data['snvs']['IKZF1_N159Y'],
                    'ZEB2_H1038R': self.data['snvs']['ZEB2_H1038R']
                })
            
            self.data['summary'] = pd.DataFrame(rows)
        
        logger.info(f"✅ Summary created: {len(self.data['summary'])} row(s)")
    
    def _flexible_confidence_and_fusion_merge(self, summary_df, classification_df, match_cols):
        """Merge with flexible Confidence and fusion orientation handling."""
        logger.info(f"🔍 DEBUG: Attempting merge with columns: {match_cols}")
        logger.info(f"🔍 DEBUG: Summary data shape: {summary_df.shape}")
        logger.info(f"🔍 DEBUG: Classification data shape: {classification_df.shape}")
        
        # Show sample data being matched
        if not summary_df.empty:
            logger.info("🔍 DEBUG: Summary data sample:")
            for col in match_cols:
                if col in summary_df.columns:
                    logger.info(f"      {col}: {summary_df[col].iloc[0]}")
        
        # Try both orientations
        matches = self._try_confidence_aware_merge(summary_df, classification_df, match_cols, swapped=False)
        if not matches.empty:
            return matches
            
        # If no matches and we have fusion columns, try swapped orientation
        if ('Gene_1_symbol(5end_fusion_partner)' in match_cols and 
            'Gene_2_symbol(3end_fusion_partner)' in match_cols):
            
            logger.info("🔄 Trying swapped fusion orientation...")
            
            # Create swapped version of summary for matching
            summary_swapped = summary_df.copy()
            summary_swapped['Gene_1_symbol(5end_fusion_partner)'] = summary_df['Gene_2_symbol(3end_fusion_partner)']
            summary_swapped['Gene_2_symbol(3end_fusion_partner)'] = summary_df['Gene_1_symbol(5end_fusion_partner)']
            
            matches = self._try_confidence_aware_merge(summary_swapped, classification_df, match_cols, swapped=True)
            if not matches.empty:
                return matches
        
        logger.info("🔍 DEBUG: No matches found in either orientation")
        return pd.DataFrame()  # Return empty if no matches found
    
    def _try_confidence_aware_merge(self, summary_df, classification_df, match_cols, swapped=False):
        """Try merge with intelligent Confidence handling."""
        orientation = "swapped" if swapped else "exact"
        logger.info(f"🔍 DEBUG: Trying {orientation} orientation merge...")
        
        if swapped:
            logger.info("🔍 DEBUG: Swapped summary data:")
            for col in match_cols:
                if col in summary_df.columns:
                    logger.info(f"      {col}: {summary_df[col].iloc[0]}")
        
        # If Confidence is not in match_cols, use standard merge
        if 'Confidence' not in match_cols:
            matches = pd.merge(
                summary_df, classification_df,
                on=match_cols, how='inner',
                suffixes=('', '_class')
            )
            logger.info(f"🔍 DEBUG: {orientation} orientation matches (no Confidence): {len(matches)}")
            return matches
        
        # Confidence-aware matching: split classification data
        non_confidence_cols = [col for col in match_cols if col != 'Confidence']
        
        # First, get potential matches on all non-Confidence columns
        potential_matches = pd.merge(
            summary_df, classification_df,
            on=non_confidence_cols, how='inner',
            suffixes=('', '_class')
        )
        
        if potential_matches.empty:
            logger.info(f"🔍 DEBUG: {orientation} orientation - no potential matches on non-Confidence columns")
            return potential_matches
        
        logger.info(f"🔍 DEBUG: {orientation} orientation - found {len(potential_matches)} potential matches on non-Confidence columns")
        
        # Now filter based on Confidence logic
        sample_confidence = summary_df['Confidence'].iloc[0]
        logger.info(f"🔍 DEBUG: Sample Confidence: {sample_confidence}")
        
        final_matches = []
        for _, match_row in potential_matches.iterrows():
            class_confidence = match_row['Confidence_class']
            
            # Check if classification confidence is empty/null
            is_empty_confidence = (pd.isnull(class_confidence) or 
                                 class_confidence == '' or 
                                 str(class_confidence).strip() == '')
            
            if is_empty_confidence:
                # Empty confidence in classification - always matches
                logger.info(f"🔍 DEBUG: Match accepted - empty Confidence in class_test")
                final_matches.append(match_row)
            elif str(class_confidence).lower() == str(sample_confidence).lower():
                # Specific confidence values must match
                logger.info(f"🔍 DEBUG: Match accepted - Confidence values match: {sample_confidence}")
                final_matches.append(match_row)
            else:
                logger.info(f"🔍 DEBUG: Match rejected - Confidence mismatch: {sample_confidence} vs {class_confidence}")
        
        result = pd.DataFrame(final_matches) if final_matches else pd.DataFrame()
        logger.info(f"🔍 DEBUG: {orientation} orientation final matches: {len(result)}")
        return result
    
    # ========== PHASE 2: SUBTYPE-AWARE MATCHING ==========
    
    def find_subtype_aware_matches(self, classification_df):
        """Apply subtype-specific matching rules."""
        subtype = self.data['allcatchr']['subtype']
        logger.info(f"🎯 SUBTYPE-AWARE MATCHING: {subtype}")
        
        # Check if we have specific rules for this subtype
        if subtype not in SUBTYPE_RULES:
            logger.info(f"⚠️ No specific rules for {subtype} - using fallback to general matching")
            return self.find_exact_matches(classification_df)
        
        rules = SUBTYPE_RULES[subtype]
        logger.info(f"📋 Rules for {subtype}: {rules['evidence_type']} evidence type")
        
        # Apply subtype-specific logic
        if rules['evidence_type'] == 'snv':
            return self._match_snv_subtype(classification_df, rules)
        elif rules['evidence_type'] == 'fusion':
            return self._match_fusion_subtype(classification_df, rules)
        elif rules['evidence_type'] == 'karyotype':
            return self._match_karyotype_subtype(classification_df, rules)
        else:
            logger.info(f"⚠️ Unknown evidence type: {rules['evidence_type']}")
            return self.find_exact_matches(classification_df)
    
    def _match_snv_subtype(self, classification_df, rules):
        """Handle SNV-based subtypes like PAX5 P80R, IKZF1 N159Y."""
        logger.info("🧬 Matching SNV-based subtype...")
        
        # Get summary data
        summary_df = self.data['summary']
        if summary_df.empty:
            return None, 'no_data'
        
        # Check if required SNV is present
        for snv_col, required_value in rules.get('required_values', {}).items():
            sample_value = self.data.get('snvs', {}).get(snv_col, False)
            if sample_value != required_value:
                logger.info(f"❌ Required SNV not met: {snv_col}={sample_value}, expected {required_value}")
                return None, 'snv_mismatch'
        
        # Build matching columns based on rules
        match_cols = ['ALLCatchR']
        
        # Add required SNV columns
        for snv_col in rules.get('required_values', {}):
            if snv_col in ['PAX5_P80R', 'IKZF1_N159Y', 'ZEB2_H1038R']:
                match_cols.append(snv_col)
        
        # Handle Confidence based on policy
        if rules['confidence_policy'] == 'restricted':
            # Only allow specific confidence values
            sample_confidence = self.data['allcatchr']['confidence']
            allowed = rules.get('allowed_confidence', [])
            if sample_confidence not in allowed:
                logger.info(f"❌ Confidence {sample_confidence} not in allowed values: {allowed}")
                return None, 'confidence_mismatch'
            match_cols.append('Confidence')
        elif rules['confidence_policy'] != 'any':
            match_cols.append('Confidence')
        
        # Add other standard columns
        match_cols.extend(['Ph-pos', 'karyotype_classifier'])
        
        # Filter classification data for entries without fusions (SNV-only)
        classification_for_match = classification_df[
            classification_df['Gene_1_symbol(5end_fusion_partner)'].isnull()
        ]
        
        logger.info(f"🔍 Matching SNV subtype on columns: {match_cols}")
        logger.info(f"🔍 Classification entries without fusions: {len(classification_for_match)}")
        
        # Use flexible confidence matching
        matches = self._flexible_confidence_and_fusion_merge(
            summary_df.head(1), classification_for_match, match_cols
        )
        
        if not matches.empty:
            logger.info(f"✅ SNV subtype matches found: {len(matches)}")
            if len(matches) == 1:
                return matches, 'exact_single'
            else:
                return matches, 'exact_multiple'
        else:
            logger.info("❌ No SNV subtype matches found")
            return None, 'no_exact'
    
    def _match_fusion_subtype(self, classification_df, rules):
        """Handle fusion-based subtypes like Ph-like, BCL2/MYC."""
        logger.info("🔗 Matching fusion-based subtype...")
        
        # Check if we have required fusions
        if not self.data['fusions']:
            logger.info("❌ No fusions found for fusion-based subtype")
            return None, 'no_fusions'
        
        # Check for secondary drivers if policy requires it
        if rules.get('secondary_driver_policy') == 'manual_curation':
            secondary_drivers = self._detect_secondary_drivers(rules)
            if secondary_drivers:
                logger.info(f"⚠️ Secondary drivers detected: {secondary_drivers} → Manual curation")
                return None, 'secondary_drivers'
        
        # Use the existing fusion matching logic but with subtype-specific rules
        return self.find_exact_matches(classification_df)
    
    def _match_karyotype_subtype(self, classification_df, rules):
        """Handle karyotype-based subtypes like Hyperdiploid."""
        logger.info("🧬 Matching karyotype-based subtype...")
        
        # Check if required karyotype is present
        for karyo_col, required_value in rules.get('required_values', {}).items():
            sample_value = self.data['karyotype'].get('prediction', 'other')
            if sample_value != required_value:
                logger.info(f"❌ Required karyotype not met: {sample_value}, expected {required_value}")
                return None, 'karyotype_mismatch'
        
        logger.info(f"✅ Required karyotype confirmed: {sample_value}")
        
        # Handle fusion filtering based on policy
        if rules.get('fusion_policy') == 'conditional':
            # Filter out low-read fusions first
            self._filter_low_read_fusions(rules)
            
            # Check if any high-read fusions remain
            if self.data['fusions']:
                high_read_fusions = []
                for fusion in self.data['fusions']:
                    try:
                        reads = int(fusion['spanning_reads'])
                        high_read_fusions.append(f"{fusion['gene_1']}::{fusion['gene_2']} (reads: {reads})")
                    except:
                        continue
                
                if high_read_fusions:
                    logger.info(f"⚠️ High-read fusions detected after filtering: {high_read_fusions} → Manual curation")
                    return None, 'high_read_fusions'
        
        # For karyotype-based subtypes, match without fusion requirements
        # Get summary data  
        summary_df = self.data['summary']
        if summary_df.empty:
            return None, 'no_data'
        
        # Match on karyotype-specific columns (exclude fusion columns)
        match_cols = ['ALLCatchR', 'Ph-pos', 'karyotype_classifier']
        
        # Add confidence if policy requires it
        if rules.get('confidence_policy') != 'any':
            match_cols.append('Confidence')
        
        # Filter classification data for entries without fusions (karyotype-only)
        classification_for_match = classification_df[
            classification_df['Gene_1_symbol(5end_fusion_partner)'].isnull()
        ]
        
        logger.info(f"🔍 Matching karyotype subtype on columns: {match_cols}")
        logger.info(f"🔍 Classification entries without fusions: {len(classification_for_match)}")
        
        # Use flexible confidence matching (but no fusion columns)
        matches = self._flexible_confidence_and_fusion_merge(
            summary_df.head(1), classification_for_match, match_cols
        )
        
        if not matches.empty:
            logger.info(f"✅ Karyotype subtype matches found: {len(matches)}")
            if len(matches) == 1:
                return matches, 'exact_single'
            else:
                return matches, 'exact_multiple'
        else:
            logger.info("❌ No karyotype subtype matches found")
            return None, 'no_exact'
    
    def _detect_secondary_drivers(self, rules):
        """Detect secondary driver fusions that might cause conflicts."""
        secondary_drivers = []
        expected_genes = set(rules.get('expected_fusions', []))
        
        for fusion in self.data['fusions']:
            fusion_genes = {fusion['gene_1'], fusion['gene_2']}
            # If fusion contains genes not expected for this subtype
            if not fusion_genes.intersection(expected_genes):
                secondary_drivers.append(f"{fusion['gene_1']}::{fusion['gene_2']}")
        
        return secondary_drivers
    
    def _check_fusion_read_thresholds(self, rules):
        """Check if any fusions exceed read count thresholds."""
        threshold = rules.get('fusion_read_threshold', 3)
        high_read_fusions = []
        
        for fusion in self.data['fusions']:
            try:
                reads = int(fusion['spanning_reads'])
                if reads >= threshold:
                    high_read_fusions.append(f"{fusion['gene_1']}::{fusion['gene_2']} (reads: {reads})")
            except:
                continue
        
        return high_read_fusions
    
    def _filter_low_read_fusions(self, rules):
        """Remove low-read fusions based on thresholds."""
        threshold = rules.get('fusion_read_threshold', 3)
        original_count = len(self.data['fusions'])
        
        filtered_fusions = []
        for fusion in self.data['fusions']:
            try:
                reads = int(fusion['spanning_reads'])
                if reads < threshold:
                    logger.info(f"🔧 Filtering low-read fusion: {fusion['gene_1']}::{fusion['gene_2']} (reads: {reads} < {threshold})")
                else:
                    filtered_fusions.append(fusion)
            except:
                # Keep fusion if read count can't be parsed
                filtered_fusions.append(fusion)
        
        self.data['fusions'] = filtered_fusions
        logger.info(f"🔧 Fusions filtered: {original_count} → {len(filtered_fusions)}")
        
        # Recreate summary after filtering
        if filtered_fusions != self.data['fusions']:
            self.create_summary_dataframe()
    
    # ========== FALLBACK: GENERAL MATCHING ==========
    
    def find_exact_matches(self, classification_df):
        """Try simple 1:1 exact matching against Class_test.csv with smart ZEB2 handling."""
        # Rename columns for consistency
        classification_df = classification_df.rename(columns={
            'karyotype classifier': 'karyotype_classifier',
            'PAX5 P80R': 'PAX5_P80R',
            'IKZF1 N159Y': 'IKZF1_N159Y', 
            'ZEB2 H1038R': 'ZEB2_H1038R'
        })
        
        summary_df = self.data['summary']
        subtype = self.data['allcatchr']['subtype']
        
        # Define base matching columns (exclude ZEB2 for smart handling)
        base_match_cols = [
            'ALLCatchR', 'Ph-pos', 'Confidence', 'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)', 'karyotype_classifier',
            'PAX5_P80R', 'IKZF1_N159Y'
        ]
        
        # Smart ZEB2 matching based on subtype
        if subtype == 'CEBP':
            # CEBP subtype: ZEB2_H1038R must be True (exact match required)
            match_cols = base_match_cols + ['ZEB2_H1038R']
            logger.info("🎯 CEBP subtype: ZEB2_H1038R must be True (exact matching)")
        else:
            # All other subtypes: ZEB2_H1038R can be True or False (ignore in matching)
            match_cols = base_match_cols
            logger.info("🎯 Non-CEBP subtype: ZEB2_H1038R ignored in matching")
        
        # Remove None values for matching (no fusions case)
        summary_for_match = summary_df.copy()
        classification_for_match = classification_df.copy()
        
        # Handle NaN/None values in fusion columns
        if summary_for_match['Gene_1_symbol(5end_fusion_partner)'].isnull().all():
            # No fusions - match only on non-fusion columns
            match_cols = [col for col in match_cols if 'Gene_' not in col]
            summary_for_match = summary_for_match.dropna(subset=['ALLCatchR'])
            classification_for_match = classification_df[
                classification_df['Gene_1_symbol(5end_fusion_partner)'].isnull()
            ]
        
        # Confidence will be handled flexibly in the merge function
        
        # For non-CEBP subtypes with ZEB2=False in data, allow matching with ZEB2=True in class_test
        if subtype != 'CEBP' and not summary_for_match['ZEB2_H1038R'].iloc[0]:
            logger.info("🔧 ZEB2=False in data, allowing match with ZEB2=True in classification DB")
        
        # Perform exact merge on selected columns, but handle fusion orientation flexibility
        matches = self._flexible_confidence_and_fusion_merge(summary_for_match, classification_for_match, match_cols)
        
        if not matches.empty:
            logger.info(f"✅ Exact matches found: {len(matches)}")
            if len(matches) == 1:
                logger.info("✅ Single exact match - classification successful!")
                return matches, 'exact_single'
            else:
                logger.info("⚠️ Multiple exact matches - need complex rules")
                return matches, 'exact_multiple'
        else:
            logger.info("⚠️ No exact matches - need complex rules")
            return None, 'no_exact'
    
    # ========== PHASE 3: COMPLEX RULES ==========
    
    def apply_complex_rules(self, classification_df, matches=None, match_type=None):
        """Apply complex rules for multiple matches or inconsistencies."""
        logger.info("🔧 Applying complex rules...")
        
        # Apply DUX4 specific rules when needed
        if match_type == 'exact_multiple' and matches is not None:
            resolved_matches = self._resolve_multiple_matches(matches)
            if resolved_matches is not None:
                return resolved_matches, 'exact_single'
            else:
                return None, 'multiple_different_diagnoses'
        elif match_type == 'no_exact':
            return self._resolve_no_matches(classification_df)
        else:
            return matches, match_type
    
    def _resolve_multiple_matches(self, matches):
        """Resolve multiple exact matches - check if they have same diagnosis."""
        logger.info("🔧 Resolving multiple matches...")
        
        # Check if all matches have the same WHO-HAEM5 and ICC diagnosis
        unique_who = matches['WHO-HAEM5'].nunique()
        unique_icc = matches['ICC'].nunique()
        
        if unique_who == 1 and unique_icc == 1:
            # All matches have same diagnosis - accept (multiple fusions for same diagnosis allowed)
            logger.info(f"✅ Multiple matches with same diagnosis: {matches.iloc[0]['WHO-HAEM5']}")
            return matches.head(1)  # Take first match, they're all the same diagnosis
        else:
            # Different diagnoses - manual curation required
            logger.info(f"❌ Multiple matches with different diagnoses: WHO-HAEM5={unique_who}, ICC={unique_icc}")
            return None
    
    def _resolve_no_matches(self, classification_df):
        """Resolve case with no exact matches using flexible rules."""
        logger.info("🔧 Resolving no matches with flexible rules...")
        
        # Apply DUX4 rules here
        subtype = self.data['allcatchr']['subtype']
        
        # Filter DUX4 fusions for non-DUX4 subtypes
        if subtype != 'DUX4':
            valid_fusions = []
            for fusion in self.data['fusions']:
                is_dux4_fusion = (fusion['gene_1'] == 'DUX4' or fusion['gene_2'] == 'DUX4')
                if is_dux4_fusion:
                    # DUX4 fusion in non-DUX4 subtype - check spanning reads
                    try:
                        spanning_reads = int(fusion['spanning_reads'])
                        if spanning_reads > 2:
                            valid_fusions.append(fusion)
                            logger.info(f"✅ DUX4 fusion accepted: {fusion['gene_1']}::{fusion['gene_2']} (reads: {spanning_reads})")
                        else:
                            logger.info(f"❌ DUX4 fusion rejected: {fusion['gene_1']}::{fusion['gene_2']} (reads: {spanning_reads} <= 2)")
                    except:
                        logger.info(f"❌ DUX4 fusion rejected: {fusion['gene_1']}::{fusion['gene_2']} (invalid read count)")
                else:
                    valid_fusions.append(fusion)
            
            # Update fusions and recreate summary
            self.data['fusions'] = valid_fusions
            self.create_summary_dataframe()
            
            # Try exact matching again with filtered data
            exact_results = self.find_exact_matches(classification_df)
            if exact_results[0] is not None:
                return exact_results
        
        # No flexible matching - if no exact matches after DUX4 filtering, it's manual curation
        logger.info("❌ No exact matches found after applying rules")
        return None, 'no_match_after_rules'
    
    # ========== OUTPUT GENERATION ==========
    
    def generate_output_files(self, results, match_type, output_paths):
        """Generate all output files."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        if results is not None and not results.empty:
            # Successful classification (exact or flexible match)
            self._write_successful_classification(results, output_paths, match_type)
        else:
            self._write_manual_curation_needed(output_paths)
        
        # Generate driver fusion file
        self._generate_driver_fusion_file(output_driver)
    
    def _write_successful_classification(self, results_df, output_paths, match_type='exact'):
        """Write output for successful classification."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        # Save main results
        results_df.to_csv(output_csv, index=False)
        
        # Write text output
        entry = results_df.iloc[0]
        fusion_detail = self._format_fusion_detail(entry)
        
        with open(output_text, 'w') as f:
            f.write(
                f"Consistency between gene expression-based subtype-allocation "
                f"(ALLCatchR: {entry['ALLCatchR']} subtype, Confidence: {entry['Confidence']}) "
                f"and the genomic driver profile ({fusion_detail}, "
                f"RNA-Seq CNV karyotype classifier: {entry['karyotype_classifier']}, "
                f"subtype defining SNPs: "
                f"{'PAX5 P80R present' if entry['PAX5_P80R'] else 'PAX5 P80R absent'}, "
                f"{'IKZF1 N159Y present' if entry['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}) "
                f"supports a classification as\n\n"
                f"{entry['WHO-HAEM5']} according to WHO-HAEM5 (Alaggio R et al. Leukemia, 2022) and\n"
                f"{entry['ICC']} according to ICC (Arber D et al. Blood, 2022).\n\n"
            )
        
        # Write curation CSV
        self._write_curation_csv(entry, output_curation, "Automatic Classification")
    
    def _write_manual_curation_needed(self, output_paths):
        """Write output when manual curation is needed."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        # Create summary with default values
        summary_df = self.data['summary'].copy()
        summary_df['WHO-HAEM5'] = "IntegrateALL couldn't confirm the subtype in concordance with WHO or ICC classification"
        summary_df['ICC'] = ""
        
        # Save data
        summary_df.to_csv(output_csv, index=False)
        
        # Write curation needed text
        entry = summary_df.iloc[0]
        fusion_detail = self._format_fusion_detail(entry)
        
        with open(output_text, 'w') as f:
            f.write(
                f"IntegrateALL classification:\n\n"
                f"Gene expression-based subtype-allocation "
                f"(ALLCatchR: {entry['ALLCatchR']} subtype, Confidence: {entry['Confidence']}) "
                f"and the genomic driver profile ({fusion_detail}, "
                f"RNA-Seq CNV karyotype classifier: {entry['karyotype_classifier']}, "
                f"subtype defining SNPs: "
                f"{'PAX5 P80R present' if entry['PAX5_P80R'] else 'PAX5 P80R absent'}, "
                f"{'IKZF1 N159Y present' if entry['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}) "
                f"seem not to be consistent with an unambiguous diagnostic classification "
                f"according to WHO-HAEM5 (Alaggio R et al. Leukemia, 2022) / "
                f"ICC (Arber D et al. Blood, 2022)."
            )
        
        self._write_curation_csv(entry, output_curation, "Manual Curation")
    
    def _format_fusion_detail(self, entry):
        """Format fusion details for output text."""
        if pd.isna(entry.get('Fusioncaller')):
            return "No driver fusions identified"
        else:
            return f"{entry['Fusioncaller']}: {entry['Gene_1_symbol(5end_fusion_partner)']}::{entry['Gene_2_symbol(3end_fusion_partner)']}"
    
    def _write_curation_csv(self, entry, output_curation, classification_type):
        """Write curation CSV output matching expected format."""
        with open(output_curation, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Header matching expected format
            writer.writerow([
                "Sample_ID", "Classification", "Subtype", "Confidence", "Fusion_details",
                "Karyotype_classifier", "PAX5_P80R", "IKZF1_N159Y", "ZEB2_H1038R",
                "WHO-HAEM5", "ICC"
            ])
            
            # Format fusion details based on classification type
            if classification_type == "Automatic Classification":
                fusion_details = self._format_fusion_summary()
            else:
                fusion_details = self._format_all_fusions()
            
            writer.writerow([
                self.sample_id,
                classification_type,
                entry['ALLCatchR'],
                entry['Confidence'],
                fusion_details,
                entry['karyotype_classifier'],
                "PAX5 P80R present" if entry['PAX5_P80R'] else "PAX5 P80R absent",
                "IKZF1 N159Y present" if entry['IKZF1_N159Y'] else "IKZF1 N159Y absent",
                "ZEB2 H1038R present" if entry['ZEB2_H1038R'] else "ZEB2 H1038R absent",
                entry.get('WHO-HAEM5', ''),
                entry.get('ICC', '')
            ])
    
    def _format_fusion_summary(self):
        """Format fusion summary for automatic classification (simplified format)."""
        if not self.data['fusions']:
            return "None: None::None"
        
        # For automatic classification, use simplified format
        fusion = self.data['fusions'][0]  # Take first fusion
        callers = set()
        fusion_name = f"{fusion['gene_1']}::{fusion['gene_2']}"
        
        # Check which callers found this fusion
        for f in self.data['fusions']:
            if f['gene_1'] == fusion['gene_1'] and f['gene_2'] == fusion['gene_2']:
                callers.add(f['caller'])
        
        if len(callers) > 1:
            caller_str = " & ".join(sorted(callers))
        else:
            caller_str = fusion['caller']
        
        return f"{caller_str}: {fusion_name}"
    
    def _format_all_fusions(self):
        """Format all detected fusions for display."""
        if not self.data['fusions']:
            return "No driver fusions identified"
        
        fusion_list = []
        for fusion in self.data['fusions']:
            fusion_list.append(f"{fusion['caller']}: {fusion['gene_1']}::{fusion['gene_2']} (reads: {fusion['spanning_reads']})")
        
        return "; ".join(fusion_list)
    
    def _format_detailed_fusions(self):
        """Format detailed fusion information for manual curation."""
        if not self.data['fusions']:
            return "No driver fusions identified"
        
        detailed_list = []
        for fusion in self.data['fusions']:
            details = [
                f"{fusion['gene_1']}::{fusion['gene_2']}",
                f"Caller: {fusion['caller']}",
                f"Reads: {fusion['spanning_reads']}"
            ]
            
            if fusion.get('breakpoint_1') or fusion.get('breakpoint_2'):
                breakpoints = f"BP1: {fusion.get('breakpoint_1', 'N/A')}, BP2: {fusion.get('breakpoint_2', 'N/A')}"
                details.append(breakpoints)
            
            if fusion.get('sequence_data'):
                seq_data = str(fusion.get('sequence_data'))  # Convert to string first
                details.append(f"Seq: {seq_data[:50]}...")  # Truncate long sequences
            
            detailed_list.append(" | ".join(details))
        
        return "; ".join(detailed_list)
    
    def _generate_driver_fusion_file(self, output_driver):
        """Generate driver fusion output file."""
        if not self.data['fusions']:
            Path(output_driver).touch()
            return
        
        # Convert to DataFrame and save
        fusion_df = pd.DataFrame(self.data['fusions'])
        fusion_df = fusion_df.rename(columns={
            'gene_1': 'Gene_1_symbol(5end_fusion_partner)',
            'gene_2': 'Gene_2_symbol(3end_fusion_partner)',
            'spanning_reads': 'Unique_spanning_reads',
            'caller': 'Fusioncaller'
        })
        
        output_columns = [
            'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)',
            'Unique_spanning_reads',
            'Fusioncaller'
        ]
        
        fusion_df[output_columns].to_csv(output_driver, index=False)


def main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, 
         hotspot_dir, classification_file, output_csv, output_text, output_curation, 
         output_driver):
    """Main processing function with 3-phase approach."""
    
    logger.info(f"🔬 Processing final classification for sample: {sample}")
    logger.info("📋 Starting 3-phase classification process...")
    
    # Initialize processor
    processor = ClassificationProcessor(sample)
    
    # Load classification reference database
    try:
        classification_df = pd.read_csv(classification_file)
        logger.info(f"✅ Loaded classification database: {len(classification_df)} entries")
    except Exception as e:
        logger.error(f"❌ Failed to load classification file: {e}")
        return
    
    # ========== PHASE 1: DATA COLLECTION ==========
    logger.info("📊 PHASE 1: Collecting data...")
    processor.collect_allcatchr_data(allcatchr_file)
    processor.collect_karyotype_data(karyotype_file)
    processor.apply_karyotype_filtering()  # Apply iAMP21 filtering after ALLCatchR is loaded
    processor.collect_snv_data(hotspot_dir)
    processor.collect_fusion_data(fusioncatcher_file, arriba_file, classification_df)
    processor.create_summary_dataframe()
    
    # ========== PHASE 2: SUBTYPE-AWARE MATCHING ==========
    logger.info("🎯 PHASE 2: Attempting subtype-aware matching...")
    results, match_type = processor.find_subtype_aware_matches(classification_df)
    
    # ========== PHASE 3: COMPLEX RULES (if needed) ==========
    if match_type != 'exact_single':
        logger.info("🔧 PHASE 3: Applying complex rules...")
        results, match_type = processor.apply_complex_rules(classification_df, results, match_type)
    
    # ========== OUTPUT GENERATION ==========
    logger.info("📝 Generating output files...")
    output_paths = (output_csv, output_text, output_curation, output_driver)
    processor.generate_output_files(results, match_type, output_paths)
    
    if results is not None and not results.empty:
        logger.info(f"✅ Classification successful: {match_type}")
    else:
        logger.info("⚠️ Manual curation required")
    
    logger.info(f"📁 Output files generated: Final_classification/{sample}_*")


if __name__ == "__main__":
    if len(sys.argv) != 12:
        print("Usage: python make_final_classification_restructured.py <sample> <allcatchr_file> "
              "<karyotype_file> <fusioncatcher_file> <arriba_file> <hotspot_dir> "
              "<classification_file> <output_csv> <output_text> <output_curation> <output_driver>")
        sys.exit(1)
    
    args = sys.argv[1:]
    main(*args)