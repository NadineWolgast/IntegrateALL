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
                for _, row in fc_df.iterrows():
                    gene1, gene2, spanning_reads = row[0], row[1], row[5]
                    if self._is_driver_fusion_in_classtest(gene1, gene2, classification_df):
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
        # Normalize IGH genes
        gene1_norm = 'IGH@' if gene1.startswith('IGH') else gene1
        gene2_norm = 'IGH@' if gene2.startswith('IGH') else gene2
        
        # Check if both genes are driver genes
        if not (gene1_norm in DRIVER_FUSION_GENES and gene2_norm in DRIVER_FUSION_GENES):
            return False
        
        # Check if combination exists in Class_test.csv (both orientations)
        exists = (
            ((classification_df['Gene_1_symbol(5end_fusion_partner)'] == gene1_norm) &
             (classification_df['Gene_2_symbol(3end_fusion_partner)'] == gene2_norm)).any() or
            ((classification_df['Gene_1_symbol(5end_fusion_partner)'] == gene2_norm) &
             (classification_df['Gene_2_symbol(3end_fusion_partner)'] == gene1_norm)).any()
        )
        
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
    
    # ========== PHASE 2: SIMPLE MATCHING ==========
    
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
        
        # For non-CEBP subtypes with ZEB2=False in data, allow matching with ZEB2=True in class_test
        if subtype != 'CEBP' and not summary_for_match['ZEB2_H1038R'].iloc[0]:
            logger.info("🔧 ZEB2=False in data, allowing match with ZEB2=True in classification DB")
        
        # Perform exact merge on selected columns
        matches = pd.merge(
            summary_for_match, classification_for_match,
            on=match_cols, how='inner',
            suffixes=('', '_class')
        )
        
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
            return self._resolve_multiple_matches(matches)
        elif match_type == 'no_exact':
            return self._resolve_no_matches(classification_df)
        else:
            return matches
    
    def _resolve_multiple_matches(self, matches):
        """Resolve multiple exact matches."""
        logger.info("🔧 Resolving multiple matches...")
        
        # For now, take the first match
        # TODO: Implement more sophisticated rules
        return matches.head(1)
    
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
                        if spanning_reads >= 3:
                            valid_fusions.append(fusion)
                            logger.info(f"✅ DUX4 fusion accepted: {fusion['gene_1']}::{fusion['gene_2']} (reads: {spanning_reads})")
                        else:
                            logger.info(f"❌ DUX4 fusion rejected: {fusion['gene_1']}::{fusion['gene_2']} (reads: {spanning_reads} < 3)")
                    except:
                        logger.info(f"❌ DUX4 fusion rejected: {fusion['gene_1']}::{fusion['gene_2']} (invalid read count)")
                else:
                    valid_fusions.append(fusion)
            
            # Update fusions and recreate summary
            self.data['fusions'] = valid_fusions
            self.create_summary_dataframe()
            
            # Try exact matching again with filtered data
            return self.find_exact_matches(classification_df)
        
        return None, 'no_match_after_rules'
    
    # ========== OUTPUT GENERATION ==========
    
    def generate_output_files(self, results, match_type, output_paths):
        """Generate all output files."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        if results is not None and not results.empty:
            self._write_successful_classification(results, output_paths)
        else:
            self._write_manual_curation_needed(output_paths)
        
        # Generate driver fusion file
        self._generate_driver_fusion_file(output_driver)
    
    def _write_successful_classification(self, results_df, output_paths):
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
        """Write enhanced curation CSV output with multiple fusions and extended hotspots."""
        with open(output_curation, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Enhanced header with additional columns
            writer.writerow([
                "Classification", "Subtype", "Confidence", "All_Fusions", "Fusion_Details",
                "Karyotype_classifier", "PAX5_P80R", "IKZF1_N159Y", "ZEB2_H1038R",
                "Extended_Hotspots", "WHO-HAEM5", "ICC"
            ])
            
            # Format all fusions (not just the first one)
            all_fusions = self._format_all_fusions()
            fusion_details = self._format_detailed_fusions() if classification_type == "Manual Curation" else ""
            
            # Format extended hotspots
            extended_hotspots = "; ".join(self.data.get('extended_hotspots', []))
            
            writer.writerow([
                classification_type,
                entry['ALLCatchR'],
                entry['Confidence'],
                all_fusions,
                fusion_details,
                entry['karyotype_classifier'],
                "PAX5 P80R present" if entry['PAX5_P80R'] else "PAX5 P80R absent",
                "IKZF1 N159Y present" if entry['IKZF1_N159Y'] else "IKZF1 N159Y absent",
                "ZEB2 H1038R present" if entry['ZEB2_H1038R'] else "ZEB2 H1038R absent",
                extended_hotspots if extended_hotspots else "None detected",
                entry.get('WHO-HAEM5', ''),
                entry.get('ICC', '')
            ])
    
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
                details.append(f"Seq: {fusion.get('sequence_data')[:50]}...")  # Truncate long sequences
            
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
    
    # ========== PHASE 2: SIMPLE MATCHING ==========
    logger.info("🎯 PHASE 2: Attempting exact matching...")
    results, match_type = processor.find_exact_matches(classification_df)
    
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