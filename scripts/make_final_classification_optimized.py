#!/usr/bin/env python3
"""
Optimized IntegrateALL final classification script.
Integrates ALLCatchR predictions, fusion data, hotspot mutations, and karyotype predictions
to generate final WHO-HAEM5/ICC classifications with improved performance.
"""

import csv
import os
import pandas as pd
import numpy as np
import sys
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from functools import partial

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Driver fusion gene list - optimized as set for O(1) lookup
DRIVER_FUSION_GENES = {
    "ABL1", "ABL2", "ACIN1", "ADAMTSL5", "AFF1", "AMPH", "ANTXR1", "ARID1B", "ATF7IP",
    "ATXN7L3", "AUTS2", "BCL2", "BCL6", "BCL9", "BCR", "BMP2K", "BRD9", "C7ORF72",
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
    """Optimized processor for final classification analysis."""
    
    def __init__(self, sample_id):
        self.sample_id = sample_id
        
    def load_karyotype_data(self, karyotype_file):
        """Load karyotype prediction data efficiently."""
        try:
            df = pd.read_csv(karyotype_file)
            if df.empty:
                return "", ""
            return df.iloc[0]["Prediction"], df.iloc[0]["Score"]
        except Exception as e:
            logger.warning(f"Error loading karyotype data: {e}")
            return "", ""
    
    def load_allcatchr_data(self, allcatchr_file):
        """Load ALLCatchR prediction data efficiently."""
        try:
            df = pd.read_csv(allcatchr_file, sep='\t')
            if df.empty:
                return "", "", ""
            row = df.iloc[0]
            return (row["Prediction"], row["Confidence"], 
                   row.get("BCR_ABL1_maincluster_pred", ""))
        except Exception as e:
            logger.warning(f"Error loading ALLCatchR data: {e}")
            return "", "", ""
    
    def check_hotspot_mutations(self, hotspot_dir):
        """Check for relevant hotspot mutation files efficiently."""
        hotspot_patterns = {
            "PAX5_P80R": False,
            "IKZF1_N159Y": False, 
            "ZEB2_H1038R": False
        }
        
        try:
            hotspot_path = Path(hotspot_dir)
            if not hotspot_path.exists():
                return hotspot_patterns
                
            # Use pathlib for efficient file listing
            for file_path in hotspot_path.iterdir():
                filename = file_path.name
                if filename.startswith("PAX5_P80R"):
                    hotspot_patterns["PAX5_P80R"] = True
                elif filename.startswith("IKZF1_N159Y"):
                    hotspot_patterns["IKZF1_N159Y"] = True
                elif filename.startswith("ZEB2_H1038"):
                    hotspot_patterns["ZEB2_H1038R"] = True
                    
        except Exception as e:
            logger.warning(f"Error checking hotspot files: {e}")
            
        return hotspot_patterns
    
    def load_fusion_data(self, fusioncatcher_file, arriba_file):
        """Load and process fusion data from both callers efficiently."""
        fusions = []
        
        # Load FusionCatcher data
        try:
            fc_df = pd.read_csv(fusioncatcher_file, sep='\t', skiprows=1, header=None)
            if not fc_df.empty and len(fc_df.columns) >= 6:
                for _, row in fc_df.iterrows():
                    fusions.append((row[0], row[1], 'FusionCatcher', row[5]))
        except Exception as e:
            logger.warning(f"Error loading FusionCatcher data: {e}")
        
        # Load Arriba data
        try:
            arriba_df = pd.read_csv(arriba_file, sep='\t', skiprows=1, header=None)
            if not arriba_df.empty and len(arriba_df.columns) >= 12:
                for _, row in arriba_df.iterrows():
                    fusions.append((row[0], row[1], 'Arriba', row[11]))
        except Exception as e:
            logger.warning(f"Error loading Arriba data: {e}")
        
        return fusions
    
    def filter_driver_fusions(self, fusions, classification_df, subtype):
        """Filter fusions for driver genes with vectorized operations."""
        if not fusions:
            return []
        
        # Convert to DataFrame for vectorized operations
        fusion_df = pd.DataFrame(fusions, columns=[
            'gene_1', 'gene_2', 'caller', 'spanning_reads'
        ])
        
        # Normalize IGH genes
        fusion_df['gene_1'] = fusion_df['gene_1'].str.replace(r'^IGH.*', 'IGH@', regex=True)
        fusion_df['gene_2'] = fusion_df['gene_2'].str.replace(r'^IGH.*', 'IGH@', regex=True)
        
        # Filter for driver genes using vectorized operations
        driver_mask = (fusion_df['gene_1'].isin(DRIVER_FUSION_GENES) & 
                      fusion_df['gene_2'].isin(DRIVER_FUSION_GENES))
        driver_fusions = fusion_df[driver_mask].copy()
        
        if driver_fusions.empty:
            return []
        
        # Check against classification database using merge operation
        classification_df = classification_df.rename(columns={
            'Gene_1_symbol(5end_fusion_partner)': 'gene_1',
            'Gene_2_symbol(3end_fusion_partner)': 'gene_2'
        })
        
        # Merge to find matching fusions
        matched_fusions = pd.merge(
            driver_fusions, classification_df[['gene_1', 'gene_2']], 
            on=['gene_1', 'gene_2'], how='inner'
        )
        
        # Apply DUX4 filtering logic
        if subtype != 'DUX4':
            dux4_mask = ~((matched_fusions['gene_1'] == 'DUX4') | 
                         (matched_fusions['gene_2'] == 'DUX4'))
            matched_fusions = matched_fusions[dux4_mask]
        
        return matched_fusions.values.tolist()
    
    def create_classification_dataframe(self, allcatchr_data, karyotype_data, 
                                       hotspot_data, filtered_fusions):
        """Create classification DataFrame with optimized structure."""
        subtype, confidence, bcr_abl1 = allcatchr_data
        karyotype, score = karyotype_data
        
        if not filtered_fusions:
            # No fusions case
            return pd.DataFrame([{
                'ALLCatchR': subtype,
                'Ph-pos': bcr_abl1,
                'Confidence': confidence,
                'Gene_1_symbol(5end_fusion_partner)': None,
                'Gene_2_symbol(3end_fusion_partner)': None,
                'Fusioncaller': None,
                'Unique_spanning_reads': None,
                'karyotype_classifier': karyotype,
                'PAX5_P80R': hotspot_data["PAX5_P80R"],
                'IKZF1_N159Y': hotspot_data["IKZF1_N159Y"],
                'ZEB2_H1038R': hotspot_data["ZEB2_H1038R"]
            }])
        else:
            # Multiple fusions case - use pandas DataFrame constructor
            fusion_df = pd.DataFrame(filtered_fusions, columns=[
                'Gene_1_symbol(5end_fusion_partner)',
                'Gene_2_symbol(3end_fusion_partner)',
                'Fusioncaller',
                'Unique_spanning_reads'
            ])
            
            # Add constant columns using vectorized assignment
            fusion_df = fusion_df.assign(
                ALLCatchR=subtype,
                **{'Ph-pos': bcr_abl1},
                Confidence=confidence,
                karyotype_classifier=karyotype,
                PAX5_P80R=hotspot_data["PAX5_P80R"],
                IKZF1_N159Y=hotspot_data["IKZF1_N159Y"],
                ZEB2_H1038R=hotspot_data["ZEB2_H1038R"]
            )
            
            # Reorder columns
            column_order = [
                'ALLCatchR', 'Ph-pos', 'Confidence',
                'Gene_1_symbol(5end_fusion_partner)', 'Gene_2_symbol(3end_fusion_partner)',
                'Fusioncaller', 'Unique_spanning_reads', 'karyotype_classifier',
                'PAX5_P80R', 'IKZF1_N159Y', 'ZEB2_H1038R'
            ]
            
            return fusion_df[column_order]
    
    def match_classifications(self, data_df, classification_df):
        """Match data against classification database using vectorized operations."""
        # Rename columns for consistency
        classification_df = classification_df.rename(columns={
            'karyotype classifier': 'karyotype_classifier',
            'PAX5 P80R': 'PAX5_P80R',
            'IKZF1 N159Y': 'IKZF1_N159Y', 
            'ZEB2 H1038R': 'ZEB2_H1038R'
        })
        
        # Define matching columns based on whether fusions exist
        has_fusions = not data_df['Gene_1_symbol(5end_fusion_partner)'].isnull().all()
        
        if has_fusions:
            merge_cols = [
                'ALLCatchR', 'Ph-pos', 'Gene_1_symbol(5end_fusion_partner)',
                'Gene_2_symbol(3end_fusion_partner)', 'karyotype_classifier',
                'PAX5_P80R', 'IKZF1_N159Y', 'ZEB2_H1038R'
            ]
        else:
            merge_cols = [
                'ALLCatchR', 'Confidence', 'karyotype_classifier',
                'PAX5_P80R', 'IKZF1_N159Y', 'ZEB2_H1038R'
            ]
        
        # Perform vectorized merge operation
        matched_df = pd.merge(
            data_df, classification_df,
            on=merge_cols, how='inner',
            suffixes=('', '_class')
        )
        
        return matched_df if not matched_df.empty else None
    
    def generate_output_files(self, results_df, data_df, output_paths):
        """Generate all output files efficiently."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        if results_df is not None and not results_df.empty:
            self._write_successful_classification(results_df, output_paths)
        else:
            self._write_manual_curation_needed(data_df, output_paths)
        
        # Generate driver fusion file
        self._generate_driver_fusion_file(data_df, output_driver)
    
    def _write_successful_classification(self, results_df, output_paths):
        """Write output for successful automatic classification."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        # Save main results
        results_df.to_csv(output_csv, index=False)
        
        if len(results_df) == 1:
            self._write_single_classification_text(results_df.iloc[0], output_text, output_curation)
        else:
            self._write_multiple_classification_text(results_df, output_text, output_curation)
    
    def _write_manual_curation_needed(self, data_df, output_paths):
        """Write output when manual curation is needed."""
        output_csv, output_text, output_curation, output_driver = output_paths
        
        # Add default classification values
        data_df = data_df.copy()
        data_df['WHO-HAEM5'] = "IntegrateALL couldn't confirm the subtype in concordance with WHO or ICC classification"
        data_df['ICC'] = ""
        
        # Save data
        data_df.to_csv(output_csv, index=False)
        
        # Write curation needed text
        self._write_curation_text(data_df.iloc[0], output_text, output_curation)
    
    def _write_single_classification_text(self, entry, output_text, output_curation):
        """Write text output for single successful classification."""
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
        
        self._write_curation_csv(entry, output_curation, "Automatic Classification")
    
    def _write_multiple_classification_text(self, results_df, output_text, output_curation):
        """Write text output for multiple classifications."""
        # Implementation for multiple entries (simplified)
        entry = results_df.iloc[0]
        self._write_single_classification_text(entry, output_text, output_curation)
    
    def _write_curation_text(self, entry, output_text, output_curation):
        """Write text output when manual curation is needed."""
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
        """Write curation CSV output."""
        with open(output_curation, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            writer.writerow([
                "Classification", "Subtype", "Confidence", "Fusion_details",
                "Karyotype_classifier", "PAX5_P80R", "IKZF1_N159Y",
                "ZEB2_H1038R", "WHO-HAEM5", "ICC"
            ])
            
            fusion_detail = self._format_fusion_detail(entry)
            
            writer.writerow([
                classification_type,
                entry['ALLCatchR'],
                entry['Confidence'],
                fusion_detail,
                entry['karyotype_classifier'],
                "PAX5 P80R present" if entry['PAX5_P80R'] else "PAX5 P80R absent",
                "IKZF1 N159Y present" if entry['IKZF1_N159Y'] else "IKZF1 N159Y absent",
                "ZEB2 H1038R present" if entry['ZEB2_H1038R'] else "ZEB2 H1038R absent",
                entry.get('WHO-HAEM5', ''),
                entry.get('ICC', '')
            ])
    
    def _generate_driver_fusion_file(self, data_df, output_driver):
        """Generate driver fusion output file efficiently."""
        data_df = data_df.copy()
        data_df['Unique_spanning_reads'] = pd.to_numeric(
            data_df['Unique_spanning_reads'], errors='coerce'
        )
        
        # Apply filtering logic
        if 'DUX4' in data_df['ALLCatchR'].values:
            # DUX4 subtype - keep all high-quality fusions
            mask = ~(data_df['Unique_spanning_reads'].isna() | 
                    (data_df['Unique_spanning_reads'] < 3))
        else:
            # Non-DUX4 subtypes - exclude DUX4 fusions
            mask = ~((data_df['Gene_1_symbol(5end_fusion_partner)'] == 'DUX4') |
                    (data_df['Gene_2_symbol(3end_fusion_partner)'] == 'DUX4') |
                    data_df['Unique_spanning_reads'].isna() |
                    (data_df['Unique_spanning_reads'] < 3))
        
        filtered_df = data_df[mask]
        
        if filtered_df.empty:
            # Create empty file
            Path(output_driver).touch()
        else:
            # Save filtered fusion data
            output_columns = [
                'Gene_1_symbol(5end_fusion_partner)',
                'Gene_2_symbol(3end_fusion_partner)',
                'Unique_spanning_reads',
                'Fusioncaller'
            ]
            filtered_df[output_columns].to_csv(output_driver, index=False)


def main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, 
         hotspot_dir, classification_file, output_csv, output_text, output_curation, 
         output_driver):
    """Main processing function with optimized workflow."""
    
    logger.info(f"🔬 Processing final classification for sample: {sample}")
    
    # Initialize processor
    processor = ClassificationProcessor(sample)
    
    # Load classification reference database
    try:
        classification_df = pd.read_csv(classification_file)
        logger.info(f"✅ Loaded classification database: {len(classification_df)} entries")
    except Exception as e:
        logger.error(f"❌ Failed to load classification file: {e}")
        return
    
    # Load all input data efficiently
    logger.info("📊 Loading input data...")
    allcatchr_data = processor.load_allcatchr_data(allcatchr_file)
    karyotype_data = processor.load_karyotype_data(karyotype_file)
    hotspot_data = processor.check_hotspot_mutations(hotspot_dir)
    fusion_data = processor.load_fusion_data(fusioncatcher_file, arriba_file)
    
    logger.info(f"🧬 Found {len(fusion_data)} total fusions")
    
    # Filter for driver fusions
    filtered_fusions = processor.filter_driver_fusions(
        fusion_data, classification_df, allcatchr_data[0]
    )
    
    logger.info(f"🎯 Found {len(filtered_fusions)} driver fusions")
    
    # Create classification DataFrame
    data_df = processor.create_classification_dataframe(
        allcatchr_data, karyotype_data, hotspot_data, filtered_fusions
    )
    
    # Match against classification database
    results_df = processor.match_classifications(data_df, classification_df)
    
    # Generate output files
    output_paths = (output_csv, output_text, output_curation, output_driver)
    processor.generate_output_files(results_df, data_df, output_paths)
    
    if results_df is not None and not results_df.empty:
        logger.info(f"✅ Automatic classification successful: {len(results_df)} matches found")
    else:
        logger.info("⚠️ Manual curation required - no automatic classification match")
    
    logger.info(f"📝 Output files generated in Final_classification/{sample}_*")


if __name__ == "__main__":
    if len(sys.argv) != 12:
        print("Usage: python make_final_classification_optimized.py <sample> <allcatchr_file> "
              "<karyotype_file> <fusioncatcher_file> <arriba_file> <hotspot_dir> "
              "<classification_file> <output_csv> <output_text> <output_curation> <output_driver>")
        sys.exit(1)
    
    args = sys.argv[1:]
    main(*args)