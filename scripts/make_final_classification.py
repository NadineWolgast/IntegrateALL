import csv
import os
import pandas as pd

no_driver_subtypes_list = ["Hyperdiploid", "Low hypodiploid", "Near haploid", "iAMP21", "IKZF1 N159Y", "PAX5 P80R",
                           "PAX5alt", "Ph-like", "ETV6::RUNX1-like", "ZNF384", "KMT2A", "CEBP"]

no_driver_karyotype_list = ["Hyperdiploid", "Low hypodiploid", "Near haploid", "iAMP21", "other"]

driver_fusion_list = [
    "BCR",
    "ABL1",
    "ABL2",
    "RCSD1",
    "ZC3HAV1",
    "EBF1",
    "PDGFRB",
    "NUP214",
    "CSF1R",
    "SSBP2",
    "ETV6",
    "ZMIZ1",
    "RANBP2",
    "ATF7IP",
    "SNX2",
    "PAG1",
    "MEF2D",
    "ZEB2",
    "CENPC",
    "LSM14A",
    "TBL1XR1",
    "FIP1L1",
    "GATAD2A",
    "LYN",
    "NCOR1",
    "EXOSC2",
    "FOXP1",
    "MYO18B",
    "NUP153",
    "SFPQ",
    "SNX1",
    "ZMYND8",
    "EPOR",
    "IGK",
    "LAIR1",
    "PCM1",
    "PPFIBP1",
    "TERF2",
    "TPR",
    "JAK2",
    "IL2RB",
    "MYH9",
    "IL7-R",
    "OFD1",
    "STRN3",
    "USP25",
    "ZNF274",
    "MYB",
    "TYK2",
    "RFX3",
    "SMU1",
    "WDR37",
    "ZNF340",
    "NTRK3",
    "GOPC",
    "ROS1",
    "PTK2B",
    "TMEM2",
    "CBL",
    "KANK1",
    "DGKH",
    "ZFAND3",
    "IKZF1",
    "AFF1",
    "KMT2A",
    "MLLT1",
    "MLLT3",
    "MLLT10",
    "USP2",
    "DCPS",
    "EPS15",
    "TNS3",
    "UBASH3B",
    "RUNX1",
    "ELMO1",
    "AMPH",
    "C7ORF72",
    "CASC15",
    "CD163",
    "EBF1",
    "ERC1",
    "EXTL1",
    "FAM136A",
    "FOXO3",
    "QSOX1",
    "SLC30A7",
    "SRRM1",
    "TMTC1",
    "RNFT2",
    "ZPBP",
    "TCF3",
    "PBX1",
    "HLF",
    "TCF4",
    "BCL2",
    "IGH@",
    "MYC",
    "DUX4",
    "PCGF5",
    "BCL9",
    "HNRNPUL1",
    "PYGO2",
    "DAZAP1",
    "FOXJ2",
    "HNRNPM",
    "SS18",
    "EP300",
    "ZNF384",
    "TAF15",
    "TCF3",
    "EWSR1",
    "SMARCA2",
    "ARID1B",
    "CLTC",
    "CREBBP",
    "DDX42",
    "NIPBL",
    "ZNF362",
    "NUTM1",
    "SLC12A6",
    "ACIN1",
    "CUX1",
    "BRD9",
    "ZNF618",
    "UBTF",
    "ATXN7L3",
    "PAX5",
    "NOL4L",
    "AUTS2",
    "ZNF521",
    "CBFA2T3",
    "DACH1",
    "CBFA2T2",
    "NCOA5",
    "ADAMTSL5",
    "ANTXR1",
    "BMP2K",
    "LEF1",
    "DACH2",
    "DBX1",
    "DMRTA2",
    "ESRRA",
    "FKBP15",
    "FOXP2",
    "ID4",
    "MBNL1",
    "MEIS2",
    "MPRIP",
    "PML",
    "RHOXF2B",
    "TAF3",
    "TMPRSS9",
    "WDR5",
    "ZNF276",
    "CEBPA",
    "CEBPE",
    "CEBPB",
    "CEBPD",
    "HOXA9",
    "MED12"
]


def get_karyotype_and_probabilities(karyotype_file):
    with open(karyotype_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            return row["Prediction"], row["Score"]
    return "", ""


def get_allcatchr_data(allcatchr_file):
    with open(allcatchr_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        allcatchr_data = next(reader)
        return allcatchr_data["Prediction"], allcatchr_data["Confidence"], allcatchr_data["BCR_ABL1_maincluster_pred"]


def check_hotspot_files(hotspot_dir):
    hotspot_files = os.listdir(hotspot_dir)
    relevant_files = {
        "PAX5_P80R": False,
        "IKZF1_N159Y": False,
        "ZEB2_H1038R": False
    }

    for file in hotspot_files:
        if file.startswith("PAX5_P80R"):
            relevant_files["PAX5_P80R"] = True
        elif file.startswith("IKZF1_N159Y"):
            relevant_files["IKZF1_N159Y"] = True
        elif file.startswith("ZEB2_H1038R"):
            relevant_files["ZEB2_H1038R"] = True

    return relevant_files


def filter_fusions(fusion_genes, unique_genes, df):
    filtered_fusions = []

    for gene_1, gene_2 in fusion_genes:
        if gene_1.startswith("IGH"):
            gene_1 = "IGH@"
        if gene_2.startswith("IGH"):
            gene_2 = "IGH@"
        if gene_1 in unique_genes and gene_2 in unique_genes:
            # Check if the pair exists in the df
            match_df = df[
                (df['Gene_1_symbol(5end_fusion_partner)'] == gene_1) &
                (df['Gene_2_symbol(3end_fusion_partner)'] == gene_2)
            ]
            if not match_df.empty:
                filtered_fusions.append((gene_1, gene_2))
    return filtered_fusions


def filter_subtype(subgruppe, unique_subtypes):
    filtered_subtype = []
    if subgruppe in unique_subtypes:
        filtered_subtype.append(subgruppe)

    return filtered_subtype


def filter_karyotype(karyotype, unique_karyotypes):
    filtered_karyotype = []
    if karyotype in unique_karyotypes:
        filtered_karyotype.append(karyotype)

    return filtered_karyotype


def gather_data(allcatchr_file, karyotype_file, fusioncatcher_file, hotspot_dir):
    subgruppe, confidence, bcr_abl1_maincluster_pred = get_allcatchr_data(allcatchr_file)
    karyotype, score = get_karyotype_and_probabilities(karyotype_file)
    relevant_files = check_hotspot_files(hotspot_dir)

    fusion_genes = []
    with open(fusioncatcher_file, 'r') as f:
        next(f)  # Überspringen der Header-Zeile
        for line in f:
            parts = line.strip().split('\t')
            fusion_genes.append((parts[0], parts[1]))

    return subgruppe, confidence, bcr_abl1_maincluster_pred, karyotype, relevant_files, fusion_genes


def all_false_first_two(relevant_files):
    keys = list(relevant_files.keys())
    return relevant_files[keys[0]] is False and relevant_files[keys[1]] is False


def check_conditions(subgruppe, confidence, bcr_abl1_maincluster_pred, karyotype, relevant_files, fusion_genes, df):
    results = []
    filtered_fusions = filter_fusions(fusion_genes, driver_fusion_list, df)
    print("filtered_fusions: ", filtered_fusions, "subgruppe: ", subgruppe, "confidence: ", confidence,
          "bcr_abl1_maincluster_pred: ", bcr_abl1_maincluster_pred, "karyotype: ", karyotype,
          "relevant_files: ", relevant_files)
    filtered_subtype = subgruppe if subgruppe in df['ALLCatchR'].values else None
    filtered_karyotype = karyotype if karyotype in df['karyotype classifier'].values else None
    #print("filtered_subtype", filtered_subtype)
    #print("filtered_karyotype", filtered_karyotype)
    #print("bcr_abl1_maincluster_pred", bcr_abl1_maincluster_pred)

    if filtered_subtype and filtered_karyotype and not filtered_fusions:
        print("Keine fusionen")
        df_filtered = df[
            (df['Gene_1_symbol(5end_fusion_partner)'].isna()) & (df['Gene_2_symbol(3end_fusion_partner)'].isna())]

        for index, row in df_filtered.iterrows():
            #print("row", row)
            if (row['ALLCatchR'] == filtered_subtype and row['karyotype classifier'] == filtered_karyotype
                    and row['Ph-pos'] == bcr_abl1_maincluster_pred):
                print(row['ALLCatchR'], filtered_subtype, row['karyotype classifier'], filtered_karyotype)
                if (filtered_subtype == "PAX5 P80R" and relevant_files["PAX5_P80R"]
                        and not relevant_files["IKZF1_N159Y"]):
                    row['Confidence'] = confidence
                    row['Ph-pos'] = bcr_abl1_maincluster_pred
                    row['PAX5_P80R'] = relevant_files["PAX5_P80R"]
                    row['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
                    row['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]
                    results.append(row)
                elif subgruppe == "IKZF1 N159Y" and relevant_files["IKZF1_N159Y"] and not relevant_files["PAX5_P80R"]:
                    row['Confidence'] = confidence
                    row['Ph-pos'] = bcr_abl1_maincluster_pred
                    row['PAX5_P80R'] = relevant_files["PAX5_P80R"]
                    row['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
                    row['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]
                    results.append(row)
                elif (subgruppe == "CEBP" and relevant_files["ZEB2_H1038R"] and not relevant_files["PAX5_P80R"]
                      and not relevant_files["IKZF1_N159Y"] and confidence == row['Confidence']):
                        row['Confidence'] = confidence
                        row['Ph-pos'] = bcr_abl1_maincluster_pred
                        row['PAX5_P80R'] = relevant_files["PAX5_P80R"]
                        row['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
                        row['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]
                        results.append(row)
                elif (filtered_subtype in ["PAX5alt", "Ph-like", "ETV6::RUNX1-like", "ZNF384", "KMT2A", "DUX4", "BCL2/MYC"] and
                      confidence == "high-confidence" and not relevant_files["PAX5_P80R"] and not relevant_files["IKZF1_N159Y"]):
                    row['Confidence'] = confidence
                    row['Ph-pos'] = bcr_abl1_maincluster_pred
                    row['PAX5_P80R'] = relevant_files["PAX5_P80R"]
                    row['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
                    row['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]
                    results.append(row)
                elif (filtered_subtype in ["Hyperdiploid", "Low hypodiploid", "Near haploid", "iAMP21", "KMT2A"] and
                      not relevant_files["PAX5_P80R"] and not relevant_files["IKZF1_N159Y"]):
                    row['Confidence'] = confidence
                    row['Ph-pos'] = bcr_abl1_maincluster_pred
                    row['PAX5_P80R'] = relevant_files["PAX5_P80R"]
                    row['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
                    row['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]
                    results.append(row)
                else:
                    print("no match")
            else:
                print("Different subtype and karyotype")
    else:
        print("Mit fusionen")
        for index, row in df.iterrows():
            all_conditions_met = True
            gene_1 = row['Gene_1_symbol(5end_fusion_partner)']
            gene_2 = row['Gene_2_symbol(3end_fusion_partner)']
            if gene_1 and gene_2:
                match_found = any(
                    (fusion_gene_1 == gene_1 and fusion_gene_2 == gene_2) for fusion_gene_1, fusion_gene_2 in
                    filtered_fusions
                )
                if not match_found:
                    all_conditions_met = False
            if row['ALLCatchR'] != subgruppe or row['karyotype classifier'] != karyotype:
                all_conditions_met = False
            if relevant_files["PAX5_P80R"] and row['PAX5 P80R'] != 'yes':
                all_conditions_met = False
            if relevant_files["IKZF1_N159Y"] and row['IKZF1 N159Y'] != 'yes':
                all_conditions_met = False
            if row['Ph-pos'] != bcr_abl1_maincluster_pred:
                all_conditions_met = False
            if all_conditions_met:
                row['Confidence'] = confidence
                row['Ph-pos'] = bcr_abl1_maincluster_pred
                row['PAX5 P80R'] = relevant_files["PAX5_P80R"]
                row['IKZF1 N159Y'] = relevant_files["IKZF1_N159Y"]
                row['ZEB2 H1038R'] = relevant_files["ZEB2_H1038R"]
                results.append(row)
    return results



def main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, hotspot_dir, classification_file, output_csv,
         output_text):
    subgruppe, confidence, bcr_abl1_maincluster_pred, karyotype, relevant_files, fusion_genes = gather_data(
        allcatchr_file, karyotype_file, fusioncatcher_file, hotspot_dir)
    df = pd.read_csv(classification_file, sep='\t')
    filtered_fusions = filter_fusions(fusion_genes, driver_fusion_list, df)

    results = check_conditions(subgruppe, confidence, bcr_abl1_maincluster_pred, karyotype, relevant_files,
                               fusion_genes, df)

    def format_fusion(result):
        gene_1 = result['Gene_1_symbol(5end_fusion_partner)']
        gene_2 = result['Gene_2_symbol(3end_fusion_partner)']
        if pd.isna(gene_1) or pd.isna(gene_2):
            return "no fusion"
        return f"{gene_1}::{gene_2}"

    # Ergebnisse in die Ausgabedatei schreiben
    if results:
        output_df = pd.DataFrame(results)
        output_df.to_csv(output_csv, index=False)
        print(f"Results written to {output_csv}")

        with open(output_text, 'w') as f:
            for result in results:
                fusion_text = format_fusion(result)
                f.write(f"Consistency between gene expression-based subtype-allocation (ALLCatchR: {subgruppe} subtype, {confidence} confidence) "
                        f"and the genomic driver profile (fusioncatcher: {fusion_text}, "
                        f"RNA-Seq CNV karyotype classifier: {karyotype}, subtype defining SNPs: "
                        f"{'PAX5 P80R present' if relevant_files['PAX5_P80R'] else 'absent'}, "
                        f"{'IKZF1 N159Y present' if relevant_files['IKZF1_N159Y'] else 'absent'}) "
                        f"supports a classification as\n\n"
                        f"{result['WHO-HAEM5']} according to WHO-HAEM5 (Alaggio R et al. Leukemia, 2022) and\n"
                        f"{result['ICC']} according to ICC (Arber D et al. Blood, 2022).\n\n")
        print(f"Detailed results written to {output_text}")

    else:
        with open(output_csv, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(["Message"])
            writer.writerow(
                ["IntegrateALL couldn't confirm the subtype in concordance with WHO or ICC classification"])
        print(f"No matches found. Message written to {output_csv}")

        with open(output_text, 'w') as f:
            f.write(f"IntegretALL classification:\n\n"
                    f"Gene expression-based subtype-allocation (ALLCatchR: {subgruppe} subtype, {confidence} confidence) "
                    f"and the genomic driver profile (fusioncatcher: no fusion, "
                    f"RNA-Seq CNV karyotype classifier: {karyotype}, subtype defining SNPs: "
                    f"{'PAX5 P80R present' if relevant_files['PAX5_P80R'] else 'absent'}, "
                    f"{'IKZF1 N159Y present' if relevant_files['IKZF1_N159Y'] else 'absent'}) "
                    f"seem not to be consistent with an unambiguous diagnostic classification according to WHO-HAEM5 "
                    f"(Alaggio R et al. Leukemia, 2022) / ICC (Arber D et al. Blood, 2022).")
        print(f"Detailed message written to {output_text}")


if __name__ == "__main__":
    import sys

    sample = sys.argv[1]
    allcatchr_file = sys.argv[2]
    karyotype_file = sys.argv[3]
    fusioncatcher_file = sys.argv[4]
    hotspot_dir = sys.argv[5]
    classification_file = sys.argv[6]
    output_csv = sys.argv[7]
    output_text = sys.argv[8]
    main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, hotspot_dir, classification_file, output_csv,
         output_text)
