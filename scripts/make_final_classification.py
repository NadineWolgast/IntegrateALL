import csv
import os
import pandas as pd

no_driver_subtypes_list = ["Hyperdiploid", "Low hypodiploid", "Near haploid", "iAMP21", "IKZF1 N159Y", "PAX5 P80R",
                           "PAX5alt", "Ph-like", "ETV6::RUNX1-like", "ZNF384", "KMT2A", "CEBP"]

no_driver_karyotype_list = ["Hyperdiploid", "Low hypodiploid", "Near haploid", "iAMP21", "other"]

driver_fusion_list = [
    "ABL1",
    "ABL2",
    "ACIN1",
    "ADAMTSL5",
    "AFF1",
    "AMPH",
    "ANTXR1",
    "ARID1B",
    "ATF7IP",
    "ATXN7L3",
    "AUTS2",
    "BCL2",
    "BCL6",
    "BCL9",
    "BCR",
    "BCR",
    "BMP2K",
    "BRD9",
    "C7ORF72",
    "CASC15",
    "CBFA2T2",
    "CBFA2T3",
    "CBL",
    "CD163",
    "CEBPA",
    "CEBPB",
    "CEBPD",
    "CEBPE",
    "CENPC",
    "CLTC",
    "CREBBP",
    "CRLF2",
    "CSF1R",
    "CUX1",
    "DACH1",
    "DACH2",
    "DAZAP1",
    "DBX1",
    "DCPS",
    "DDX42",
    "DGKH",
    "DMRTA2",
    "DUX4",
    "EBF1",
    "ELMO1",
    "ELN",
    "EP300",
    "EPOR",
    "EPS15",
    "ERC1",
    "ESRRA",
    "ETV6",
    "EWSR1",
    "EXOSC2",
    "EXTL1",
    "FAM136A",
    "FBRSL1",
    "FIP1L1",
    "FKBP15",
    "FOXJ2",
    "FOXO3",
    "FOXP2",
    "GATA2DA",
    "GATAD2A",
    "GOPC",
    "HLF",
    "HNRNPM",
    "HNRNPUL1",
    "HOXA9",
    "ID4",
    "IGH@",
    "IGK",
    "IKZF1",
    "IL2RB",
    "IL7-R",
    "JAK2",
    "KANK1",
    "KMT2A",
    "LAIR1",
    "LEF1",
    "LSM14A",
    "LYN",
    "MBNL1",
    "MED12",
    "MEF2D",
    "MEIS2",
    "MLLT1",
    "MLLT10",
    "MLLT3",
    "MPRIP",
    "MYB",
    "MYC",
    "MYH9",
    "MYO18B",
    "NCOA5",
    "NCOR1",
    "NIPBL",
    "NOL4L",
    "NTRK3",
    "NUP153",
    "NUP214",
    "NUTM1",
    "OFD1",
    "P2RY8",
    "PAG1",
    "PAX5",
    "PBX1",
    "PCGF5",
    "PCM1",
    "PDGFRA",
    "PDGFRB",
    "PML",
    "PPFIBP1",
    "PTK2B",
    "PYGO2",
    "QSOX1",
    "RANBP2",
    "RCSD1",
    "RFX3",
    "RHOXF2B",
    "RNFT2",
    "ROS1",
    "RUNX1",
    "SFPQ",
    "SLC12A6",
    "SLC30A7",
    "SMARCA2",
    "SNX1",
    "SNX2",
    "SNX29",
    "SRRM1",
    "SS18",
    "SSBP2",
    "STIM2",
    "STRN3",
    "TAF15",
    "TAF3",
    "TBL1XR1",
    "TCF3",
    "TCF4",
    "TERF2",
    "THADA",
    "TMEM2",
    "TMPRSS9",
    "TMTC1",
    "TNIP1",
    "TNS3",
    "TPR",
    "TYK2",
    "UBASH3B",
    "UBTF",
    "USP2",
    "USP25",
    "WDR37",
    "WDR5",
    "ZC3HAV1",
    "ZEB2",
    "ZFAND3",
    "ZMIZ1",
    "ZMYND8",
    "ZNF274",
    "ZNF276",
    "ZNF340",
    "ZNF362",
    "ZNF384",
    "ZNF521",
    "ZNF618",
    "ZPBP"
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

    for gene_1, gene_2, caller, unique_spanning_reads in fusion_genes:
        if gene_1.startswith("IGH"):
            gene_1 = "IGH@"
        if gene_2.startswith("IGH"):
            gene_2 = "IGH@"
        if gene_1 in unique_genes and gene_2 in unique_genes and int(unique_spanning_reads) >= 3:
        #if gene_1 in unique_genes and gene_2 in unique_genes:
            print(gene_1, gene_2, )
            # Check if the pair exists in the df
            match_df = df[
                (df['Gene_1_symbol(5end_fusion_partner)'] == gene_1) &
                (df['Gene_2_symbol(3end_fusion_partner)'] == gene_2)
                ]
            if not match_df.empty:
                filtered_fusions.append((gene_1, gene_2, caller, unique_spanning_reads))
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


def gather_data(allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir):
    subgruppe, confidence, bcr_abl1_maincluster_pred = get_allcatchr_data(allcatchr_file)
    karyotype, score = get_karyotype_and_probabilities(karyotype_file)
    relevant_files = check_hotspot_files(hotspot_dir)

    fusion_genes = []
    with open(fusioncatcher_file, 'r') as f:
        next(f)  # Überspringen der Header-Zeile
        for line in f:
            parts = line.strip().split('\t')
            fusion_genes.append((parts[0], parts[1], 'fusioncatcher', parts[5])) 

    with open(arriba_file, 'r') as f:
        next(f)  # Überspringen der Header-Zeile
        for line in f:
            parts = line.strip().split('\t')
            fusion_genes.append((parts[0], parts[1], 'arriba', parts[11]))

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

    if subgruppe and karyotype and not filtered_fusions:
        df_filtered = df[
            (df['Gene_1_symbol(5end_fusion_partner)'].isna()) & (df['Gene_2_symbol(3end_fusion_partner)'].isna())]

        for index, row in df_filtered.iterrows():
            if (row['ALLCatchR'] == filtered_subtype and row['karyotype classifier'] == filtered_karyotype
                    and row['Ph-pos'] == bcr_abl1_maincluster_pred):
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


def create_fusion_table(filtered_fusions):
    table_data = []

    for fusion in filtered_fusions:
        gene1, gene2, fusioncaller, unique_spanning_reads = fusion
        table_data.append({
            'Gene1': gene1,
            'Gene2': gene2,
            'Fusioncaller': fusioncaller,
            'unique_spanning_reads': unique_spanning_reads
        })

    fusion_table = pd.DataFrame(table_data)

    return fusion_table

def main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir, classification_file,
         output_csv, output_text, output_driver):
    subgruppe, confidence, bcr_abl1_maincluster_pred, karyotype, relevant_files, fusion_genes = gather_data(
        allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir)
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

    def get_filtered_fusion_genes(filtered_fusions, matched_fusions):
        filtered_gene_pairs = set()
        for fusion in filtered_fusions:
            gene_1, gene_2 = fusion[0], fusion[1]
            if not pd.isna(gene_1) and not pd.isna(gene_2):
                pair = (gene_1, gene_2)
                if pair not in matched_fusions:
                    filtered_gene_pairs.add(pair)
        return list(filtered_gene_pairs)

    matched_fusions = set()

    # Ergebnisse in die Ausgabedatei schreiben
    if results:
        output_df = pd.DataFrame(results)
        output_df.to_csv(output_csv, index=False)
        print(f"Results written to {output_csv}")

        with open(output_text, 'w') as f:
            for result in results:
                fusion_text = format_fusion(result)
                if fusion_text != "no fusion":
                    gene_1, gene_2 = result['Gene_1_symbol(5end_fusion_partner)'], result['Gene_2_symbol(3end_fusion_partner)']
                    matched_fusions.add((gene_1, gene_2))
                f.write(f"Consistency between gene expression-based subtype-allocation (ALLCatchR: {subgruppe} subtype, {confidence} confidence) "
                        f"and the genomic driver profile (fusioncatcher/arriba: {fusion_text}, "
                        f"RNA-Seq CNV karyotype classifier: {karyotype}, subtype defining SNPs: "
                        f"{'PAX5 P80R present' if relevant_files['PAX5_P80R'] else 'absent'}, "
                        f"{'IKZF1 N159Y present' if relevant_files['IKZF1_N159Y'] else 'absent'}) "
                        f"supports a classification as\n\n"
                        f"{result['WHO-HAEM5']} according to WHO-HAEM5 (Alaggio R et al. Leukemia, 2022) and\n"
                        f"{result['ICC']} according to ICC (Arber D et al. Blood, 2022).\n\n")

            # Wenn es mehr als ein Fusionspaar gibt, schreiben Sie die gefilterten Fusionsgene, die keinen Match hatten
            if len(filtered_fusions) > 1:
                fusion_table = create_fusion_table(filtered_fusions)
                fusion_table.to_csv(output_driver, index=False)
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
                    f"and the genomic driver profile (fusioncatcher/arriba: no fusion, "
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
    arriba_file = sys.argv[5]
    hotspot_dir = sys.argv[6]
    classification_file = sys.argv[7]
    output_csv = sys.argv[8]
    output_driver = sys.argv[10]
    main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir, classification_file,
         output_csv, output_text, output_driver)

