import csv
import os
import pandas as pd

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
        elif file.startswith("ZEB2_H1038"):
            relevant_files["ZEB2_H1038R"] = True

    return relevant_files


def filter_fusions(fusion_genes, unique_genes, df, subgruppe):
    filtered_fusions = []
    for gene_1, gene_2, caller, unique_spanning_reads in fusion_genes:
        # Anpassen von IGH
        if gene_1.startswith("IGH"):
            gene_1 = "IGH@"
        if gene_2.startswith("IGH"):
            gene_2 = "IGH@"

        # Überprüfen, ob beide Gene in unique_genes enthalten sind
        if gene_1 in unique_genes and gene_2 in unique_genes:
            # Überprüfen, ob das Gen-Paar im Classification DF vorhanden ist
            match_df = df[
                (df['Gene_1_symbol(5end_fusion_partner)'] == gene_1) &
                (df['Gene_2_symbol(3end_fusion_partner)'] == gene_2)
                ]
            if not match_df.empty and (subgruppe == 'DUX4' or ('DUX4' not in [gene_1, gene_2])):
                filtered_fusions.append((gene_1, gene_2, caller, unique_spanning_reads))

    return filtered_fusions


def gather_data(allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir, classification_file):
    # Extrahiere die notwendigen Informationen aus den Eingabedateien
    subgruppe, confidence, bcr_abl1_maincluster_pred = get_allcatchr_data(allcatchr_file)
    karyotype, score = get_karyotype_and_probabilities(karyotype_file)
    relevant_files = check_hotspot_files(hotspot_dir)

    fusion_genes = []
    # Lese Fusionen aus FusionCatcher-Datei ein
    with open(fusioncatcher_file, 'r') as f:
        next(f)  # Überspringen der Header-Zeile
        for line in f:
            parts = line.strip().split('\t')
            fusion_genes.append((parts[0], parts[1], 'FusionCatcher', parts[5]))

    # Lese Fusionen aus Arriba-Datei ein
    with open(arriba_file, 'r') as f:
        next(f)  # Überspringen der Header-Zeile
        for line in f:
            parts = line.strip().split('\t')
            fusion_genes.append((parts[0], parts[1], 'Arriba', parts[11]))

    # Lese die Klassifikationsdatei ein
    df = pd.read_csv(classification_file, sep=',')  # Hier ',' als Delimiter verwenden

    # Filtere die Fusionen
    filtered_fusions = filter_fusions(fusion_genes, driver_fusion_list, df, subgruppe)

    if not filtered_fusions:
        # Erstelle einen DataFrame mit einer Zeile und allen notwendigen Werten
        result_df = pd.DataFrame([{
            'ALLCatchR': subgruppe,
            'Ph-pos': bcr_abl1_maincluster_pred,
            'Confidence': confidence,
            'Gene_1_symbol(5end_fusion_partner)': None,
            'Gene_2_symbol(3end_fusion_partner)': None,
            'Fusioncaller': None,
            'Unique_spanning_reads': None,
            'karyotype_classifier': karyotype,
            'PAX5_P80R': relevant_files["PAX5_P80R"],
            'IKZF1_N159Y': relevant_files["IKZF1_N159Y"],
            'ZEB2_H1038R': relevant_files["ZEB2_H1038R"]
        }])
    else:
        # Erstellen des DataFrames für den Rückgabewert
        result_df = pd.DataFrame(filtered_fusions, columns=[
            'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)',
            'Fusioncaller',
            'Unique_spanning_reads'
        ])

        # Hinzufügen der konstanten Spalten
        result_df['ALLCatchR'] = subgruppe
        result_df['Ph-pos'] = bcr_abl1_maincluster_pred
        result_df['Confidence'] = confidence
        result_df['karyotype_classifier'] = karyotype
        result_df['PAX5_P80R'] = relevant_files["PAX5_P80R"]
        result_df['IKZF1_N159Y'] = relevant_files["IKZF1_N159Y"]
        result_df['ZEB2_H1038R'] = relevant_files["ZEB2_H1038R"]

        # Spalten neu anordnen
        result_df = result_df[[
            'ALLCatchR',
            'Ph-pos',
            'Confidence',
            'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)',
            'Fusioncaller',
            'Unique_spanning_reads',
            'karyotype_classifier',
            'PAX5_P80R',
            'IKZF1_N159Y',
            'ZEB2_H1038R'
        ]]

    return result_df


import pandas as pd

def check_conditions(data_df, df_classification):
    comparison_df = df_classification
    matched_rows = []

    if not data_df['Gene_1_symbol(5end_fusion_partner)'].isnull().all():
        for _, data_row in data_df.iterrows():
            for _, class_row in comparison_df.iterrows():
                match_found = (
                    data_row['ALLCatchR'] == class_row['ALLCatchR'] and
                    data_row['Ph-pos'] == class_row['Ph-pos'] and
                    data_row['Gene_1_symbol(5end_fusion_partner)'] == class_row['Gene_1_symbol(5end_fusion_partner)'] and
                    data_row['Gene_2_symbol(3end_fusion_partner)'] == class_row['Gene_2_symbol(3end_fusion_partner)'] and
                    data_row['karyotype_classifier'] == class_row['karyotype classifier'] and
                    data_row['PAX5_P80R'] == class_row['PAX5 P80R'] and
                    data_row['IKZF1_N159Y'] == class_row['IKZF1 N159Y'] and
                    data_row['ZEB2_H1038R'] == class_row['ZEB2 H1038R']
                )

                if match_found:
                    data_row['WHO-HAEM5'] = class_row['WHO-HAEM5']
                    data_row['ICC'] = class_row['ICC']
                    matched_rows.append(data_row)

    else:
        for _, data_row in data_df.iterrows():
            for _, class_row in comparison_df.iterrows():
                match_found = (
                    data_row['ALLCatchR'] == class_row['ALLCatchR'] and
                    data_row['Confidence'] == class_row['Confidence'] and
                    data_row['karyotype_classifier'] == class_row['karyotype classifier'] and
                    data_row['PAX5_P80R'] == class_row['PAX5 P80R'] and
                    data_row['IKZF1_N159Y'] == class_row['IKZF1 N159Y'] and
                    data_row['ZEB2_H1038R'] == class_row['ZEB2 H1038R']
                )

                if match_found:
                    data_row['WHO-HAEM5'] = class_row['WHO-HAEM5']
                    data_row['ICC'] = class_row['ICC']
                    matched_rows.append(data_row)

    if matched_rows:
        return pd.DataFrame(matched_rows)
    else:
        return None


def main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir, classification_file,
         output_csv, output_text, output_curation, output_driver):
    # Lade die Klassifikationsdatei
    df = pd.read_csv(classification_file, sep=',')

    # Extrahiere die Daten
    data_df = gather_data(
        allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir, classification_file)

    # Führe die Überprüfung der Bedingungen durch
    results = check_conditions(data_df, df)
    output_df = pd.DataFrame(results)
    data_df['Unique_spanning_reads'] = pd.to_numeric(data_df['Unique_spanning_reads'], errors='coerce')
    filtered_data_df = data_df[data_df['Unique_spanning_reads'] > 2]
    helper_df = []
    if (len(filtered_data_df) and len(output_df)) == int(0):
        matched_rows = []
        helper_df = data_df.iloc[0]
        helper_df = helper_df[["ALLCatchR", "Confidence", "karyotype_classifier", "PAX5_P80R", "IKZF1_N159Y", "ZEB2_H1038R"]]
        for _, class_row in df.iterrows():
            match_found = (
                    helper_df['ALLCatchR'] == class_row['ALLCatchR'] and
                    helper_df['Confidence'] == class_row['Confidence'] and
                    helper_df['karyotype_classifier'] == class_row['karyotype classifier'] and
                    helper_df['PAX5_P80R'] == class_row['PAX5 P80R'] and
                    helper_df['IKZF1_N159Y'] == class_row['IKZF1 N159Y'] and
                    helper_df['ZEB2_H1038R'] == class_row['ZEB2 H1038R']
            )

            if match_found:
                helper_df['WHO-HAEM5'] = class_row['WHO-HAEM5']
                helper_df['ICC'] = class_row['ICC']
                helper_df['Gene_1_symbol(5end_fusion_partner)'] = None
                helper_df['Gene_2_symbol(3end_fusion_partner)'] = None
                helper_df['Fusioncaller'] = None
                helper_df['Unique_spanning_reads'] = None
                matched_rows.append(helper_df)
                output_df = pd.DataFrame(matched_rows)
    else:
        print(" ")

    # Prüfe, ob `output_df` leer ist
    if output_df.empty:
        # Ergänze die `WHO-HAEM5` und `ICC` Spalten mit Standardwerten
        data_df['WHO-HAEM5'] = "IntegrateALL couldn't confirm the subtype in concordance with WHO or ICC classification"
        data_df['ICC'] = ""  # Leerer Eintrag für `ICC`
        # Speichere die Daten und bearbeite die Fusions- und Textdateien
        data_df.to_csv(output_csv, index=False)
        handle_driver_fusions(data_df, output_driver)
        write_no_match(output_text, output_curation, data_df)
    elif (not output_df.empty) and (len(output_df) != len(filtered_data_df)):

        if len(filtered_data_df) > len(output_df):
            filtered_data_df[
                'WHO-HAEM5'] = "IntegrateALL couldn't confirm the subtype in concordance with WHO or ICC classification"
            filtered_data_df['ICC'] = ""  # Leerer Eintrag für `ICC`
            # Speichere die Daten und bearbeite die Fusions- und Textdateien
            filtered_data_df.to_csv(output_csv, index=False)
            handle_driver_fusions(filtered_data_df, output_driver)
            write_no_match(output_text, output_curation, filtered_data_df)
        else:
            if len(output_df) > 1:
                # Wenn mehrere Einträge vorhanden sind, handle die Mehrfacheinträge
                handle_multiple_entries(output_df, output_text, output_curation)
                output_df.to_csv(output_csv, index=False)
                handle_driver_fusions(data_df, output_driver)
            else:
                # Wenn nur ein Eintrag vorhanden ist, verarbeite den Einzelnen
                write_single_entry_text(output_df.iloc[0], output_text, output_curation)
                output_df.to_csv(output_csv, index=False)
                handle_driver_fusions(data_df, output_driver)


    else:
        # Wenn `output_df` die gleiche Länge wie `data_df` hat, verarbeite es weiter
        if len(output_df) > 1:
            # Wenn mehrere Einträge vorhanden sind, handle die Mehrfacheinträge
            handle_multiple_entries(output_df, output_text, output_curation)
            output_df.to_csv(output_csv, index=False)
            handle_driver_fusions(data_df, output_driver)
        else:
            # Wenn nur ein Eintrag vorhanden ist, verarbeite den Einzelnen
            write_single_entry_text(output_df.iloc[0], output_text, output_curation)
            output_df.to_csv(output_csv, index=False)
            handle_driver_fusions(data_df, output_driver)


def write_no_match(output_text, output_curation, output_entry_text):
    fusion_details = ""
    if output_entry_text['Fusioncaller'].isnull().all():
        # Platzhalterwerte für den Fall, dass keine Fusionen vorhanden sind
        fusion_details = "No driver fusions identified"
    else:
        # Wenn Fusionen vorhanden sind, sammle die Details
        fusion_details = ""
        for _, row in output_entry_text.iterrows():
            fusion_details += (
                f"{row['Fusioncaller']}: {row['Gene_1_symbol(5end_fusion_partner)']}::{row['Gene_2_symbol(3end_fusion_partner)']}; "
            )
        fusion_details = fusion_details.rstrip("; ")

    output_entry_text = output_entry_text.iloc[0]

    with open(output_text, 'w') as f:
        f.write(f"IntegrateALL classification:\n\n"
                f"Gene expression-based subtype-allocation (ALLCatchR: {output_entry_text['ALLCatchR']} subtype, Confidence: {output_entry_text['Confidence']}) "
                f"and the genomic driver profile ({fusion_details}, "
                f"RNA-Seq CNV karyotype classifier: {output_entry_text['karyotype_classifier']}, subtype defining SNPs: "
                f"{'PAX5 P80R present' if output_entry_text['PAX5_P80R'] else 'PAX5 P80R absent'}, "
                f"{'IKZF1 N159Y present' if output_entry_text['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}) "
                f"seem not to be consistent with an unambiguous diagnostic classification according to WHO-HAEM5 "
                f"(Alaggio R et al. Leukemia, 2022) / ICC (Arber D et al. Blood, 2022).")

    with open(output_curation, 'w') as c:
        c.write(
            "Classification, Subtype,Confidence,Fusion_details,Karyotype_classifier,PAX5_P80R,IKZF1_N159Y,ZEB2_H1038R,WHO-HAEM5,ICC\n")
        c.write(f"{'Manual Curation'}, "
            f"{output_entry_text['ALLCatchR']}, {output_entry_text['Confidence']}, {fusion_details}, {output_entry_text['karyotype_classifier']}, "
            f"{'PAX5 P80R present' if output_entry_text['PAX5_P80R'] else 'PAX5 P80R absent'}, "
            f"{'IKZF1 N159Y present' if output_entry_text['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}, "
            f"{'ZEB2 H1038R present' if output_entry_text['ZEB2_H1038R'] else 'ZEB2 H1038R absent'}, "
            f"{''}, {''}\n"
            )


def handle_multiple_entries(output_df, output_text, output_curation):
    if output_df['ICC'].nunique() == 1:
        write_multiple_entry_text(output_df, output_text, output_curation)

    else:
        output_df['Unique_spanning_reads'] = pd.to_numeric(output_df['Unique_spanning_reads'], errors='coerce')
        filtered_output_df = output_df[output_df['Unique_spanning_reads'] > 2]

        if filtered_output_df.empty:
            write_no_match(output_text, output_df, output_curation)

        elif filtered_output_df['ICC'].nunique() != 1:
            write_multiple_entry_with_different_ICC_text(filtered_output_df, output_text, output_curation)
        else:
            selected_entry = filtered_output_df.iloc[0].copy()
            write_single_entry_text(selected_entry, output_text, output_curation)



def write_single_entry_text(entry, output_text, output_curation):
    with open(output_text, 'w') as f:
        f.write(
            f"Consistency between gene expression-based subtype-allocation (ALLCatchR: {entry['ALLCatchR']} subtype, Confidence: {entry['Confidence']}) "
            f"and the genomic driver profile ({entry['Fusioncaller']}: {entry['Gene_1_symbol(5end_fusion_partner)']}::{entry['Gene_2_symbol(3end_fusion_partner)']}, "
            f"RNA-Seq CNV karyotype classifier: {entry['karyotype_classifier']}, subtype defining SNPs: "
            f"{'PAX5 P80R present' if entry['PAX5_P80R'] else 'PAX5 P80R absent'}, "
            f"{'IKZF1 N159Y present' if entry['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}) "
            f"supports a classification as\n\n"
            f"{entry['WHO-HAEM5']} according to WHO-HAEM5 (Alaggio R et al. Leukemia, 2022) and\n"
            f"{entry['ICC']} according to ICC (Arber D et al. Blood, 2022).\n\n")

    with open(output_curation, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Schreibe den Header
        writer.writerow([
            "Classification", "Subtype", "Confidence", "Fusion_details",
            "Karyotype_classifier", "PAX5_P80R", "IKZF1_N159Y",
            "ZEB2_H1038R", "WHO-HAEM5", "ICC"
        ])

        # Schreibe die Daten
        writer.writerow([
            "Automatic Classification",
            entry['ALLCatchR'],
            entry['Confidence'],
            f"{entry['Fusioncaller']}: {entry['Gene_1_symbol(5end_fusion_partner)']}::{entry['Gene_2_symbol(3end_fusion_partner)']}",
            entry['karyotype_classifier'],
            "PAX5 P80R present" if entry['PAX5_P80R'] else "PAX5 P80R absent",
            "IKZF1 N159Y present" if entry['IKZF1_N159Y'] else "IKZF1 N159Y absent",
            "ZEB2 H1038R present" if entry['ZEB2_H1038R'] else "ZEB2 H1038R absent",
            entry['WHO-HAEM5'],
            entry['ICC']
        ])




def write_multiple_entry_with_different_ICC_text(output_df, output_text, output_curation):
    fusion_details = ""
    for _, row in output_df.iterrows():
        fusion_details += (
            f"{row['Fusioncaller']}: {row['Gene_1_symbol(5end_fusion_partner)']}::{row['Gene_2_symbol(3end_fusion_partner)']}; "
        )
    fusion_details = fusion_details.rstrip("; ")

    output_entry_text = output_df.iloc[0]

    with open(output_text, 'w') as f:
        f.write(f"IntegrateALL classification:\n\n"
                f"Gene expression-based subtype-allocation (ALLCatchR: {output_entry_text['ALLCatchR']} subtype, Confidence: {output_entry_text['Confidence']}) "
                f"and the genomic driver profile ({fusion_details}, "
                f"RNA-Seq CNV karyotype classifier: {output_entry_text['karyotype_classifier']}, subtype defining SNPs: "
                f"{'PAX5 P80R present' if output_entry_text['PAX5_P80R'] else 'PAX5 P80R absent'}, "
                f"{'IKZF1 N159Y present' if output_entry_text['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}) "
                f"seem not to be consistent with an unambiguous diagnostic classification according to WHO-HAEM5 "
                f"(Alaggio R et al. Leukemia, 2022) / ICC (Arber D et al. Blood, 2022).")

    with open(output_curation, 'w') as c:
        c.write(
            "Classification, Subtype,Confidence,Fusion_details,Karyotype_classifier,PAX5_P80R,IKZF1_N159Y,ZEB2_H1038R,WHO-HAEM5,ICC\n")
        c.write(f"{'Manual Curation'}, "
            f"{output_entry_text['ALLCatchR']}, {output_entry_text['Confidence']}, {fusion_details}, {output_entry_text['karyotype_classifier']}, "
            f"{'PAX5 P80R present' if output_entry_text['PAX5_P80R'] else 'PAX5 P80R absent'}, "
            f"{'IKZF1 N159Y present' if output_entry_text['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}, "
            f"{'ZEB2 H1038R present' if output_entry_text['ZEB2_H1038R'] else 'ZEB2 H1038R absent'}, "
             f"{''}, {''}\n"
            )



def write_multiple_entry_text(output_df, output_text, output_curation):
    if output_df['ICC'].nunique() == 1:
        selected_entry = output_df.iloc[0].copy()
        entry = selected_entry

        # Überprüfen, ob es verschiedene Fusioncaller gibt und ergänze sie
        additional_fusioncallers = output_df['Fusioncaller'].unique()
        if len(additional_fusioncallers) > 1:
            # Kombiniere alle unterschiedlichen Fusioncaller zu einem String
            combined_fusioncallers = " & ".join(additional_fusioncallers)
            entry.loc['Fusioncaller'] = combined_fusioncallers

    with open(output_text, 'w') as f:
        f.write(
            f"Consistency between gene expression-based subtype-allocation (ALLCatchR: {entry['ALLCatchR']} subtype, Confidence: {entry['Confidence']}) "
            f"and the genomic driver profile ({entry['Fusioncaller']}: {entry['Gene_1_symbol(5end_fusion_partner)']}::{entry['Gene_2_symbol(3end_fusion_partner)']}, "
            f"RNA-Seq CNV karyotype classifier: {entry['karyotype_classifier']}, subtype defining SNPs: "
            f"{'PAX5 P80R present' if entry['PAX5_P80R'] else 'PAX5 P80R absent'}, "
            f"{'IKZF1 N159Y present' if entry['IKZF1_N159Y'] else 'IKZF1 N159Y absent'}) "
            f"supports a classification as\n\n"
            f"{entry['WHO-HAEM5']} according to WHO-HAEM5 (Alaggio R et al. Leukemia, 2022) and\n"
            f"{entry['ICC']} according to ICC (Arber D et al. Blood, 2022).\n\n")

    with open(output_curation, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Schreibe den Header
        writer.writerow([
            "Classification", "Subtype", "Confidence", "Fusion_details",
            "Karyotype_classifier", "PAX5_P80R", "IKZF1_N159Y",
            "ZEB2_H1038R", "WHO-HAEM5", "ICC"
        ])

        # Schreibe die Daten
        writer.writerow([
            "Automatic Classification",
            entry['ALLCatchR'],
            entry['Confidence'],
            f"{entry['Fusioncaller']}: {entry['Gene_1_symbol(5end_fusion_partner)']}::{entry['Gene_2_symbol(3end_fusion_partner)']}",
            entry['karyotype_classifier'],
            "PAX5 P80R present" if entry['PAX5_P80R'] else "PAX5 P80R absent",
            "IKZF1 N159Y present" if entry['IKZF1_N159Y'] else "IKZF1 N159Y absent",
            "ZEB2 H1038R present" if entry['ZEB2_H1038R'] else "ZEB2 H1038R absent",
            entry['WHO-HAEM5'],
            entry['ICC']
        ])



def handle_driver_fusions(data_df, output_driver):
    data_df['Unique_spanning_reads'] = pd.to_numeric(data_df['Unique_spanning_reads'], errors='coerce')
    if 'DUX4' in data_df['ALLCatchR'].values:
        # Wenn DUX4 subtype vorhanden ist, filtere ohne die DUX4-Einträge zu entfernen
        filtered_df = data_df[
            ~((data_df['Unique_spanning_reads'].isna()) |  # Schließt NaN-Werte aus
              (data_df['Unique_spanning_reads'] < 3))  # Schließt Werte < 3 aus
        ]
    else:
        # Normale Filterung, DUX4-Einträge werden ausgeschlossen
        filtered_df = data_df[
            ~((data_df['Gene_1_symbol(5end_fusion_partner)'] == 'DUX4') |
              (data_df['Gene_2_symbol(3end_fusion_partner)'] == 'DUX4') |
              (data_df['Unique_spanning_reads'].isna()) |  # Schließt NaN-Werte aus
              (data_df['Unique_spanning_reads'] < 3))  # Schließt Werte < 3 aus
        ]
    if filtered_df.empty:
        with open(output_driver, 'w') as empty_file:
            pass  # Erstelle eine leere Datei ohne Inhalt

    else:
        filtered_df = filtered_df[[
            'Gene_1_symbol(5end_fusion_partner)',
            'Gene_2_symbol(3end_fusion_partner)',
            'Unique_spanning_reads',
            'Fusioncaller'
        ]]

        # Schreibe die Daten in eine CSV-Datei
        filtered_df.to_csv(output_driver, index=False, header=True)



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
    output_text = sys.argv[9]
    output_curation = sys.argv[10]
    output_driver = sys.argv[11]
    main(sample, allcatchr_file, karyotype_file, fusioncatcher_file, arriba_file, hotspot_dir, classification_file,
         output_csv, output_text, output_curation, output_driver)
