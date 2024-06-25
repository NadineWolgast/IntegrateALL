import csv
import pickle
import sys
from collections import Counter

import numpy as np
import pandas as pd


in_file = sys.argv[1]
outfile = sys.argv[2]
sample_id = sys.argv[3]

np.set_printoptions(suppress=True)

def generate_dummy_output(outfile, sample_id):
    # Erstelle den Dummy-Output
    with open(outfile, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Schreibe die Header
        writer.writerow(['Sample', 'Score', 'Prediction', 'Hyperdiploid', 'Low hypodiploid', 'Near haploid', 'iAMP21', 'Other'])
        writer.writerow([sample_id, '0%', 'unclassified', '0%', '0%', '0%', '0%', '0%'])

def generate_prediction_statement(predictions):
    # Zähle die Häufigkeit der einzelnen Vorhersagen
    prediction_counts = Counter(predictions)
    print("prediction_counts", prediction_counts)
    # Finde den am häufigsten vorhergesagten Subtyp
    most_common_subtype = prediction_counts.most_common(1)[0][0]
    print("most_common_subtype", most_common_subtype)
    print(prediction_counts[most_common_subtype])
    # Berechne den Prozentsatz der Vorhersagen für den häufigsten Subtyp
    confidence_percentage = (prediction_counts[most_common_subtype] / len(predictions)) * 100
    # Erstelle die Aussage
    prediction_statement = f"Predicted subtype is {most_common_subtype} with confidence {confidence_percentage:.2f}%"
    return prediction_statement


with open('scripts/ensemble_classifier_250524.pkl', 'rb') as model_file:
    model_data = pickle.load(model_file)
    rf_classifier = model_data['model']
    scaler = model_data['scaler']
    label_encoder = model_data['label_encoder']

new = pd.read_csv(in_file)

# Prüfen der Spaltenanzahl
if new.shape[1] < 185:
    print(f"Die Eingabedatei {in_file} hat weniger als 185 Spalten. Generiere Dummy-Output.")
    generate_dummy_output(outfile, sample_id)
else:
    new_data = new.iloc[:, 1:]

    # Skalierung der Daten mit dem geladenen Scaler
    scaled_data = scaler.transform(new_data)

    pred = rf_classifier.predict(scaled_data)
    pred_pro = rf_classifier.predict_proba(scaled_data)

    # Extrahiere das vorhergesagte Label
    predicted_label = label_encoder.inverse_transform([pred[0]])[0]
    index_predicted_label = list(rf_classifier.classes_).index(pred[0])
    probability = pred_pro[0][index_predicted_label]

    # Convert probabilities to whole numbers as percentages
    print(probability, pred_pro[0])
    pred_pro_percent = [f"{int(prob * 100)}%" for prob in pred_pro[0]]

    with open(outfile, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Schreibe die Header
        writer.writerow(['Sample', 'Score', 'Prediction', 'Hyperdiploid', 'Low hypodiploid', 'Near haploid', 'iAMP21', 'Other'])
        writer.writerow([sample_id, f"{int(probability * 100)}%", predicted_label] + pred_pro_percent)
