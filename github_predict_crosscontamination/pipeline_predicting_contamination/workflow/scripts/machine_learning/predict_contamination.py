"""
Author: Evy Kuenen
Date: 15-5-2025
Functionality: This script predicts the labels of the inputted prediction file without labels.
Usage:
    Required directory structure:
                    predict_contamination.py needs to be in scripts/machine_learning
                    prediction file will be retrieved from
                    "..\\..\\output_from_scripts\\test_train_data_ml\\prediction.csv"
    Required files: prediction.csv
    Calling script: "python 3 predict_contamination.py"
"""
import sys
import base64
import os
import webbrowser
from io import BytesIO
import joblib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from pathlib import Path


def process_before_predict(df):
    """
    Prepares the dataframe for prediction by selecting specific columns.

    Parameters:
    df (pandas.DataFrame): The input dataframe containing the data to
        be processed.

    Returns:
    pandas.DataFrame: A new dataframe with only the selected columns.
    """
    used_columns = ['ani', "correlation", "fragmentation", "average contig length"]

    # Get only needed columns
    df_new = df[used_columns].copy()
    print(df["average contig length"])
    print(df.describe())
    return df_new


def predict(df, model_path):
    """
    Predicts the pathotype of samples in the dataframe using a
    pre-trained model.

    Parameters:
    df (pandas.DataFrame): The input dataframe containing the samples to
        be predicted.
    model_path (str): The file path to the pre-trained model.
    unique_variants_path (str): The file path to the CSV file containing
        unique variants.

    Returns:
    str: An HTML string representing the results table with predicted
        frequencies and percentages.
    """
    df_2 = process_before_predict(df)
    model = joblib.load(model_path)

    proba = model.predict_proba(df_2.values)
    predictions = model.predict(df_2.values)

    result_df = pd.DataFrame({
        "sample": df["sample"].values,
        "virus": df["virus"].values,
        "host": df["host"].values,
        "prediction": predictions,
        "prob_class_0": proba[:, 0],
        "prob_class_1": proba[:, 1],
        "ani": df["ani"].values,
        "correlation": df["correlation"].values,
        "fragmentation": df["fragmentation"].values,
        "average contig length": df["average contig length"].values,
        "count contigs": df["count_contigs"].values,
    })

    print(result_df)
    return result_df


def visualize_prediction_html(prediction_df, template_path, output_path):
    # Plot maken
    fig, ax = plt.subplots(figsize=(7, 4))
    samples = prediction_df["sample"]
    viruses = prediction_df["virus"]
    sample_virus_labels = [f"{s}\n{v}" for s, v in zip(samples, viruses)]
    probs = prediction_df["prob_class_1"].replace(0, 1e-6)
    predictions = prediction_df["prediction"]
    colors = ["#e14b31" if pred == 1 else "#a7d5ed" for pred in predictions]

    ax.bar(sample_virus_labels, probs, color=colors)
    ax.set_title("Probability of contamination per sample-virus pair")
    ax.set_xlabel("Sample\nVirus")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
    ax.set_ylabel("Probability class 1 (contaminated)")
    plt.setp(ax.get_xticklabels(), rotation=20, fontsize=8)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.setp(ax.get_xticklabels(), rotation=20)

    legend_elements = [
        Patch(facecolor="#e14b31", label='Contaminated (1)'),
        Patch(facecolor="#a7d5ed", label='Clean (0)')
    ]
    ax.legend(handles=legend_elements, title="Classification", loc='lower left', bbox_to_anchor=(1, 0.5))

    buf = BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format="png")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    img_html = f'<img src="data:image/png;base64,{encoded}" style="max-width:100%;">'

    # Tabel
    column_renames = {
        "prob_class_0": "P(clean) ",
        "prob_class_1": "P(contaminated) ",
        "correlation": "Host-virus correlation ",
        "sample": "Sample ",
        "virus": "Virus ",
        "host": "Host ",
        'prediction': "Prediction ",
        'ani': "ANI ",
        "fragmentation": "Fragmentation ",
        "average contig length": "Average contig length ",
        "count_contigs": "Count contigs ",
    }
    table_html = "<table border='1'><thead><tr>" + \
                 "".join([f"<th>{column_renames.get(col, col)}</th>" for col in prediction_df.columns]) + \
                 "</tr></thead><tbody>"

    for _, row in prediction_df.iterrows():
        table_html += "<tr>"
        is_contaminated = row.get("prediction") == 1
        row_style = " style='background-color: #de6e56; color: white;'" if is_contaminated \
            else " style='background-color: #63bff0; color: white;'"
        table_html += f"<tr{row_style}>"
        for col, val in row.items():
            if col == "prediction":
                val = "gecontamineerd" if val == 1 else "clean"
            elif isinstance(val, float):
                val = f"{val:.3f}"
            table_html += f"<td>{val}</td>"
        table_html += "</tr>"

    table_html += "</tbody></table>"

    # read template and fill in with information from script
    with open(template_path, "r") as f:
        template = f.read()

    html_filled = template.format(img_html=img_html, table_html=table_html)

    with open(output_path, "w") as f:
        f.write(html_filled)

    print(f"HTML popup gegenereerd als '{output_path}'")
    file_path = os.path.abspath(output_path)
    webbrowser.open_new_tab(f"file://{file_path}")


def main():
    # input
    df_path = sys.argv[1]
    model_path = sys.argv[2]
    
    # output
    template_path = Path(sys.argv[3])
    output_path = sys.argv[4]

    df = pd.read_csv(df_path)
    prediction_df = predict(df, model_path)
    visualize_prediction_html(prediction_df, template_path, output_path)

main()