import os
import re
from datetime import datetime

import sys
import glob
import pandas as pd
from bs4 import BeautifulSoup

# Helper function to find the next sibling text after a specific label
def extract_from_p_tags(soup, startswith="I am blue"):
    results = []

    for p in soup.find_all('p'):
        text = p.get_text(separator=' ').strip()
        if text.startswith(startswith):
            parts = text.split(':', 1)
            if len(parts) > 1:
                results.append(parts[1].strip(".").strip())

    if len(results) > 0:
        return results[0]
    return None


def extract_qc_settings(html_file):

    cohort_name_from_file = re.search("Report_DataQc_(.+).html$", html_file).group(1)

    last_modified = datetime.fromtimestamp(os.path.getmtime(html_file)).isoformat()

    df_pandas = pd.read_html(html_file, flavor='bs4', thousands ='.')

    # Load the HTML file
    with open(html_file, 'r', encoding='utf-8') as file:
        soup = BeautifulSoup(file, 'html.parser')

    # Extract the required fields
    cohort = extract_from_p_tags(soup, "This report is for cohort:")
    expression_platform = extract_from_p_tags(soup, "Expression platform for this cohort is:")
    post_qc_sample_size = extract_from_p_tags(soup, "Post-QC sample size for this cohort is:")
    lof_s_threshold = extract_from_p_tags(soup, "LOF S-threhold for genetic outliers:")
    threshold_sd_genetic = extract_from_p_tags(soup, "Threshold of SD difference from mean genetic PC1 and PC2 is:")
    threshold_sd_expression = extract_from_p_tags(soup, "Threshold of SD difference from mean gene expression PC1 and PC2 is:")
    angle_contamination = extract_from_p_tags(soup, "Angle of contamination area in sex chromosome gene expression plot:")
    slope_sex_chromosome = extract_from_p_tags(soup, "Slope in sex chromosome gene expression plot:")
    # Create a DataFrame
    data = {
        "File": [html_file],
        "Last Modified": [last_modified],
        "Cohort": [cohort_name_from_file],
        "Expression Platform": [expression_platform],
        "Post-QC Sample Size": [post_qc_sample_size],
        "LOF S-threshold": [lof_s_threshold],
        "Threshold SD Genetic": [threshold_sd_genetic],
        "Threshold SD Expression": [threshold_sd_expression],
        "Angle Contamination": [angle_contamination],
        "Slope Sex Chromosome": [slope_sex_chromosome],
        "Nr of samples before QC": [df_pandas[0].loc[2, "Nr. of eQTL samples"]],
        "Nr of samples after QC": [df_pandas[1].loc[df_pandas[1].shape[0]-1, "Nr. of samples"]],
        "Nr of variants before QC": [df_pandas[0].loc[0, "Nr. of SNPs"]],
        "Nr of variants after QC": [df_pandas[0].loc[df_pandas[1].shape[0]-1, "Nr. of SNPs"]],
        "Nr of genes before QC": [df_pandas[1].loc[0, "Nr. of features"]],
        "Nr of genes after QC": [df_pandas[1].loc[df_pandas[1].shape[0]-1, "Nr. of features"]]
    }
    return pd.DataFrame(data)

renaming_dict = {
    "cruchagalab_MATT_EUR": "Cruchagalab_EUR",
    "cruchagalab_MATT_AFR": "Cruchagalab_AFR",
    "YFS": "YoungFinns",
    "LIFEHEART": "LIFE_HEART",
    "LIFEADULT": "LIFE_ADULT"
}

files = glob.glob("manuscript_freeze_2025-06-20/Report_DataQc_*.html")

concatenated_df = pd.concat([extract_qc_settings(html_file) for html_file in files])
concatenated_df.Cohort = concatenated_df.Cohort.replace(renaming_dict)

# Save the DataFrame to a CSV file
concatenated_df.to_csv('cohort_qc_info.tsv', sep="\t", index=False)

print("The cohort QC information has been extracted and saved to cohort_qc_info.csv.")
