#!/usr/bin/env python
"""
Author: Evy Kuenen
Date: 26-2-2025
Functionality: This script makes a training and testing file for the machine learning model.

Usage:
    Required directory structure:
                    make_train_and_test_file.py needs to be in scripts/make_train_and_test_file directory
                     from the ..\\..\\output from scripts directory: ani scores will be retrieved from
                    ani_score_output\\ani_scores.csv
                    fragmentation will be retrieved from fragmentation_output\\fragmentation.csv
                    host virus correlation will be retrieved from
                    correlation_analysis_output\\host_virus_correlation_matrix.csv
                    labels will be retrieved from sample_ids_contaminated_samples\\samples_contaminated.csv
                    count and length will be retrieved from count_contigs\\count_contigs.csv
                    and count_contigs\\length_contigs.csv
    Required files: ani_scores.csv, fragmentation.csv, host_virus_correlation_matrix.csv, samples_contaminated.csv,
                    count_contigs.csv and length_contigs.csv
    Calling script: "python 3 make_train_and_test_file.py"
"""
import numpy as np
import pandas as pd
import sys

def read_fragmentation(path_fragmentation):
    try:
        with open(path_fragmentation, "r", encoding="utf-8") as file:
            for line in file:
                splitted_line_fragmentation = line.replace('\n', '').split(',')
                if not '' in splitted_line_fragmentation:
                    return splitted_line_fragmentation
    except:
        print('cant open file!')


def read_host_virus_correlation(path_host_virus_correlation, host, virus):
    try:
        df = pd.read_csv(path_host_virus_correlation, index_col=0, engine='python')
        correlation = df.loc[virus, host]
    except Exception as e:
        print('host not found or was unknown in file.', e)
        correlation = 0
    return correlation


def make_contamination_set(path_contamination):
    contaminated_sample_ids = set()
    try:
        with open(path_contamination, "r", encoding="utf-8") as file:
            for line in file:
                line.strip()
                line_without_space = line.replace('\n', '').replace('[', '').replace(']', '').replace('\'', '')
                contaminated_sample_ids.add(line_without_space)
        return contaminated_sample_ids
    except:
        print('cant open file!')


def load_dict_from_file(path):
    data_dict = {}
    with open(path, 'r', encoding='utf-8') as f:
         for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 4:
                sample, _, virus, value = parts[:4]
                virus = virus.replace('[\'', '').replace('\']', '')
                data_dict[(sample, virus)] = value
    return data_dict


def get_train_line(path_ani_scores, path_host_virus_correlation, path_fragmentation,
                   path_contig_count, path_length_contigs):
    train_data = []

    contig_counts_dict = load_dict_from_file(path_contig_count)
    contig_lengths_dict = load_dict_from_file(path_length_contigs)
    fragmentation_dict = load_dict_from_file(path_fragmentation)
   
    with open(path_ani_scores, "r", encoding="utf-8") as file:
        for line in file:
            parts = line.strip().split(',')
            if '' in parts or '"' in parts:
                continue

            sample, host, virus = parts[:3]
            virus = virus.replace('[\'', '').replace('\']', '').strip()
            ani = round(float(parts[3]), 2)
            correlation = round(read_host_virus_correlation(path_host_virus_correlation, host, virus), 3)
            sample = sample.strip()
            virus = virus.strip()
            print(str(sample))
            print(virus)
            contig_count = int(contig_counts_dict.get((sample, virus), 1))
            contig_length = int(contig_lengths_dict.get((sample, virus), 1))
            fragmentation = round(float(fragmentation_dict.get((sample, virus), 1)), 2)

           
            if contig_count > 0 and contig_length > 0:
               average_contig_length = round(float(contig_length / contig_count))
               train_line = [sample, host, virus, ani, correlation, fragmentation,
                                  average_contig_length, contig_count]
               print(train_line)
               train_data.append(tuple(train_line))

    train_data = list(set(train_data))
    return train_data



def split_test_train_data(train_data, prediction_file):
    """
    Splits data into train and test sets after augmentation.

    Args:
        train_data (list): List of original training data.
        train_file (str): Path to save train data.
        test_file (str): Path to save test data.
    """
    df = pd.DataFrame(train_data,
                      columns=["sample", "host", "virus", "ani", "correlation", "fragmentation",
                               "average contig length", "contig_count"])

    df.to_csv(prediction_file, index_label=False)



def main():
    # input
    path_ani_scores = sys.argv[1]
    path_fragmentation = sys.argv[2]
    path_host_virus_correlation = sys.argv[3]
    path_contig_count = sys.argv[4]
    path_length_contigs = sys.argv[5]
    # output
    prediction_file = sys.argv[6]


    read_fragmentation(path_fragmentation)
    train_data = get_train_line(path_ani_scores, path_host_virus_correlation, path_fragmentation,
                                 path_contig_count, path_length_contigs)
    split_test_train_data(train_data, prediction_file)



main()