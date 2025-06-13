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
    Calling script: "python make_train_and_test_file.py"
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split


def read_fragmentation(path_fragmentation):
    """
    read fragmentation from file with fragmentation and split line into sample, virus and fragmentation score.
    :param path_fragmentation: path to file with fragmentation
    :return splited_line_fragmentation: list with sample, virus, fragmentation score
    """
    try:
        with open(path_fragmentation, "r", encoding="utf-8") as file:
            for line in file:
                splitted_line_fragmentation = line.replace('\n', '').split(',')
                if not '' in splitted_line_fragmentation:
                    return splitted_line_fragmentation
    except:
        print('cant open file!')


def read_host_virus_correlation(path_host_virus_correlation, host, virus):
    """
    Read host virus correlation from correlation matrix
    :param path_host_virus_correlation: path to matrix with host virus correlation
    :param host: name of host
    :param virus: name of virus
    :return: correlation score
    """
    try:
        host = host.strip().lower()
        virus = virus.strip().lower()
        df = pd.read_csv(path_host_virus_correlation, index_col=0, engine='python')
        df.columns = df.columns.str.strip().str.lower()
        df.index = df.index.str.strip().str.lower()
        correlation = df.loc[virus, host]
    except:
        print('host not found or was unknown in file.')
        correlation = 0
    return correlation


def make_contamination_set(path_contamination):
    """
    get label if sample and virus in sample is contamination and put in list.
    :param path_contamination: path to file with contaminated samples and viruses
    :return: set with ids of contaminated samples
    """
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
    """
    make dictionary from files with sample and virus as key and score as value.
    :param path: path that is given to function
    :return data dict: dictionary with as key sample and virus and as value value from file.
    """
    data_dict = {}
    with open(path, 'r', encoding='utf-8') as f:
         for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 4:
                sample, _, virus, value = parts[:4]
                virus = virus.replace('[\'', '').replace('\']', '')
                data_dict[(sample, virus)] = value
    return data_dict


def get_train_line(path_ani_scores, path_host_virus_correlation, path_fragmentation, contaminated_sample_ids,
                   path_contig_count, path_length_contigs):
    """
    make dataset for train test and validation
    :param path_ani_scores: path to ani scores
    :param path_host_virus_correlation: path to host virus correlation matrix
    :param path_fragmentation: path to file with fragmentation scores
    :param contaminated_sample_ids: set with contaminated samples
    :param path_contig_count: path to file with contig counts
    :param path_length_contigs: path to file with contig lengths
    :return train_data: dataframe with sample, virus, ani, corr, frag, average length, count, contamination
    """
    train_data = []

    contig_counts_dict = load_dict_from_file(path_contig_count)
    contig_lengths_dict = load_dict_from_file(path_length_contigs)
    fragmentation_dict = load_dict_from_file(path_fragmentation)

    contaminated_dict = {}
    contaminated_samples = []
    for item in contaminated_sample_ids:
        sample = item.replace("\',)", "").replace("(\'", "").replace(")", "").replace("(", "").split(',')[0]
        contaminated_samples.append(sample)

    for line in contaminated_sample_ids:
        line = line.replace("\',)", "").replace("(\'", "").replace(")", "").replace("(", "")
        if ',' in line:
            sample, virus = map(str.strip, line.split(','))
            if len(virus) > 0:
                sample = sample.strip()
                virus = virus.strip()
                contaminated_dict[sample] = virus
        else:
            contaminated_dict[line] = None

    with open(path_ani_scores, "r", encoding="utf-8") as file:
        for line in file:
            parts = line.strip().split(',')
            if '' in parts or '"' in parts:
                continue

            sample, host, virus = parts[:3]
            virus = virus.replace('[\'', '').replace('\']', '').strip()
            sample = sample.strip()
            virus = virus.strip()
            host = host.strip()
            ani = round(float(parts[3]), 4)
            correlation = round(read_host_virus_correlation(path_host_virus_correlation, host, virus), 4)
            print(str(sample))
            print(virus)
            print(host)
            contig_count = int(contig_counts_dict.get((sample, virus), 1))
            contig_length = int(contig_lengths_dict.get((sample, virus), 1))
            fragmentation = round(float(fragmentation_dict.get((sample, virus), 1)), 4)

            # Labels toekennen
            if sample in contaminated_dict:
                expected_virus = contaminated_dict[sample]
                label = 1 if expected_virus is None or virus == expected_virus else 0

                if contig_count > 0 and contig_length > 0:
                    average_contig_length = round(float(contig_length / contig_count))
                    train_line = [sample, host, virus, ani, correlation, fragmentation,
                                  average_contig_length, contig_count, label]
                    print(train_line)
                    train_data.append(tuple(train_line))
            if '000000' in sample:
                label = 1
                if contig_count > 0 and contig_length > 0:
                    average_contig_length = round(float(contig_length / contig_count))
                    train_line = [sample, host, virus, ani, correlation, fragmentation,
                                  average_contig_length, contig_count, label]
                    print(train_line)
                    train_data.append(tuple(train_line))
            if sample in str(contaminated_samples) and sample not in contaminated_dict and '000000' not in sample:
                label = 1
                if contig_count > 0 and contig_length > 0:
                    average_contig_length = round(float(contig_length / contig_count))
                    train_line = [sample, host, virus, ani, correlation, fragmentation,
                                  average_contig_length, contig_count, label]
                    print(train_line)
                    train_data.append(tuple(train_line))

            if sample not in str(contaminated_samples)and '000000' not in sample:
                label = 0
                if contig_count > 0 and contig_length > 0:
                    average_contig_length = round(float(contig_length / contig_count))
                    train_line = [sample, host, virus, ani, correlation, fragmentation,
                                  average_contig_length, contig_count, label]
                    print(train_line)
                    train_data.append(tuple(train_line))

    train_data = list(set(train_data))
    return train_data


def get_feature_correlation(train_data):
    """
    get feature correlation if features are correlated.
    :param train_data: dataset for train, test and validation
    :return correlation matrix: matrix of how much features are correlated
    """
    df = pd.DataFrame(train_data, columns=[
        'sample', 'host', 'virus', 'ani', 'correlation',
        'fragmentation', 'average_contig_length', 'contig_count', 'label'
    ])

    numeric_df = df[['ani', 'correlation', 'fragmentation', 'average_contig_length', 'contig_count']]
    correlation_matrix = numeric_df.corr()

    # Plotting with matplotlib
    fig, ax = plt.subplots(figsize=(8, 6))
    cax = ax.matshow(correlation_matrix, cmap='coolwarm')

    # Add colorbar
    fig.colorbar(cax)

    # Set axis ticks and labels
    ticks = np.arange(len(correlation_matrix.columns))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(correlation_matrix.columns, rotation=45, ha="left")
    ax.set_yticklabels(correlation_matrix.columns)

    # Annotate with correlation values
    for i in range(len(correlation_matrix.columns)):
        for j in range(len(correlation_matrix.columns)):
            value = correlation_matrix.iloc[i, j]
            ax.text(j, i, f"{value:.2f}", va='center', ha='center', color='black')

    plt.title('Feature Correlation Matrix')
    plt.tight_layout()
    plt.show()



def augment_data(train_data):
    """
    Augments the dataset by adding synthetic variations.

    Args:
        train_data (list): Original dataset in list format.
        factor (int): The multiplication factor for dataset size.
        noise_level (float): Standard deviation for Gaussian noise.

    Returns:
        pd.DataFrame: Augmented dataset
    """
    df = pd.DataFrame(train_data,
                      columns=["sample", "host", "virus", "ani", "correlation", "fragmentation",
                               "average contig length", 'contig_count', "contamination"])

    augmented_data = df.copy()
    augment_factor = 2
    noise_level_small = 0.01

    for _ in range(augment_factor):  # Generate more data
        temp = df.copy()

        # Add small Gaussian noise to numerical columns
        temp["ani"] = temp["ani"].astype(float) + np.random.normal(0, noise_level_small, len(temp))
        temp["correlation"] = temp["correlation"].astype(float) + np.random.normal(0, noise_level_small, len(temp))
        temp["fragmentation"] = temp["fragmentation"].astype(float) + np.random.normal(0, noise_level_small, len(temp))
        temp["average contig length"] = temp["average contig length"].astype(float) + np.random.normal(0, noise_level_small, len(temp))
        augmented_data = pd.concat([augmented_data, temp], ignore_index=True)

        # Ensure values stay within realistic bounds (-1 or 0 to 1 for probabilities)
        temp["ani"] = temp["ani"].clip(0, 1)
        temp["correlation"] = temp["correlation"].clip(-1, 1)
        temp["fragmentation"] = temp["fragmentation"].clip(0, 1)
        augmented_data = pd.concat([augmented_data, temp], ignore_index=True)
        print(augmented_data)
    return augmented_data


def split_train_test_data(train_data, train_file, test_file, val_file):
    """
    Splits data into train, validation, and test sets after augmentation.

    Args:
        train_data (list): List of original training data.
        train_file (str): Path to save train data.
        test_file (str): Path to save test data.
        val_file (str): Path to save validation data.
    """
    df = pd.DataFrame(train_data,
                      columns=["sample", "host", "virus", "ani", "correlation", "fragmentation",
                               "average contig length", 'contig_count', "contamination"])

    # train (80%) and temp (20%)
    train_set, temp_set = train_test_split(
        df, test_size=0.2, random_state=42, stratify=df["contamination"]
    )

    # temp -> test (10%) and validation (10%)
    test_set, val_set = train_test_split(
        temp_set, test_size=0.5, random_state=42, stratify=temp_set["contamination"]
    )

    augmented_train_set = augment_data(train_set.values.tolist())
    augmented_train_set.to_csv(train_file, index=False)
    test_set.to_csv(test_file, index=False)
    val_set.to_csv(val_file, index=False)

    print(f"Originele Data: {len(df)} rijen")
    print(f"Train Set (geaugmenteerd): {len(augmented_train_set)} rijen (origineel: {len(train_set)})")
    print(f"Validation Set: {len(val_set)} rijen")
    print(f"Test Set: {len(test_set)} rijen")


def main():
    # input
    path_ani_scores = "..\\..\\output_from_scripts\\ani_score_output\\ani_scores.csv"
    path_fragmentation = "..\\..\\output_from_scripts\\fragmentation_output\\fragmentation.csv"
    path_host_virus_correlation = "..\\..\\output_from_scripts\\correlation_analysis_output\\host_virus_correlation_matrix.csv"
    path_contamination = '..\\..\\output_from_scripts\\sample_ids_contaminated_samples\\samples_contaminated.csv'
    path_contig_count = '..\\..\\output_from_scripts\\count_contigs\\count_contigs.csv'
    path_length_contigs = '..\\..\\output_from_scripts\\count_contigs\\length_contigs.csv'

    # output
    train_file = "..\\..\\output_from_scripts\\test_train_data_ml\\train_set.csv"
    test_file = "..\\..\\output_from_scripts\\test_train_data_ml\\test_set.csv"
    val_file = "..\\..\\output_from_scripts\\test_train_data_ml\\val_set.csv"

    read_fragmentation(path_fragmentation)
    contaminated_sample_ids = make_contamination_set(path_contamination)
    train_data = get_train_line(path_ani_scores, path_host_virus_correlation, path_fragmentation,
                                contaminated_sample_ids, path_contig_count, path_length_contigs)
    get_feature_correlation(train_data)
    split_train_test_data(train_data, train_file, test_file, val_file)


main()