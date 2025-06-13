#!/usr/bin/env python
"""
Author:         Evy Kuenen
Date:           20-2-2025
Functionality:  This script finds the correlation between hosts and viruses from KRONA files and stores them in a
                correlation matrix in the output directory.
Usage:
    Required directory structure:
                    get_host_virus_correlation.py needs to be in scripts/getting_features directory
                    krona_data will be retrieved from
                    ..\\..\\output_from_scripts\\getting_sample_host_virus_from_krona_files\\krona_data.csv
    Required files: krona_data.csv
    Calling script: "python 3 get_host_virus_correlation.py"
"""

import os
import csv
import pandas as pd
import numpy as np
from collections import defaultdict


def make_output_dir(directory_name):
    """
    Makes output directory where output files of correlation analysis get stored.
    :param directory_name: path to directory for saving correlation results.
    """
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as different_error:
        print(f"An error occurred: {different_error}")


def read_file(normalized_data):
    """
    read krona information from file and save in 2d list.
    :param normalized_data: path to file with krona information
    :return normalized_data_list: 2d list with krona information
    """
    normalized_data_list = []
    try:
        with open(normalized_data, "r", encoding='utf-8') as filecontent:
            for line in filecontent:
                splitted_line = line.split(',')
                if splitted_line[1] != 'unknown':
                    normalized_data_list.append(splitted_line)
    except:
        print('cant open file!')
    return normalized_data_list


def make_host_virus_table(normalized_data_list):
    """
    Make matrix with occurrence of hosts and viruses per sample where host is known. If a host or a virus occur put 1 in
    matrix otherwise 0. Also save file in output directory.
    :param normalized_data_list: 2d list with samplename, hostname without unknown hosts with one name per host, virus
    name, tax id of virus and level of taxonomy that is species.
    :return hosts: set with hostnames.
    """
    hosts = set()
    hosts_list = []
    viruses = set()
    samples = {}

    for row in normalized_data_list:
        if len(row) == 5:
            sample_id, host_name, virus_name, tax_id, rank = row
            viruses.add(virus_name)
            hosts.add(host_name)
            hosts_list.append(host_name)
            if sample_id not in samples:
                samples[sample_id] = {"hosts": set(), "viruses": set()}  # make dictionary with sample and their host
            samples[sample_id]["hosts"].add(host_name)                   # and viruses.
            samples[sample_id]["viruses"].add(virus_name)
            all_entities = sorted(hosts | viruses)
        if len(row) == 6:
            sample_id, host_name, virus_name, tax_id, contig, rank = row
            viruses.add(virus_name)
            hosts.add(host_name)
            hosts_list.append(host_name)
            if sample_id not in samples:
                samples[sample_id] = {"hosts": set(),
                                      "viruses": set()}  # make dictionary with sample and their host and viruses.
            samples[sample_id]["hosts"].add(host_name)
            samples[sample_id]["viruses"].add(virus_name)
            all_entities = sorted(hosts | viruses)
    print("virus count:", len(viruses))
    return hosts_list, hosts, samples, all_entities, viruses


def save_occurrance_matrix(all_entities, host_virus_occurrence_matrix, samples):
    """
    save for every sample what viruses and hosts occurred.
    :param all_entities: list with all hosts and viruses
    :param host_virus_occurrence_matrix: file with for every sample what viruses and hosts occurred.
    :param samples: name of sample
    """
    try:
        with open(host_virus_occurrence_matrix, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(["Sample"] + all_entities)
            for sample_id, entities in samples.items():
                row = [sample_id] + [1 if entity in (entities["hosts"] | entities["viruses"]) else 0 for entity in
                                     all_entities]
                writer.writerow(row)
    except:
        print('cant open file!')


def get_host_virus_frequencies(samples):
    """
    Bereken hoe vaak een host en een virus voorkomen in samples en hoe vaak ze samen voorkomen.
    :param samples: dictionary van sample_id -> {'hosts': set(), 'viruses': set()}
    :param hosts: set van alle hostnamen
    :return: lijst van tuples: (host, virus, host_count, virus_count, co_occurrence_count)
    """
    host_counts = defaultdict(int)
    virus_counts = defaultdict(int)
    co_occurrence_counts = defaultdict(int)

    for sample in samples.values():
        sample_hosts = sample["hosts"]
        sample_viruses = sample["viruses"]

        for host in sample_hosts:
            host_counts[host] += 1

        for virus in sample_viruses:
            virus_counts[virus] += 1

        for host in sample_hosts:
            for virus in sample_viruses:
                co_occurrence_counts[(host, virus)] += 1
    return host_counts, virus_counts, co_occurrence_counts


def make_dataframe_metadata(hosts_list, host_counts, viruses, virus_counts, co_occurrence_counts):
    """
    combine host, virus, hostcount, viruscount and co occurrence in a dataframe.
    :param hosts_list: list with hosts
    :param host_counts: list with host counts
    :param viruses: list with viruses
    :param virus_counts: list with virus counts
    :param co_occurrence_counts: list with co occurrance count
    :return results: dataframe with combined host, virus, hostcount, viruscount and co occurrence.
    """
    results = []
    for host in hosts_list:
        h_count = host_counts.get(host, 0)
        for virus in viruses:
            v_count = virus_counts.get(virus, 0)
            co_count = co_occurrence_counts.get((host, virus), 0)
            results_line = (host, virus, round(h_count), round(v_count), round(co_count))
            results.append(results_line)
    return results


def get_correlation_host_virus(host_virus_occurrence_matrix, hosts, hosts_list, path_host_virus_correlation_matrix,
                               path_1000_smallest_correlations, path_1000_largest_correlation, samples, viruses):
    """
    Get correlation between host and virus from occurence per sample matrix. Save file with correlations and top 1000
    best correlations.
    param host_virus_occurrence_matrix: 2d list with sample name and occurrence of virus and hosts.
    param hosts: set with hostnames
    """
    dataframe_host_virus = pd.read_csv(host_virus_occurrence_matrix)
    df_features = dataframe_host_virus.set_index(dataframe_host_virus.columns[0])
    correlation_matrix = df_features.corr()
    correlation_matrix.to_csv(path_host_virus_correlation_matrix)

    all_features = set(correlation_matrix.columns)
    viruses_set = all_features - hosts
    print(hosts)
    correlation_flat = correlation_matrix.where(np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool))
    correlation_flat = correlation_flat.stack().reset_index()
    correlation_flat.columns = ["virus", "host", "Correlation"]

    host_virus_correlations = correlation_flat[
        (correlation_flat["virus"].isin(hosts) & correlation_flat["host"].isin(viruses_set)) |
        (correlation_flat["virus"].isin(viruses_set) & correlation_flat["host"].isin(hosts))].copy()

    host_virus_correlations["host"] = host_virus_correlations.apply(
        lambda row: row["virus"] if row["virus"] in hosts else row["host"], axis=1)
    host_virus_correlations["virus"] = host_virus_correlations.apply(
        lambda row: row["host"] if row["virus"] in hosts else row["virus"], axis=1)

    host_counts, virus_counts, co_occurrence_counts = get_host_virus_frequencies(samples)
    frequencies = make_dataframe_metadata(hosts_list, host_counts, viruses, virus_counts, co_occurrence_counts)
    freq_df = pd.DataFrame(frequencies, columns=["host", "virus", "HostCount", "VirusCount", "CoOccurrence"])
    print("count samples:", len(samples))

    combined_df = pd.merge(
        host_virus_correlations,
        freq_df,
        on=["host", "virus"],
        how="left",
        suffixes=("", "_dup")  # dont save double virus or host names
    )

    int_cols = ["HostCount", "VirusCount", "CoOccurrence"]
    for col in int_cols:
        if col in combined_df.columns:
            if pd.api.types.is_float_dtype(combined_df[col]):
                if combined_df[col].dropna().apply(float.is_integer).all():
                    combined_df[col] = combined_df[col].fillna(0).astype(int)

    combined_df = combined_df.drop_duplicates(subset=["host", "virus"])

    top_1000_smallest_host_virus_correlations = combined_df.nsmallest(1000, "Correlation")
    top_1000_smallest_host_virus_correlations.to_csv(path_1000_smallest_correlations, index=False)

    top_1000_host_virus_correlations = combined_df.nlargest(1000, "Correlation")
    top_1000_host_virus_correlations.to_csv(path_1000_largest_correlation, index=False)


def main():
    # input
    normalized_data = "..\\..\\output_from_scripts\\getting_sample_host_virus_from_krona_files\\krona_data.csv"

    # output
    directory_name = "..\\..\\output_from_scripts\\correlation_analysis_output"
    host_virus_occurrence_matrix = "..\\..\\output_from_scripts\\correlation_analysis_output\\host_virus_occurrence_per_sample_matrix.csv"
    path_host_virus_correlation_matrix = "..\\..\\output_from_scripts\\correlation_analysis_output\\host_virus_correlation_matrix.csv"
    path_1000_smallest_correlations = "..\\..\\output_from_scripts\\correlation_analysis_output\\top_1000_smallest_host_virus_correlations.csv"
    path_1000_largest_correlation = "..\\..\\output_from_scripts\\correlation_analysis_output\\top_1000_host_virus_correlations.csv"

    make_output_dir(directory_name)

    normalized_data_list = read_file(normalized_data)
    hosts_list, hosts, samples, all_entities, viruses = make_host_virus_table(normalized_data_list)
    save_occurrance_matrix(all_entities, host_virus_occurrence_matrix, samples)
    get_correlation_host_virus(host_virus_occurrence_matrix, hosts, hosts_list, path_host_virus_correlation_matrix,
                               path_1000_smallest_correlations, path_1000_largest_correlation, samples, viruses)


if __name__ == "__main__":
    main()
