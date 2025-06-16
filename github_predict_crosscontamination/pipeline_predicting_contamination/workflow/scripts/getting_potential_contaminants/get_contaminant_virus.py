#!/usr/bin/env python
"""
Author: Evy Kuenen
Date: 26-2-2025
Functionality: This script gets the potential contaminants by getting all viruses per sample
and grouping them per batch. If 2 or more samples in a batch have a virus they get put in a file
with their samplename, hostname, viruses and the count of viruses.

Usage:
    Required directory structure:
                    get_contaminant_virus.py needs to be in scripts/getting_potential_contaminants directory
                    The file describing the krona data needs to be in
                    ..\\..\\output_from_scripts\\getting_sample_host_virus_from_krona_files\\krona_data.csv
                    The path to krona files with contis needs to be in
                    ..\\..\\test_data\\data\\clc_output_with_krona_virus_taxon_contigs_output
                    The path to extra contigs from krona files that were not foldable needs to be in
                    ..\\..\\test_data\\data\\extra_krona_raporten_aanvraag_extra_contigs
    Required files: krona_data.csv
                    files in ..\\..\\test_data\\data\\clc_output_with_krona_virus_taxon_contigs_output
                    files in ..\\..\\test_data\\data\\extra_krona_raporten_aanvraag_extra_contigs
    Calling script: "python 3 get_contaminant_virus.py"
"""
import ast
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path
import re
from bs4 import BeautifulSoup


def make_output_dir(directory_name):
    """
    Makes output directory where output files of finding potential contaminant virus get stored.
    """
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


def group_viruses_by_sample(normalized_data):
    """
    Reads CSV file and groups viruses on batch.
    :param normalized_data: path to CSV file with sample name, host name and viruses found on host.
    :return samples: dictionary with sample names and set of viruses found on them.
    """
    try:
        with open(normalized_data, "r", encoding='utf-8') as filecontent:
            normalized_data_list = []
            for line in filecontent:
                splitted_line = line.split(',')
                normalized_data_list.append(splitted_line)
    except:
        print('cant open file!')
    viruses = set()
    samples = {}

    for row in normalized_data_list:
        if len(row) == 5:
            sample_id, host_name, virus_name, tax_id, rank = row
            viruses.add(virus_name)
            if sample_id not in samples:
                samples[sample_id] = {"viruses": set()}  # make dictionary with sample and their viruses
            samples[sample_id]["viruses"].add(virus_name)
        else:
            sample_id, host_name, virus_name, tax_id, rank, *contigs = row
            viruses.add(virus_name)
            if sample_id not in samples:
                samples[sample_id] = {"viruses": set()}  # make dictionary with sample and their viruses
            samples[sample_id]["viruses"].add(virus_name)

    return samples


def make_list_viruses_per_sample(samples):
    """
    Convert dictionary with sample as key and viruses as value to list.
    :param samples: dictionary with sample names and set of viruses found on them.
    :return: list with sample names and set of viruses found on them.
    """
    filtered_samples = [
        [sample_id, list(samples[sample_id]["viruses"])]
        for sample_id in samples
    ]
    return filtered_samples


def compare_to_batch(filtered_samples):
    """
    Groups samples by batch and identifies potential contaminants based on shared viruses.
    :param filtered_samples: list with sample names and set of viruses found on them.
    :return contaminant_samples: dictionary mapping sample to viruses that appear in at least two samples per batch.
    """
    sample_with_potential_contaminant = []
    grouped_samples = defaultdict(list)
    for sample_id, viruses in filtered_samples:
        parts = sample_id.split('-')
        first_part = parts[0]
        middle_part = parts[1]
        key = (first_part, middle_part)
        grouped_samples[key].append((sample_id, set(viruses)))

    for key, virus_sets in grouped_samples.items():
        virus_count = defaultdict(int)
        for virus_set in virus_sets:
            for virus in virus_set:
                if str(virus).count('{') == 1:  # looks if it is not sample name because virus sets between {}
                    for species in virus:
                        virus_count[species] += 1
        # get only viruses that occur in at least 2 samples
        common_viruses = [virus for virus, count in virus_count.items() if count >= 2]

        for sample_id, viruses in filtered_samples:
            parts = sample_id.split('-')
            first_part = parts[0]
            middle_part = parts[1]
            key_sample = (first_part, middle_part)
            if key == key_sample:
                matched_viruses = [virus for virus in viruses if
                                   virus in common_viruses]

                if matched_viruses:
                    sample_with_potential_contaminant.append([sample_id, matched_viruses])

    contaminant_samples = {}
    for item in sample_with_potential_contaminant:
        sample_name = item[0]
        viruses = item[1]
        if sample_name not in contaminant_samples:
            contaminant_samples[sample_name] = {"viruses": set()}
            contaminant_samples[sample_name]["viruses"].update(viruses)
    return contaminant_samples


def get_samples_potential_contamination_with_host(contaminant_samples, normalized_data):
    """
    Matches samples with potential contaminants to their hosts.
    :param contaminant_samples: dictionary mapping sample to viruses that appear in at least two samples per batch.
    :param normalized_data: path to krona info
    :return sample_contaminants_hosts: list with sample, host and viruses that occur on that host and at least one
    other sample in the batch and the count of how many viruses there are.
    """
    sample_contaminants_hosts = []
    sample_contaminants_virus = []

    for sample in contaminant_samples.items():
        virus_list = []
        sample_name = sample[0]
        viruses = list(sample[1].values())
        for item in viruses:
            virus_list.append(list(item))
        virus_sample = [sample_name, virus_list]
        sample_contaminants_virus.append(virus_sample)

    host_samples = {}
    try:
        with open(normalized_data, "r",
                  encoding='utf-8') as host_data:
            for line in host_data:
                splitted_line = line.split(',')
                sample_name = splitted_line[0]
                hosts = splitted_line[1]
                if sample_name not in host_samples:
                    host_samples[sample_name] = {"hosts": set()}
                host_samples[sample_name]["hosts"].add(hosts)
    except:
        print('cant open file!')

    for sample in host_samples.items():  # match virus and host on sample
        host = list(sample[1].values())
        sample_name = sample[0]
        for host_name in host:
            for line in sample_contaminants_virus:
                if sample_name == line[0]:
                    for virus_list in line[1]:
                        for host in list(host_name):
                            host_virus_sample_item = [sample_name, host, virus_list, len(virus_list)]
                            sample_contaminants_hosts.append(host_virus_sample_item)
    return sample_contaminants_hosts


def filter_duplicates(sample_contaminants_hosts):
    """
    Filters out duplicate samples where the host is marked as 'unknown' and saves sample, host, viruses and count of
    viruses in file.
    :param sample_contaminants_hosts: list with sample, host and viruses that occur on that host and at least one
    other sample in the batch and the count of how many viruses there are.
    """
    count_samples = []
    for row in sample_contaminants_hosts:
        count_samples.append(row[0])

    duplicates = []
    non_duplicates = []
    for item in count_samples:
        if item not in non_duplicates:
            non_duplicates.append(item)
        else:
            duplicates.append(item)

    for item in sample_contaminants_hosts:
        for duplicate in duplicates:
            if item[0] == duplicate:
                if item[1] == 'unknown':
                    sample_contaminants_hosts.remove(item)
    cleaned_sample_contaminants_host = sample_contaminants_hosts
    return cleaned_sample_contaminants_host


def save_potential_contaminats(path_potential_contaminants, cleaned_sample_contaminants_host):
    """
    Save information about potential_contaminants in csv file.

    :param path_potential_contaminants: path to file to save lines with sample, host, viruses that occur in more than
    one sample in batch.
    :param cleaned_sample_contaminants_host: list with sample, host, list with viruses.
    :return:
    """
    try:
        with open(path_potential_contaminants, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["sample", "host", "viruses", "count_viruses"])
            writer.writerows(cleaned_sample_contaminants_host)
    except:
        print('cant open file!')


def add_contigs(path_potential_contaminants, clc_output_with_krona_virus_taxon_contigs_output):
    """
    Add contigs from KRONA data to viruses and save to csv file.

    :param path_potential_contaminants: path to file with potential contaminated samples.
    :param clc_output_with_krona_virus_taxon_contigs_output: path to directory with contig names
    :return contaminants_with_contigs: list with sample, host, list with viruses and
    their contigs.
    """
    path_to_contig_names = Path(clc_output_with_krona_virus_taxon_contigs_output)
    sample_pattern = re.compile(r"(\d{6,}-\d{3,}-\d{3,})")
    contaminants_with_contigs = []
    files_no_contigs = set()
    files_total = set()

    with open(path_potential_contaminants, "r", newline="", encoding="utf-8") as file:
        for line in file:
            line_parts = line.split('[')
            sample_id = line_parts[0].split(',')[0].strip()
            host = line_parts[0].split(',')[1]
            if ']' in line:
                virus_names = line_parts[1].split(']')[0]
                cleaned_virus_names = virus_names.replace("[", "").replace("\\", "").replace("'", "")\
                    .replace('"', "").strip()
                virus_list = [v.strip() for v in cleaned_virus_names.split(',')]
                contig_list = []
                for virus in virus_list:
                    virus_contigs = [[virus]]
                    for year_dir in path_to_contig_names.iterdir():
                        if year_dir.is_dir() and year_dir.name.isdigit():  # loop through year directories
                            for krona_contig_file in year_dir.iterdir():
                                lines_no_contig = 0
                                total_lines_file = 0
                                if krona_contig_file.is_file():  # loop through files
                                    match = sample_pattern.match(krona_contig_file.name)
                                    if match:
                                        extracted_id = match.group(1)  # Extract matched sample ID
                                        if extracted_id == sample_id:
                                            files_total.add(sample_id)  # needs set because loops file for every virus
                                            with krona_contig_file.open("r", encoding="utf-8") as file:
                                                for line_content in file:
                                                    total_lines_file = total_lines_file+1

                                                    if virus.lower() in line_content.lower():
                                                        # even if pep says last part not needed it is needed
                                                        contigs = re.findall(r"\[(.*?)\]", line_content)
                                                        if contigs:
                                                            virus_contigs.extend(contigs)
                                                    if 'N/A' in line_content.strip():
                                                        lines_no_contig = lines_no_contig+1
                                            if (total_lines_file-1) == lines_no_contig:
                                                files_no_contigs.add(sample_id)
                                                virus_contigs.extend(contigs)
                                    if not match:
                                        try:
                                            with open(krona_contig_file, "r", encoding='utf-8') as filecontent:
                                                lines = filecontent.readlines()
                                                for line in lines[1:2]:
                                                    sample = re.search(sample_pattern, line).group(1)
                                                    if sample == sample_id:
                                                        files_total.add(sample_id)
                                                        for line in lines[1:]:
                                                            total_lines_file = total_lines_file + 1
                                                            line_items = line.strip().split("\t")
                                                            virus_krona = line_items[0]
                                                            if virus_krona.lower() == virus.lower():
                                                                contigs = str(line_items[3]).replace('[', '')\
                                                                    .replace(']', '')

                                                                virus_contigs.append(contigs)
                                                                if 'N/A' in line.strip():
                                                                    lines_no_contig = lines_no_contig + 1
                                                        if (total_lines_file - 1) == lines_no_contig:
                                                            files_no_contigs.add(sample_id)
                                                            virus_contigs.extend(contigs)
                                        except:
                                            continue

                    contig_list.extend(virus_contigs)
                contaminants_with_contigs_line = sample_id, host, contig_list
                print(contaminants_with_contigs_line)
                contaminants_with_contigs.append(contaminants_with_contigs_line)
    print(files_no_contigs)
    print("without contigs:", len(files_no_contigs))
    print("total:", len(files_total))
    return contaminants_with_contigs


def save_potential_contaminant_contigs(potential_contaminants_with_contigs, contaminants_with_contigs):
    """
    Save information about save_potential_contaminats_with_foldable_krona_contigs in csv file.
    :param potential_contaminants_with_contigs: path to file to save lines with contigs
    :param contaminants_with_contigs: list with sample, host, list with viruses and
    their contigs.
    """
    try:
        with open(potential_contaminants_with_contigs, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(["sample", "host", "viruses with contigs", "count_viruses"])
            writer.writerows(contaminants_with_contigs)
    except:
        print('cant open file!')


def make_list_from_multiple_dimension_list(virus_data):
    """
    Make from list with multiple dimensions a 2d list.
    :param virus_data: multiple dimensions list of viruses and their contigs.
    :return [virus_data]: 2d list with viruses and their contigs
    """
    if isinstance(virus_data, list):
        return virus_data
    try:
        while isinstance(virus_data, str):
            virus_data = ast.literal_eval(virus_data)  # Keep decoding until it's a list
        if isinstance(virus_data, list):
            return virus_data  # Successfully parsed as a list
    except (ValueError, SyntaxError):
        pass  # If parsing fails, fallback to treating it as a single-item list
    return [virus_data]


def fill_in_last_contigs(path_to_extra_contig_names, potential_contaminants_with_contigs):
    """
    Fill in all contigs of KRONA files that were not foldable and match them with their viruses in every sample and put
    them in a csv file.
    :param path_to_extra_contig_names: path to not foldable KRONA files
    :param potential_contaminants_with_contigs path to file with contaminants with contigs
    """
    lines_all_contigs = []
    sample_ids = []
    with open(potential_contaminants_with_contigs, "r", encoding="utf-8") as file:
        next(file)
        for line in file:
            if 'contig' not in line:
                sample_id, host, virus_data = line.strip().split(',', 2)
                print(host)
                print(sample_id)
                sample_ids.append(sample_id)
                sample_ids.append(sample_id)
                virus_list = make_list_from_multiple_dimension_list(virus_data)
                flattened_virus_list = [item for sublist in virus_list for item in sublist]
                for year_dir in path_to_extra_contig_names.iterdir():
                    if year_dir.is_dir() and year_dir.name.isdigit():
                        if year_dir.name == '2023' or year_dir.name.isdigit() == '2024' \
                                or year_dir.name.isdigit() == '2025':  # loop through year directories
                            for batch_dir in year_dir.iterdir():
                                if batch_dir.is_dir():
                                    krona_paths = [
                                        batch_dir / "5. Externe krona rapporten",
                                        batch_dir / "5. Externe KRONA rapporten",
                                        batch_dir / "5. Klikbare KRONA rapporten",
                                        batch_dir / "5. Klikbare krona rapporten",
                                        batch_dir / "5. Externe KRONA rapporten",
                                        batch_dir / "5. Extern KRONA rapporten",
                                        batch_dir / "5. Extern KRONA rapport",
                                        batch_dir / "5. Externe KRONA rapport",
                                        batch_dir / "5. Extern krona rapporten",
                                        batch_dir
                                    ]
                                    for krona_path in krona_paths:
                                        if krona_path.exists() and krona_path.is_dir():
                                            for sample_dir in krona_path.iterdir():
                                                file_id = sample_dir.name.split('_')[0]
                                                if sample_id == file_id and sample_dir.is_dir():
                                                    members = []
                                                    contigs_of_members = []
                                                    for sample_file in sample_dir.iterdir():
                                                        if sample_file.is_file() and "krona.html" in sample_file.name:
                                                            virus_members = extract_members_from_krona(
                                                                sample_file, flattened_virus_list)
                                                            for virus, member in virus_members.items():
                                                                member_contig = [virus, member]
                                                                members.append(member_contig)
                                                    for sample_file in sample_dir.iterdir():
                                                        if sample_file.is_dir():
                                                            for item in members:
                                                                for member in item:
                                                                    if isinstance(member, list):
                                                                        for item in member:
                                                                            for contig_files in sample_dir.iterdir():
                                                                                if '.files' in contig_files.name:
                                                                                    for contig_file in contig_files.iterdir():
                                                                                        if item in contig_file.name:
                                                                                            with contig_file.open("r", encoding="utf-8") as extra_contig_info:
                                                                                                for line in extra_contig_info:
                                                                                                    line = line.replace('data(\'', '').replace('\\n\\', '').replace('\')', '')
                                                                                                    contigs_of_members.append(line)

                                                    if len(members) > 0:
                                                        member = members[0][0]
                                                        contigs_of_members.insert(0, [member])
                                                        new_line_contigs = [sample_id, host, [contigs_of_members]]
                                                        lines_all_contigs.append(new_line_contigs)
                    if year_dir.is_dir() and year_dir.name.isdigit():
                        if year_dir.name == '2022' or year_dir.name == '2021':
                            for batch_dir in year_dir.iterdir():
                                if batch_dir.is_dir():
                                    file_id = batch_dir.name.split('_')[0]
                                    if sample_id == file_id and batch_dir.is_dir():
                                        members = []
                                        contigs_of_members = []
                                        for sample_file in batch_dir.iterdir():
                                            if sample_file.is_file() and "krona.html" in sample_file.name:
                                                virus_members = extract_members_from_krona(sample_file, flattened_virus_list)
                                                for virus, member in virus_members.items():
                                                    member_contig = [virus, member]
                                                    members.append(member_contig)
                                        for sample_file in batch_dir.iterdir():
                                            if sample_file.is_dir():
                                                for item in members:
                                                    list_members = item[1]
                                                    if isinstance(list_members, list):
                                                        for member in list_members:
                                                            for contig_files in batch_dir.iterdir():
                                                                if '.files' in contig_files.name:
                                                                    for contig_file in contig_files.iterdir():
                                                                        if member in contig_file.name:
                                                                            with contig_file.open("r",
                                                                                                  encoding="utf-8") as extra_contig_info:
                                                                                for line in extra_contig_info:
                                                                                    line = line.replace('data(\'',
                                                                                                        '').replace(
                                                                                        '\\n\\', '').replace('\')', '')
                                                                                    contigs_of_members.append(line)
                                        if len(members) > 0:
                                            member = members[0][0]
                                            contigs_of_members.insert(0, [member])
                                            new_line_contigs = [sample_id, host, [contigs_of_members]]
                                            lines_all_contigs.append(new_line_contigs)

            if 'contig' in line:
                print(line)
                lines_all_contigs.append([line])
    return lines_all_contigs


def save_potential_contaminats_with_foldable_krona_contigs(potential_contaminants_with_non_foldable_krona_contigs, lines_all_contigs):
    """
    Save information about save_potential_contaminats_with_foldable_krona_contigs in csv file.
    :param lines_all_contigs: sample_host_virus_aniscores: list with sample, host, list with viruses and
    their contigs.
    :param potential_contaminants_with_non_foldable_krona_contigs: path to file to save lines with contigs
    """
    try:
        with open(potential_contaminants_with_non_foldable_krona_contigs, "w", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerows(lines_all_contigs)
    except:
        print('cant open file!')


def extract_members_from_krona(file_path, virus_list):
    """
    Get all viruses per sample.
    :param file_path: path to file with virus and contig in KRONA html files
    :param virus_list: list with all viruses
    :return virus_members:  list with all viruses per sample
    """
    with file_path.open("r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, "html.parser")
        virus_members = {}
        for node in soup.find_all("node", recursive=True):
            virus_name = node.get("name")
            if virus_name in virus_list:
                members_tag = node.find_all("members")
                for members in members_tag:
                    members_val = members.find_all("val")
                    for val in members_val:
                        if val.text.strip():
                            if virus_name not in virus_members:
                                virus_members[virus_name] = []
                            virus_members[virus_name].append(val.text.strip())
        return virus_members


def main():
    # input
    normalized_data = sys.argv[1]
    path_to_extra_contig_names = Path(sys.argv[2])
    clc_output_with_krona_virus_taxon_contigs_output = sys.argv[3]
    # output
    directory_name = sys.argv[4]
    path_potential_contaminants = sys.argv[5]
    potential_contaminants_with_non_foldable_krona_contigs = sys.argv[6]
    potential_contaminants_with_contigs = sys.argv[7]

    make_output_dir(directory_name)
    samples = group_viruses_by_sample(normalized_data)
    filtered_samples = make_list_viruses_per_sample(samples)
    contaminant_samples = compare_to_batch(filtered_samples)
    sample_contaminants_hosts = get_samples_potential_contamination_with_host(contaminant_samples, normalized_data)
    cleaned_sample_contaminants_host = filter_duplicates(sample_contaminants_hosts)
    save_potential_contaminats(path_potential_contaminants, cleaned_sample_contaminants_host)
    contaminants_with_contigs = add_contigs(path_potential_contaminants, clc_output_with_krona_virus_taxon_contigs_output)
    save_potential_contaminant_contigs(potential_contaminants_with_contigs, contaminants_with_contigs)
    lines_all_contigs = fill_in_last_contigs(path_to_extra_contig_names, potential_contaminants_with_contigs)
    save_potential_contaminats_with_foldable_krona_contigs(potential_contaminants_with_non_foldable_krona_contigs,
                                                           lines_all_contigs)
    


if __name__ == "__main__":
    main()
