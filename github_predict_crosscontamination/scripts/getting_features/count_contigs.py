#!/usr/bin/env python
"""
Author:         Evy Kuenen
Date:           30-3-2025
Functionality:  This script gets the count of contigs of the potential contaminants and the total length of the contigs.

Usage:
    Required directory structure:
                    count_contigs.py needs to be in scripts/getting_features directory
                    potential_contaminats_with_non_foldable_krona_contigs will be retrieved from
                    ..\\..\\output_from_scripts\\potential_contaminants_output\\potential_contaminants_with_non_foldable_krona_contigs.csv
    Required files: potential_contaminants_with_non_foldable_krona_contigs.csv
                    directory with fasta files of all contigs of all samples
    Calling script: "python 3 count_contigs.py"

"""
import csv
from pathlib import Path
from Bio import SeqIO


def parse_contig_file(contig_ids_potential_contaminants):
    """
    Reads the file containing contig IDs and parses each line.

    :param contig_ids_potential_contaminants: Path to the file containing contig IDs.
    :return parsed_data: A list of tuples containing (sample_id, host, splitted_line).
    """
    parsed_data = []
    with open(contig_ids_potential_contaminants, "r", encoding="utf-8", errors="ignore") as file:
        for line in file:
            sample_id, host, splitted_line = parse_file(line)
            if sample_id != 1 and host != 1 and splitted_line != 1:
                parsed_data.append((sample_id, host, splitted_line))
    return parsed_data


def get_contig_sequences(sample_id, virus_contigs_id, path_to_contig_seq):
    """
    Retrieves the sequences of a specific virus within a sample.

    :param sample_id: The sample ID of interest.
    :param virus_contigs_id: List of virus contig IDs.
    :param path_to_contig_seq: Path to the directory containing contig sequences.
    :return contigs_virus_seq, contigs_virus: A set of contig sequences and a list of contig names.
    """
    contigs_virus_seq = set()
    contigs_virus = []

    for file in Path(path_to_contig_seq).iterdir():
        if str(sample_id) in str(file):
            for seq_record in SeqIO.parse(file, "fasta"):
                for contigname in virus_contigs_id:
                    contig_number = contigname.split('_')[-1]
                    contig_number_fasta = seq_record.id.strip().split('_')[-1]
                    if contig_number == contig_number_fasta:
                        contigs_virus_seq.add(str(seq_record.seq).strip())
                        contigs_virus.append(contigname.strip())
    return contigs_virus_seq, contigs_virus


def analyze_virus_contig_associations(contig_ids_potential_contaminants, path_to_contig_seq):
    """
    Analyzes virus-associated contigs and their occurrences across different samples and batches.
    :param contig_ids_potential_contaminants: A dataset containing potential contaminant contig IDs.
    :param path_to_contig_seq: The file path to the contig sequence data.
    :return sample_host_virus_fragmentation: A list of fragmented host-virus relations based on batch comparisons.
    """
    parsed_data = parse_contig_file(contig_ids_potential_contaminants)
    sample_host_virus_count = []
    sample_host_virus_length = []

    for sample_id, host, splitted_line in parsed_data:
        viruses_per_sample = []
        for item in splitted_line:
            cleaned_item = item.replace('\"', '').replace('\'', '').replace(']', '')
            if len(cleaned_item) > 0 and not str(cleaned_item[0]).isdigit():
                splitted_item = cleaned_item.split(',')
                virus = splitted_item[:1]
                virus_contigs_id = splitted_item[1:]

                contigs_virus_seq, contigs_virus = get_contig_sequences(
                    sample_id, virus_contigs_id, path_to_contig_seq)

                viruses_per_sample.append([virus, contigs_virus_seq])

        for virus_seqs in viruses_per_sample:
            seqs_virus = set(virus_seqs[1])
            if len(seqs_virus) >= 1:
                sample_host_virus_count_line, sample_host_virus_length_line = analyze_count_length(
                    seqs_virus, sample_id, virus_seqs, host)
                if sample_host_virus_count_line is not None and sample_host_virus_length_line is not None:
                    sample_host_virus_count.append(sample_host_virus_count_line)
                    print(sample_host_virus_count_line)
                    sample_host_virus_length.append(sample_host_virus_length_line)
                    print(sample_host_virus_length_line)
    return sample_host_virus_count, sample_host_virus_length


def analyze_count_length(seqs_virus, sample_id, virus_seqs, host):
    """
    Count contigs and their length.
    :param seqs_virus: Set of sequences from the sample.
    :param sample_id: The sample ID.
    :param virus_seqs: Virus sequences information.
    :param host: The host information.
    :return: A list containing [sample_id, host, virus, fragmentation_percentage].
    """
    seqs_virus = set(map(str, seqs_virus))
    try:
        sample_host_virus_count_line, sample_host_virus_length_line = count_contigs(sample_id, seqs_virus, virus_seqs,
                                                                                    host)
        return sample_host_virus_count_line, sample_host_virus_length_line
    except:
        print('error calculating fragmentation')
        sample_host_virus_count_line = [sample_id, host, virus_seqs[0], 0]
        sample_host_virus_length_line = [sample_id, host, virus_seqs[0], 0]
        return sample_host_virus_count_line, sample_host_virus_length_line


def match_virus_contigs_on_contigs_fasta(contig, virus):
    """
    Matches a virus with its corresponding contigs from a FASTA file.
    :param contig: A string representation of contig information.
    :param virus: A string representation of a virus identifier.
    :return: A list of virus-contig matches if a match is found, otherwise None.
    """
    virus_contig = contig.replace('\"', '').replace('\'', '').replace(']', '').replace('\\', '').split(',')
    if len(virus_contig) >= 1 and 'contig' in str(virus_contig):
        virus = str(virus).replace('\'', '').replace('[', '').replace(']', '')
        if str(virus_contig[0]).lower() == str(virus).lower():
            return virus_contig


def get_contig_info(contig_line, batch, sample_id):
    """
    Get contigs per virus.
    :param contig_line: A string representing contig data, formatted with batch numbers and sample identifiers.
    :param batch: A tuple containing batch numbers (batch ID and sub-batch ID).
    :param sample_id: The sample ID to exclude from the filtering.
    :return contigs_other_viruses: A list of contigs from other viruses in the same batch, if found.
    """
    contig_line = contig_line.strip()
    if len(contig_line.split('-')) > 2:
        batch_number = contig_line.replace('\"', '').replace('\'', '').split('-')[0], \
            contig_line.replace('\"', '').replace('\'', '').split('-')[1]
        if batch_number == batch and contig_line.split(',')[0] != sample_id:
            contigs_other_viruses = contig_line.split('[')
            return contigs_other_viruses


def filter_sample_right_length(sample_items, file):
    """
    Filters sample IDs to ensure they have the correct format and length.
    :param sample_items: A list containing sample data, where the first item is the sample ID.
    :param file: A string representing the file in which the sample ID should appear.
    :return file: The file if the sample ID matches the expected length criteria, otherwise None.
    """
    if len(sample_items) >= 3:
        sample_items[2] = sample_items[2].split('_')[0]
        sample_items[0] = sample_items[0].replace(' ', '')
        if len(sample_items[0]) == 6 or len(sample_items[0]) == 7:
            sample_items[0] = "".join(filter(str.isdigit, sample_items[0]))
            sample = f"{sample_items[0]}-{sample_items[1]}-{sample_items[2].split('_')[0]}"
            if str(sample) in str(file):
                return file


def count_contigs(sample_id, seqs_virus, virus_seqs, host):
    '''
    Calculate the fragmentation percentage for every virus based on the longest contig sequence virus of a sample
    devided by the longest contig sequence of the batch of that virus.
    :param seqs_batch_virus: list with sequences of all of a same virus found in a batch
    :param seqs_virus: list with sequences of same virus found in sample
    :return fragmentation_percentage: float with percentage of fragmentation
    '''
    try:
        count_contigs = len(seqs_virus)
        length_contigs = 0
        for item in seqs_virus:
            length_contigs = length_contigs + len(item)
        print(f"count for {virus_seqs[0]}: {count_contigs}")
        print(f"length for {virus_seqs[0]}: {length_contigs}")
        sample_host_virus_count_line = [sample_id, host, virus_seqs[0],
                                        count_contigs]
        print(sample_host_virus_count_line)
        sample_host_virus_length_line = [sample_id, host, virus_seqs[0], length_contigs]
        print(sample_host_virus_length_line)
        return sample_host_virus_count_line, sample_host_virus_length_line
    except:
        print('error calculating count')
        sample_host_virus_count_line = [sample_id, host, virus_seqs[0], 0]
        sample_host_virus_length_line = [sample_id, host, virus_seqs[0], 0]
        return sample_host_virus_count_line, sample_host_virus_length_line


def parse_file(line):
    """
    Parses a line from the contig file.

    :param line: A single line from the file.
    :return: (sample_id, host, splitted_line) or (1, 1, 1) if parsing fails.
    """
    line = line.strip()
    splitted_line = line.split('[')
    sample_and_host = splitted_line[0].split(',')
    if len(sample_and_host) > 1:
        sample_id = sample_and_host[0].replace('\"', '')
        host = sample_and_host[1]
        if sample_id is not None and host is not None and splitted_line is not None:
            return sample_id, host, splitted_line
    else:
        sample_id = 1
        host = 1
        splitted_line = 1
        return sample_id, host, splitted_line


def save_count(sample_host_virus_count, count_file):
    """
    Save information about countin csv file.
    :param sample_host_virus_count: list with sample, host, list with viruses and
    their count
    :param count_file: path to file with count of contigs
    """
    try:
        with open(count_file, "w", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerows(sample_host_virus_count)
    except:
        print('cant open file!')


def save_length(sample_host_virus_length_line, length_file):
    """
    Save information about count of contigs in csv file.
    :param sample_host_virus_length_line: list with sample, host, list with viruses and
    their count
    :param length_file: path to file with length of contigs.
    """
    try:
        with open(length_file, "w", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerows(sample_host_virus_length_line)
    except:
        print('cant open file!')


def main():
    # input
    contig_ids_potential_contaminants = "/MolbioStorage/6.Temp/Temp_Evy/potential_contaminants_with_non_foldable_krona_contigs.csv"
    path_to_contig_seq = Path("/MolbioStorage/6.Temp/Temp_Evy/contigs")

    # output
    length_file = "length_contigs.csv"
    count_file = "count_contigs.csv"

    sample_host_virus_count_line, sample_host_virus_length_line = analyze_virus_contig_associations(
        contig_ids_potential_contaminants, path_to_contig_seq)
    save_count(sample_host_virus_count_line, count_file)
    save_length(sample_host_virus_length_line, length_file)


if __name__ == "__main__":
    main()