#!/usr/bin/env python
"""
Author: Evy Kuenen
Date: 26-2-2025
Functionality: This script gets ANI score of the potential contaminants by getting the contigs of the viruses that are
potential contaminants and calculating their ani score with k-mers.

Usage:
    Required directory structure:
                    calculate_ani_score.py needs to be in scripts/getting_features directory
                    potential_contaminats_with_non_foldable_krona_contigs will be retrieved from
                    ..\\..\\output_from_scripts\\potential_contaminants_output\\potential_contaminants_with_non_foldable_krona_contigs.csv
    Required environment: pyani needs to be installed in environment with in the with minmatch=5 and mincluster=10
                    parameters adjusted in anim.pywith minmatch=5 and mincluster=10 parameters adjusted in anim.py
    Required files: potential_contaminants_with_non_foldable_krona_contigs.csv
                    directory with fasta files of all contigs of all samples
    Calling script: "python 3 calculate_ani_score.py"
"""

import os
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import csv


def calculate_ani(seqs_virus, batch_without_sample, log_file):
    """
    Calculate ani score per virus on the whole sample based on kmer size.
    :param seqs_virus: list with all contigs of virus in sample
    :param all_contig_seqs_of_sample: list all contigs in sample
    :return: ani score as float
    """

    if not seqs_virus or not batch_without_sample:
        return "no virus contigs"

    base_dir = r"/MolbioStorage/6.Temp/Temp_Evy/testing_ani/temp/"
    genome_dir = os.path.join(base_dir, "genomes")
    output_dir = os.path.join(base_dir, "pyani_output")

    os.makedirs(genome_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # make fasta files
    try:
        virus_fasta = os.path.join(genome_dir, "virus_sample.fasta")
        with open(virus_fasta, "w") as vf:
            for i, seq in enumerate(seqs_virus):
                vf.write(f">virus_contig_{i}\n{seq}\n")
    except:
        print('cant open file!')

    try:
        sample_fasta = os.path.join(genome_dir, "virus_batch.fasta")
        with open(sample_fasta, "w") as sf:
            for i, seq in enumerate(batch_without_sample):
                sf.write(f">sample_contig_{i}\n{seq}\n")
    except:
        print('cant open file!')

    # PyANI ANIm run
    try:
        result = subprocess.run(
            ["average_nucleotide_identity.py", "-mANIm", "-i", genome_dir, "-o", output_dir, "-v", "-f", "--maxmatch"],
            check=True,
            capture_output=True,
            text=True
        )
        try:
            with open(log_file, 'a') as log:
                log.write("\n=== Nieuwe PyANI Run ===\n")
                log.write(result.stderr)
        except:
            print('cant open file!')
        print("PyANI output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Fout bij uitvoeren van PyANI:", e.stderr)
        return "fout bij pyani uitvoering"

    ani_results_file = os.path.join(output_dir, "ANIm_percentage_identity.tab")
    print("Bestaat het bestand?", os.path.exists(ani_results_file))
    print("Bestanden in output dir:", os.listdir(output_dir))

    if not os.path.exists(ani_results_file):
        return "fout ophalen resultaten"

    ani_df = pd.read_csv(ani_results_file, sep="\t", index_col=0)
    print("PyANI DataFrame:\n", ani_df)

    # get ani scores
    try:
        ani_score = ani_df.loc["virus", "sample"]
        return ani_score if not pd.isna(ani_score) else "fout"
    except KeyError:
        print("Mogelijke fout: Bestandennamen komen niet overeen.")
        return "fout"


def parse_contig_file(contig_ids_potential_contaminants):
    """
    Reads the file containing contig IDs and parses each line.

    :param contig_ids_potential_contaminants: Path to the file containing contig IDs.
    :return parsed_data: A list of tuples containing (sample_id, host, splitted_line).
    """
    parsed_data = []
    try:
        with open(contig_ids_potential_contaminants, "r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                sample_id, host, splitted_line = parse_file(line)
                if sample_id != 1 and host != 1 and splitted_line != 1:
                    parsed_data.append((sample_id, host, splitted_line))
        return parsed_data
    except Exception as e:
        print('cant open file!')
        print(e)


def filter_contigs_in_same_batch(contig_ids_potential_contaminants, batch, sample_id, virus):
    """
    Searches for contigs of the same virus in other samples within the batch.
    :param contig_ids_potential_contaminants: Path to the file with contig IDs.
    :param batch: A tuple containing batch information.
    :param sample_id: The sample ID of interest.
    :param virus: The virus name.
    :return contigs_same_virus_other_sample: A set of contigs from the same virus found in other samples.
    """
    contigs_same_virus_other_sample = set()
    try:
        with open(contig_ids_potential_contaminants, "r", encoding="utf-8", errors="ignore") as file:
            for contig_line in file:
                contigs_other_viruses = filter_on_same_batch(contig_line, batch, sample_id)
                if contigs_other_viruses is not None:
                    for contig in contigs_other_viruses:
                        virus_contig = match_virus_contigs_on_contigs_fasta(contig, virus)
                        if virus_contig is not None:
                            for contig_part in virus_contig[1:]:
                                contigs_same_virus_other_sample.add(contig_part.strip())
        return contigs_same_virus_other_sample
    except:
        print('cant open file!')


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
                        # Ensure that the contig is not already added
                        if contigname not in contigs_virus:
                            contigs_virus_seq.add(str(seq_record.seq).strip())  # Save the sequence
                            contigs_virus.append(contigname.strip())  # Save the contig name
    return contigs_virus_seq, contigs_virus


def get_contigs_of_virus_in_other_samples(contigs_same_virus_other_sample, path_to_contig_seq):
    """
    Retrieves sequences of the virus in other samples within the batch.

    :param contigs_same_virus_other_sample: Set of contigs from the same virus in other samples.
    :param path_to_contig_seq: Path to the directory containing contig sequences.
    :return contigs_same_virus_other_sample_seq: A set of sequences from the virus in other samples.
    """
    contigs_same_virus_other_sample_seq = set()
    for file in Path(path_to_contig_seq).iterdir():
        for contigname in contigs_same_virus_other_sample:
            sample_items = contigname.split('-')
            # not needed if sample numbers (genomescan) are named correctly
            file_path = filter_sample_right_length(sample_items, file)
            if file_path is not None:
                for seq_record in SeqIO.parse(file_path, "fasta"):
                    contig_number = contigname.split('_')[-1]
                    contig_number_fasta = seq_record.id.strip().split('_')[-1]
                    if contig_number == contig_number_fasta:
                        contigs_same_virus_other_sample_seq.add(str(seq_record.seq).strip())
    return contigs_same_virus_other_sample_seq


def analyze_virus_contig_associations(contig_ids_potential_contaminants, path_to_contig_seq):
    """
    Analyzes virus-associated contigs and their occurrences across different samples and batches.

    :param contig_ids_potential_contaminants: A dataset containing potential contaminant contig IDs.
    :param path_to_contig_seq: The file path to the contig sequence data.
    :return sample_host_virus_ani: A list of ani score host-virus relations based on batch comparisons.
    """
    log_file = "pyani.log"
    sample_host_virus_ani = []
    parsed_data = parse_contig_file(contig_ids_potential_contaminants)
    for sample_id, host, splitted_line in parsed_data:
        viruses_per_sample = []
        viruses_per_batch = []
        for item in splitted_line:
            cleaned_item = item.replace('\"', '').replace('\'', '').replace(']', '')
            if len(cleaned_item) > 0 and not str(cleaned_item[0]).isdigit():
                splitted_item = cleaned_item.split(',')
                virus = splitted_item[:1]
                virus_contigs_id = splitted_item[1:]

                batch = sample_id.split('-')[0], sample_id.split('-')[1]
                contigs_same_virus_other_sample = filter_contigs_in_same_batch(
                    contig_ids_potential_contaminants, batch, sample_id, virus)

                contigs_virus_seq, contigs_virus = get_contig_sequences(
                    sample_id, virus_contigs_id, path_to_contig_seq)

                contigs_same_virus_other_sample_seq = get_contigs_of_virus_in_other_samples(
                    contigs_same_virus_other_sample, path_to_contig_seq)
                if contigs_virus_seq is not None and contigs_virus is not None:
                    viruses_per_sample.append([virus, contigs_virus_seq])
                    viruses_per_batch.append([virus, contigs_same_virus_other_sample_seq])

        for virus_seqs in viruses_per_sample:
            seqs_virus = set(virus_seqs[1])
            for batch_virus_seqs in viruses_per_batch:
                seqs_batch_virus = set(batch_virus_seqs[1])

                if len(seqs_virus) >= 1 and virus_seqs[0] == batch_virus_seqs[0] and len(seqs_batch_virus) >= 1:
                    sample_host_virus_ani_line = compare_batch_with_sample(
                        seqs_batch_virus, seqs_virus, sample_id, virus_seqs, host, log_file)
                    if sample_host_virus_ani_line is not None:
                        sample_host_virus_ani.append(sample_host_virus_ani_line)
    return sample_host_virus_ani


def compare_batch_with_sample(seqs_batch_virus, seqs_virus, sample_id, virus_seqs, host, log_file):
    """
    Compares sequences from a batch with those from a sample and calculates ani score.

    :param seqs_batch_virus: Set of sequences from the batch.
    :param seqs_virus: Set of sequences from the sample.
    :param sample_id: The sample ID.
    :param virus_seqs: Virus sequences information.
    :param host: The host information.
    :return: A list containing [sample_id, host, virus, ani score].
    """
    # remove seqs_virus from seqs_batch
    seqs_batch_virus = set(map(lambda x: str(x).strip().lower(), seqs_batch_virus))
    seqs_virus = set(map(lambda x: str(x).strip().lower(), seqs_virus))
    batch_without_sample = seqs_batch_virus.symmetric_difference(seqs_virus)

    if len(batch_without_sample) >= 1:
        try:
            ani_percentage = calculate_ani(batch_without_sample, seqs_virus, log_file)
            print(f"ani score for {virus_seqs[0]}: {ani_percentage}")
            sample_host_virus_ani_line = [sample_id, host, virus_seqs[0],
                                          ani_percentage]
            print(sample_host_virus_ani_line)
            return sample_host_virus_ani_line
        except:
            print('error calculating ani score')
            sample_host_virus_ani_line = [sample_id, host, virus_seqs[0], 0]
            return sample_host_virus_ani_line
    if len(batch_without_sample) == 0:
        print(f"ani score Score for {virus_seqs[0]}: 1, all seqs same")
        sample_host_virus_ani_line = [sample_id, host, virus_seqs[0], 1]
        print(sample_host_virus_ani_line)
        return sample_host_virus_ani_line


def match_virus_contigs_on_contigs_fasta(contig, virus):
    """
    Matches a virus with its corresponding contigs from a FASTA file.

    :param contig: A string representation of contig information.
    :param virus: A string representation of a virus identifier.
    :return: A list of virus-contig matches if a match is found, otherwise None.
    """
    virus_contig = contig.replace('\"', '').replace('\'', '').replace(']', '').replace('\\', '').split(',')
    if len(virus_contig) > 1 and 'contig' in str(virus_contig):
        virus = str(virus).replace('\'', '').replace('[', '').replace(']', '')
        if str(virus_contig[0]).lower() == str(virus).lower():
            return virus_contig


def filter_on_same_batch(contig_line, batch, sample_id):
    """
    Filters contigs belonging to the same batch while excluding the provided sample ID.

    :param contig_line: A string representing contig data, formatted with batch numbers and sample identifiers.
    :param batch: A tuple containing batch numbers (batch ID).
    :param sample_id: The sample ID to exclude from the filtering.
    :return contigs_other_viruses: A list of contigs from other viruses in the same batch, if found.
    """
    contig_line = contig_line.strip()
    if len(contig_line.split('-')) > 2:
        batch_number = contig_line.replace('\"', '').replace('\'', '').split('-')[0], \
            contig_line.replace('\"', '').replace('\'', '').split('-')[1]
        start_sample = contig_line.replace('\"', '').replace('\'', '').split('-')[0]
        if start_sample == '000000':
            contigs_other_viruses = contig_line.split('[')
            return contigs_other_viruses
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
    if len(sample_items) == 1:
        for item in sample_items:
            if '_' in item:
                if '000000-' in str(file):
                    return file


def parse_file(line):
    """
    Parses a line from the contig file. Split contig_ids_potential_contaminants line in host virus and contiglist.
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


def save_ani_scores(sample_host_virus_aniscores, path_ani):
    """
    Save information about ani scores in csv file.
    :param sample_host_virus_aniscores: sample_host_virus_aniscores: list with sample, host, list with viruses and
    their ani score
    """
    try:
        with open(path_ani, "w", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerows(sample_host_virus_aniscores)
    except:
        print("cant open file!")


def main():
    # input
    contig_ids_potential_contaminants = "/MolbioStorage/6.Temp/Temp_Evy/potential_contaminants_with_non_foldable_krona_contigs.csv"
    path_to_contig_seq = Path("/MolbioStorage/6.Temp/Temp_Evy/contigs")

    # output
    path_ani = "ani_scores.csv"
    sample_host_virus_aniscores = analyze_virus_contig_associations(contig_ids_potential_contaminants,
                                                                    path_to_contig_seq)
    save_ani_scores(sample_host_virus_aniscores, path_ani)


main()