#!/usr/bin/env python
"""
Date: 22-2-2025
Functionality: This script gets all fasta files on the CLC server from the samples.

Usage:
    Required directory structure:
                    getting_all_contig_seqs_from_clc.py needs to be in scripts/getting_starting_data directory
                    The fasta files needs to be retrieved from the
                    Databases/CLCArchive/1.Diagnostics_Collections/Virus_Viroids/ in the CLC server
    Required files: files in Databases/CLCArchive/1.Diagnostics_Collections/Virus_Viroids/
    Calling script: "python 3 getting_all_contig_seqs_from_clc.py"
"""

import shutil
import os
import re


def save_results(file_path, output_dir):
    """
    Save files to output dir.
    :param file_path: file that needs to be saved.
    :param output_dir: directory where file needs to be saved.
    """
    os.makedirs(output_dir, exist_ok=True)
    file_name = os.path.basename(file_path)
    destination_path = os.path.join(output_dir, file_name)
    shutil.copy(file_path, destination_path)


def main(base_dir, output_dir):
    """
    Iterates through directories to process HTML files and extract virus data.
    :param base_dir: Base directory containing year-wise subdirectories with HTML files.
    :param output_dir: directory where file needs to be saved.
    """
    year_re = re.compile(r'\d{4}$')
    folder_re = re.compile(r'^\d+-\d+')
    denovo_re = re.compile(r'.+denovo', re.IGNORECASE)
    structure_1 = ['2021', '2022', '2023', '2024', '2025']
    structure_2 = ['2016', '2020', '2017', '2018', '2019', '2021']

    for year in os.listdir(base_dir):
        print(year)
        if year_re.match(year):
            if year in structure_1:
                os.makedirs(year, exist_ok=True)
                year_path = os.path.join(base_dir, year)
                for batch in os.listdir(year_path):
                    if folder_re.match(batch):
                        batch_path = os.path.join(year_path, batch)
                        for denovo_dir in os.listdir(batch_path):
                            if denovo_re.match(denovo_dir):
                                denovo_path = os.path.join(batch_path, denovo_dir)
                                for sample in os.listdir(denovo_path):
                                    sample_path = os.path.join(denovo_path, sample)
                                    if os.path.isdir(sample_path):
                                        for file_name in os.listdir(sample_path):
                                            if 'chunked' in file_name and 'sampled' not in file_name:
                                                print(file_name)
                                                file_path = os.path.join(sample_path, file_name)
                                                if file_path:
                                                    save_results(file_path, output_dir)
        if year_re.match(year):
            if year in structure_2:
                os.makedirs(year, exist_ok=True)
                year_path = os.path.join(base_dir, year)
                for batch in os.listdir(year_path):
                    if not batch.startswith('.'):
                        batch_path = os.path.join(year_path, batch)
                        for sample_dir in os.listdir(batch_path):
                            sample_dir_path = os.path.join(batch_path, sample_dir)
                            if os.path.isdir(sample_dir_path):
                                for file_name in os.listdir(sample_dir_path):
                                    if 'chunked' in file_name and 'sampled' not in file_name:
                                        print(file_name)
                                        file_path = os.path.join(sample_dir_path, file_name)
                                        if file_path:
                                            save_results(file_path, output_dir)


if __name__ == "__main__":
    BASE_DIR = "../data/testing/Databases/CLCArchive/1.Diagnostics_Collections/Virus_Viroids/"
    OUTPUT_DIR = "\\test_data\\match_contig_seq_with_potential_contaminant\\Databases\\CLCData\\output"
    main(BASE_DIR, OUTPUT_DIR)
