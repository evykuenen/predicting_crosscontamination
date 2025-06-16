#!/usr/bin/env python
"""
Date: 22-2-2025
Functionality: This script reads all krona .hmtl files following modern clc folder structure and extracts node
information.

Usage:
    Required directory structure:
                    read_krona_files.py needs to be in scripts/getting_starting_data directory
                    The krona files needs to be retrieved from the
                    /Databases/CLCData/1.Diagnostics_Collections/Virus_Viroids/ directory on the CLC server of the NIVIP
    Required files: files in /Databases/CLCData/1.Diagnostics_Collections/Virus_Viroids/
    Calling script: "python 3 read_krona_files.py"
"""

from bs4 import BeautifulSoup
import os
import re
import sys


def process_html(file_path):
    '''
    Parses an HTML file and extracts all 'node' elements with the attribute 'name' set to 'Viruses'.

    :param file_path: Path to the HTML file.
    :return: List of 'node' elements.
    '''
    try:
        with open(file_path) as html:
            soup = BeautifulSoup(html, 'html.parser')
        return soup.find_all("node", {"name": "Viruses"})
    except:
        print('cant open file!')


def extract_node_data(viruses):
    '''
    Extracts name, taxon,contig and rank from virus-related 'node' elements

    :param file_path: Path to the HTML file.
    :param viruses: List of 'node' elements containing virus data.
    :return: List of tuples containing (name, taxon, rank).
    '''
    namelist = []
    taxalist = []
    ranklist = []
    contiglist = []

    def recursive_extract(node):
        '''
        Recursively extracts name, taxon, contig and rank information from a given node.

        :param node: A BeautifulSoup node element.
        '''
        name = node.get('name')
        if name:
            namelist.append(name)

        taxon = node.find('taxon')
        if taxon and taxon.find('val', href=True):
            taxalist.append(taxon.find('val')['href'])
        else:
            taxalist.append(None)

        contig = node.find('members')
        if contig is not None:
            contig_names = []
            vals_elements = contig.find_all('vals')

            for vals in vals_elements:
                val_texts = [val.text.strip() for val in vals.find_all('val') if val.text.strip()]
                contig_names.extend(val_texts)
            if contig_names:
                contiglist.append(contig_names)
            else:
                contiglist.append(None)
        else:
            contiglist.append(None)

        rank = node.find('rank')
        if rank and rank.find('val'):
            ranklist.append(rank.find('val').text)
        else:
            ranklist.append(None)

        for child in node.find_all('node', recursive=False):
            recursive_extract(child)

    for virus in viruses:
        recursive_extract(virus)

    return list(zip(namelist, taxalist, ranklist, contiglist))


def save_results(sample_name, year, data, output_dir):
    '''
    Saves extracted virus data to a text file within a specified year directory.

    :param sample_name: Name of the sample.
    :param year: Year of the data collection.
    :param data: List of tuples containing (name, taxon, rank).
    '''
    year_dir = os.path.join(output_dir, year)
    os.makedirs(year_dir, exist_ok=True)
    file_path = os.path.join(year_dir, f"{sample_name}.txt")
    try:
        with open(file_path, "w") as f:
            f.write("Name\tTaxon\tRank\tContig\n")
            for name, taxon, rank, contig in data:
                f.write(f"{name}\t{taxon or 'N/A'}\t{rank or 'N/A'}\t{contig or 'N/A'}\n")
    except:
        print('cant open file!')


def main(base_dir, output_dir):
    '''
    Iterates through directories to process HTML files and extract virus data.

    :param base_dir: Base directory containing year-wise subdirectories with HTML files.
    '''
    year_re = re.compile(r'\d{4}$')
    folder_re = re.compile(r'^\d+-\d+')
    denovo_re = re.compile(r'.+denovo', re.IGNORECASE)
    denov0_2_re = re.compile(r'.+de novo', re.IGNORECASE)

    for year in os.listdir(base_dir):
        if year_re.match(year):
            if year == '2023' or year == '2022' or year == '2021' or year == '2024' or year == '2025' or year == '0000':
                os.makedirs(year, exist_ok=True)
                year_path = os.path.join(base_dir, year)

                for batch in os.listdir(year_path):
                    if folder_re.match(batch):
                        batch_path = os.path.join(year_path, batch)
                        for denovo_dir in os.listdir(batch_path):
                            if denovo_re.match(denovo_dir) or denov0_2_re.match(denovo_dir):
                                denovo_path = os.path.join(batch_path, denovo_dir)
                                for sample in os.listdir(denovo_path):
                                    sample_path = os.path.join(denovo_path, sample)
                                    print(sample_path)
                                    print("\n")
                                    if os.path.isdir(sample_path):
                                        for file_name in os.listdir(sample_path):
                                            print(file_name)
                                            if file_name.endswith(".html") and 'virus' in file_name.lower():                                      
                                                print(file_name)
                                                file_path = os.path.join(sample_path, file_name)
                                                viruses = process_html(file_path)

                                                result = extract_node_data(viruses)
                                                if result:
                                                    save_results(sample, year, result, output_dir)
                                                    print(f"Results for sample '{sample}' saved to {year}/{sample}.txt")
                                    if os.path.isfile(sample_path):
                                        if sample_path.endswith(".html") and 'virus' in sample_path.lower():
                                            viruses = process_html(sample_path)
                                            result = extract_node_data(viruses)
                                            if result:
                                                save_results(sample, year, result)
                                                print(f"Results for sample '{sample}' saved to {year}/{sample}.txt")
        if year_re.match(year):
            if year == '2016' or year == '2020' or year == '2017' or year == '2018' or year == '2019' \
                    or year == '2021' or year == '2022':
                os.makedirs(year, exist_ok=True)
                year_path = os.path.join(base_dir, year)

                for batch in os.listdir(year_path):
                    if not batch.startswith('.'):
                        batch_path = os.path.join(year_path, batch)
                        for sample_dir in os.listdir(batch_path):
                            sample_dir_path = os.path.join(batch_path, sample_dir)
                            if os.path.isdir(sample_dir_path):
                                for file_name in os.listdir(sample_dir_path):
                                    if file_name.endswith(".html") and not file_name.startswith('.'):
                                        sample = re.search(r"\d+-\d+-\d+", sample_dir)
                                        file_path = os.path.join(sample_dir_path, file_name)
                                        viruses = process_html(file_path)
                                        result = extract_node_data(viruses)
                                        if result:
                                            save_results(sample.group(), year, result)
                                            print(f"Results for sample '{sample.group()}' saved to {year}/{sample.group()}.txt")


if __name__ == "__main__":
    # input
    BASE_DIR = sys.argv[1]
    
    # output
    output_dir = sys.argv[2]
    main(BASE_DIR, output_dir)


