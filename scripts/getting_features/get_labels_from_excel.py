#!/usr/bin/env python
"""
Author: Evy Kuenen
Date: 26-2-2025
Functionality: This script gets the labels for training if a sample was contaminated or not.

Usage:
    Required directory structure:
                    get_labels_from_excel.py needs to be in scripts/getting_features directory
                    The excel files where contamination is described needs to be retrieved from the
                    Molbio\inzendingen\Virologie directory in the molbiostorage and the files to convert the lims
                    number to the sample number needs to be retrieved from the Molbio\Sequentiedata\Genomescan
                    directory.
    Required files: files in Molbio\inzendingen\Virologie
                    files in Molbio\Sequentiedata\Genomescan
    Calling script: "python 3 get_labels_from_excel.py"
"""
import os
import pandas as pd
from pathlib import Path
from openpyxl.packaging.manifest import mimetypes
import re


def make_output_dir():
    """
    Makes output directory where output files of label analysis get stored.
    """
    # naar beneden
    directory_name = "..\\..\\output_from_scripts\\sample_ids_contaminated_samples"
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


def get_all_excel_files(path_excel_dir):
    """
    Gets all KRONA file names from given directory.
    :param path_krona_dir: path to directory with all KRONA files.
    :return krona_files: list with paths to KRONA txt files.
    """
    file_paths = []
    for root, _, files in os.walk(path_excel_dir):
        for file in files:
            if not root.endswith('2013') and not root.endswith('2014') and not root.endswith('2015'):
                file_paths.append(os.path.join(root, file))
    return file_paths


def detect_file_type(file):
    """
    Detect file type of file.
    :param file: file path
    :return: type of file
    """
    mime_type, _ = mimetypes.guess_type(file)
    return mime_type or "unknown"


def make_2d_list(file_paths):
    """
    make a 2d list of content of the excel files from virology.
    :param file_paths: list with paths of excel files
    :return excel_file: content of excel file in 2d list
    """
    excel_file = []
    for file in file_paths:
        try:
            file_type = detect_file_type(file)
            if "officedocument.spreadsheetml" in file_type:  # .xlsx
                df = pd.read_excel(file, sheet_name=None, engine="openpyxl")
            elif "vnd.ms-excel" in file_type:
                df = pd.read_excel(file, sheet_name=None, engine="xlrd")
            elif "vnd.oasis.opendocument.spreadsheet" in file_type:  # .ods (LibreOffice)
                df = pd.read_excel(file, sheet_name=None, engine="odf")
            elif file_type == 'unknown' and file.endswith('.xlsx'):
                df = pd.read_excel(file, sheet_name=None, engine="openpyxl")
            else:
                print(f"Skipping unsupported file: {file} (Detected type: {file_type})")
                continue

            for sheet_df in df.values():
                excel_file.append(sheet_df.values.tolist())

        except Exception as e:
            print(f"Error reading {file}: {e}")
    return excel_file


def get_lims_nmbr(excel_file):
    """
    Retrieve limsnumber from excel file
    :param excel_file: list with content of excel file
    :return lims_row: list with list with every line in excel with lims number.
    """
    lims_row = []  # lims number is item[2] in every row
    for file in excel_file:
        for line in file:
            if str(line[2]).isdigit():
                if int(line[2]) > 10:
                    lims_row.append(line)
    return lims_row


def extract_virus_names(text):
    """
    get virusname if virus is mentioned in sample with contamination.
    :param text: line of file where contamination is mentioned
    :return extracted_data: virus names if found.
    """
    regex_patterns = {
        "woord_na_besmet_met": r"besmet met (\w+) van uit monster",
        "woorden_voor_dat_gedetecteerd_is": r"([\w\s\-]+) dat gedetecteerd is",
        "woorden_na_Wel_voor_gedetecteerd": r"Wel (.*?) gedetecteerd",
        "woorden_voor_is_ook_gedetecteerd": r"([\w\s\-]+) is ook gedetecteerd",
        "woorden_na_this_might_suggest_that": r"this might suggest that (.*?) in sample",
        "woorden_na_gaven_met_voor_De_coverage": r"gaven met (.*?)\. De coverage",
        "woord_voor_4_contigs_gedetecteerd": r"([\w\s\-]+): 4 contigs gedetecteerd",
        "woorden_na_sequenties_van_voor_die_gedetecteerd": r"sequenties van (.*?) die gedetecteerd",
        "woord_voor_zeer_waarschijnlijk_contaminatie": r"([\w\s\-]+) zeer waarschijnlijk contaminatie",
        "woord_voor_gedetecteerd_waarschijnlijk_contaminatie": r"([\w\s\-]+) gedetecteerd, waarschijnlijk contaminatie",
        "woord_voor_waarschijnlijk_contaminatie_uit": r"([\w\s\-]+) zeer waarschijnlijk contaminatie uit",
        "woord_voor_zeer_waarschijnlijk_contaminatie_uit": r"([\w\s\-]+) zeer waarschijnlijk contaminatie uit",
        "woord_voor_waarschijnlijk_contminatie_uit": r"([\w\s\-]+) waarschijnlijk contminatie uit",
        "woord_voor_gedecteerd_5": r"([\w\s\-]+) gedecteerd, 5",
        "woord_voor_waarschijnlijk_contaminatie_klein": r"([\w\s\-]+) waarschijnlijk contaminatie, klein",
        "woord_voor_waarschijnlijk_contaminatie_komt": r"([\w\s\-]+) waarschijnlijk contaminatie, komt",
        "woord_voor_zeer_waarschijnlijk_contaminatie_haakjes": r"([\w\s\-]+) zeer waarschijnlijk contaminatie \(",
        "woord_voor_waarschijnlijk_contaminatie_haakjes": r"([\w\s\-]+) waarschijnlijk contaminatie \(",
        "woord_voor_gedecteerd_zeer_waarschijnlijk": r"([\w\s\-]+) gedecteerd, zeer waarschijnlijk",
        "woorden_na_1_voor_zijn_ook_gedetecteerd_in": r"1\. (.*?) zijn ook gedetecteerd in",
        "woorden_voor_contaminatie_uit_104326": r"([\w\s\-]+) contaminatie uit 104326-095-005",
        "woorden_na_virussen_uit_monster_voor_gedetecteerd": r"virussen uit monster 4630595 (.*?) gedecteerd, "
    }
    extracted_data = {}
    for key, pattern in regex_patterns.items():
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            extracted_data[key] = match.group(1)
    return extracted_data


def get_label_for_lims_nmbr(lims_row):
    """
    find if contamination is mentioned in remarks about sample.
    :param lims_row: list with list with every line in excel with lims number.
    :return lims_with_contamination: list with limsnumbers with contamination confirmed in sample.
    """
    virus_names = {
        'TMV': 'Tobacco mosaic virus',
        'BBWV 1': 'Broad bean wilt virus 1',
        'ToBRFV': 'Tomato brown rugose fruit virus',
        'PepMV': 'Pepino mosaic virus',
        'ToChV': 'Tomato chocolate virus',
        'AMV': 'Alfalfa mosaic virus',
        'CMV': 'Cucumber mosaic virus satellite RNA',
        'ToMV': 'Tomato mosaic virus'
    }
    lims_with_contamination = []
    virus_that_is_contamination = []
    for line in lims_row:
        for item in line:
            if isinstance(item, str):
                if 'contamination' in item.lower() or 'contaminant' in item.lower() \
                        or 'contaminatie' in item.lower() or 'besmet' in item.lower():
                    if 'waarschijnlijk geen' not in item.lower() and not 'geen sprake van contaminatie'in item.lower() \
                            and not 'geen contamniantie' in item.lower() and not 'no contamin' in item.lower() \
                            and not 'verschillende sequenties' in item.lower() \
                            and not 'kruisbesmetting lijkt uitgesloten' in item.lower() \
                            and not 'geen sprake is van een besmetting' in item.lower()\
                            and not 'waarschijnlijk niet om contaminatie' in item.lower()\
                            and not 'contaminatie tussen de monsters uitgesloten' in item.lower()\
                            and not 'contamination is unlikely' in item.lower():
                        lims_with_contamination.append(line)
                        result = extract_virus_names(item)
                        if result:
                            for virus in result.values():
                                cleaned_virus = virus.strip().replace(' compleet genoom gedetecteerd', '')\
                                    .replace('satellite RNA is', '').replace('zeer', '').replace('zee', '')\
                                    .replace('In overleg met Marleen besloten niet nader te analyseren', '')\
                                    .replace('is', '').replace('(2 olaten) ', '').replace('(ToMV)', '')\
                                    .replace('(PMMoV)', '').split(' en ')
                                for item in cleaned_virus:
                                    item = item.strip()
                                    if item in virus_names:
                                        virus_that_is_contamination.append([line[2], virus_names[item]])
    return virus_that_is_contamination, lims_with_contamination


def match_label_lims(virus_that_is_contamination, lims_with_contamination):
    """
    match lims number with label if it is contamination.
    :param virus_that_is_contamination: list with viruses that are contamination in sample.
    :param lims_with_contamination: list with limsnumbers with contamination confirmed in sample.
    :return virus_that_is_contamination: list with viruses that are contamination in sample.
    """
    contaminated_lims_set = set()
    for item in virus_that_is_contamination:
        if len(item) >= 2:
            limsnumber_of_contaminated_virus = item[1]
            contaminated_lims_set.add(limsnumber_of_contaminated_virus)
            for sample in lims_with_contamination:
                limsnumber = sample[2]
                if limsnumber not in contaminated_lims_set and limsnumber not in virus_that_is_contamination:
                    virus_that_is_contamination.append([limsnumber])
                    contaminated_lims_set.add(limsnumber)
    print(virus_that_is_contamination)
    return virus_that_is_contamination


def read_bcf_files( path_limsnumber_to_sample):
    """
    match lims numbers with contamination on sample number
    :param virus_that_is_contamination: list with limsnumbers with contamination confirmed in sample.
    :param path_limsnumber_to_sample: path to bcf files with limsnumber and sample number
    :return lims_number_with_sample: dictionary with sample as key and limsnumber as value
    """
    bcf_excel_file = []
    path_limsnumber_to_sample = Path(path_limsnumber_to_sample)
    for year_dir in path_limsnumber_to_sample.iterdir():
        if year_dir.is_dir() and year_dir.name.isdigit():
            for batch_dir in year_dir.iterdir():
                for BCF_dir in batch_dir.iterdir():
                    if 'bcf' in BCF_dir.name.lower():
                        for sample_file in BCF_dir.iterdir():
                            if sample_file.name.endswith('xls') or sample_file.name.endswith('xlsx'):
                                if sample_file.is_file():
                                    file_type = detect_file_type(sample_file)
                                    if "officedocument.spreadsheetml" in file_type:  # .xlsx
                                        df = pd.read_excel(sample_file, sheet_name=None, engine="openpyxl")
                                        for sheet_df in df.values():
                                            bcf_excel_file.append(sheet_df.values.tolist())
                                    elif "vnd.ms-excel" in file_type:  # Mogelijk .xls, maar kan corrupt zijn
                                        df = pd.read_excel(sample_file, sheet_name=None, engine="xlrd")
                                        for sheet_df in df.values():
                                            bcf_excel_file.append(sheet_df.values.tolist())
                                    elif "vnd.oasis.opendocument.spreadsheet" in file_type:  # .ods (LibreOffice)
                                        df = pd.read_excel(sample_file, sheet_name=None, engine="odf")
                                        for sheet_df in df.values():
                                            bcf_excel_file.append(sheet_df.values.tolist())
                                    elif file_type == 'unknown' and str(sample_file).endswith('.xlsx') \
                                            and not str(sample_file.name).startswith('~$'):
                                        df = pd.read_excel(sample_file, sheet_name=None, engine="openpyxl")
                                        for sheet_df in df.values():
                                            bcf_excel_file.append(sheet_df.values.tolist())
                                    else:
                                        print(f"Skipping unsupported file: {sample_file} (Detected type: {file_type})")
                                        continue
    return bcf_excel_file


def match_sample_and_lims_nmbr(bcf_excel_file):
    """
    match lims numbers with contamination on sample number
    :param bcf_excel_file: lines of bcf files in list
    :return lims_number_with_sample: dictionary with lims number and their sample id.
    """
    lims_number_with_sample = {}

    for item in bcf_excel_file:
        for line in item:
            if len(line) > 4:
                if '-' in str(line[2]):
                    if str(line[2].split('-')[0]).isdigit():
                        try:
                            key = line[2]
                            if str(line[3]).isdigit():
                                value = line[3]
                            elif '-' in line[3]:
                                lims_number = line[3].split('-')[0].strip()
                                value = lims_number
                                if 'vt-' in line[3].lower():
                                    if line[3].split('_')[1].isdigit():
                                        lims_number = line[3].split('_')[1]
                                        value = lims_number
                                    if not line[3].split('_')[1].isdigit():
                                        lims_number = line[3].split('_')[0]
                                        lims_number = lims_number.split(' ')[1]
                                        value = lims_number
                            elif '_' in line[3] and not '(Annico Cove)' in line[3]:
                                if 'ppc' in line[3].split('_')[0].lower():
                                    lims_number = line[3].split('_')[0].split('_')[1]
                                    value = lims_number
                                elif 'hts' in line[3].lower():
                                    if '-' in line[3].lower():
                                        lims_number = line[3].split('-')[1].strip()
                                        value = lims_number
                                    if '2024 rko' in line[3].lower():
                                        lims_number = line[3].split('_')[1]
                                        value = lims_number
                                elif len(line[3].split('_')[0]) > 6:
                                    lims_number = line[3].split('_')[0]
                                    if ' ' in lims_number and lims_number.split(' ')[1].isdigit() and \
                                            len(str(lims_number.split(' ')[1])) > 4:
                                        lims_number = lims_number.split(' ')[1]
                                        value = lims_number
                                    if ' ' in lims_number and lims_number.split(' ')[1].isdigit() and \
                                            len(str(lims_number.split(' ')[6])) > 4:
                                        lims_number = lims_number.split(' ')[6]
                                        value = lims_number
                                    if ' ' in lims_number and not lims_number.split(' ')[1].isdigit():
                                        lims_number = lims_number.split(' ')[0]
                                        value = lims_number

                                    value = lims_number
                                elif 'bkd-' in line[3].lower():
                                    lims_number = line[3].split('_')[0]
                                    lims_number = lims_number.split(' ')[1]
                                    value = lims_number
                                elif line[3].count('_') > 4:
                                    lims_number = line[3].split('_')[3]
                                    value = lims_number
                            elif '(Annico Cove)' in line[3]:
                                lims_number = line[3].split('_')[0]
                                lims_number = lims_number.split(' ')[-1]
                                value = lims_number
                            elif ' ' in line[3]:
                                lims_number = line[3].split(' ')[0]
                                value = lims_number
                            elif 'L.' in line[3]:
                                lims_number = line[3].replace('L.', '')
                                value = lims_number
                            elif 'wag' in line[3].lower():
                                lims_number = line[3].replace('WAG', '')
                                value = lims_number
                            else:
                                value = None
                            if value is not None:
                                if str(value).isdigit():
                                    if key not in lims_number_with_sample:
                                        value.replace('{', '').replace('}', '')
                                        lims_number_with_sample[key] = set()
                                    key.replace('_INS', '').replace('WAG0', '').replace('pad 402 (Annico Cove) ', '')\
                                        .replace('pad 156 (Annico Cove) ', '').replace('\n', '').replace('M ', '')
                                    if len(value) < 6:
                                        print(line)
                                        print(value)
                                        print(key)
                                    lims_number_with_sample[key].add(value)
                        except Exception as e:
                            print(f"Error reading {line}: {e}")
    return lims_number_with_sample


def save_lims_to_sample(lims_number_with_sample, path_sample_to_lims):
    """
    save dictionary with lims numbers and their corresponding sample id.
    :param lims_number_with_sample:  dictionary with lims number and sample id.
    :param path_sample_to_lims: path to file to save dictionary.
    """
    try:
        with open(path_sample_to_lims, "w", encoding="utf-8", newline="") as file:
            for item in lims_number_with_sample:
                lims_nmbr = str(lims_number_with_sample[item]).replace('{\'', '').replace('\'}', '')
                line = f"{item}, {lims_nmbr}\n"
                file.write(line)
    except:
        print('cant open file!')


def change_lims_nmbr_to_sample(lims_number_with_sample, virus_that_is_contamination):
    """
    make list with contaminated samples by changing lims number to sample number and put samples in file.
    :param lims_number_with_sample: lims_number_with_sample: dictionary with sample as key and limsnumber as value
    :param virus_that_is_contamination: list with limsnumbers with contamination confirmed in sample.
    """
    samples_with_contamination = []
    for key, value in lims_number_with_sample.items():
        for contamination_line in virus_that_is_contamination:
            single_value = next(iter(value)) if isinstance(value, set) and len(value) == 1 else value
            if len(contamination_line) == 1:
                if str(contamination_line[0]) == str(single_value):
                    samples_with_contamination.append([key])
            if len(contamination_line) == 2:
                if str(contamination_line[0]) == str(single_value):
                    samples_with_contamination.append([key, contamination_line[1]])
    return samples_with_contamination


def filter_samples_duplicates(samples_with_contamination):
    """
    filter duplicate samples if samples are already in batch with virus.
    :param samples_with_contamination: list with samples that are contamination
    :return: set with samples that are contamination
    """
    cleaned_contamination_samples = []
    samples_already_in_list = set()
    for item in samples_with_contamination:
        if len(item) > 1:
            cleaned_contamination_samples.append(item)
            samples_already_in_list.add(item[0])
        if len(item) == 1:
            for sample in item:
                if sample not in str(samples_already_in_list):
                    samples_already_in_list.add(sample)
                    cleaned_contamination_samples.append(item)
    cleaned_contamination_samples = [tuple(item) for item in cleaned_contamination_samples]
    cleaned_contamination_samples = list(set(cleaned_contamination_samples))
    return cleaned_contamination_samples


def save_contaminated_samples(cleaned_contamination_samples, path_contaminated_samples):
    """
    save samples that are contaminated.
    :param cleaned_contamination_samples: set with contaminated samples and viruses if they are mentioned in bcf file.
    :param path_contaminated_samples: path to file to save contaminated samples.
    """
    try:
        with open(path_contaminated_samples, "w", encoding="utf-8", newline="") as file:
            for item in cleaned_contamination_samples:
                file.write(str(item) + "\n")
    except:
        print('cant open file!')


def main():
    # input
    path_excel_dir = "..\\..\\test_data\\data\\label_data_virology\\Virologie\\Virologie"
    path_limsnumber_to_sample = '..\\..\\test_data\\data\\label_data_virology\\Genomescan'

    # output
    path_contaminated_samples = "..\\..\\output_from_scripts\\sample_ids_contaminated_samples\\samples_contaminated.csv"
    path_sample_to_lims = "..\\..\\output_from_scripts\\sample_ids_contaminated_samples\\sample_to_lims.csv"

    make_output_dir()
    file_paths = get_all_excel_files(path_excel_dir)
    excel_file = make_2d_list(file_paths)
    lims_row = get_lims_nmbr(excel_file)
    virus_that_is_contamination, lims_with_contamination = get_label_for_lims_nmbr(lims_row)
    virus_that_is_contamination = match_label_lims(virus_that_is_contamination, lims_with_contamination)
    bcf_excel_file = read_bcf_files(path_limsnumber_to_sample)
    lims_number_with_sample = match_sample_and_lims_nmbr(bcf_excel_file)
    save_lims_to_sample(lims_number_with_sample, path_sample_to_lims)
    samples_with_contamination = change_lims_nmbr_to_sample(lims_number_with_sample, virus_that_is_contamination)
    cleaned_contamination_samples = filter_samples_duplicates(samples_with_contamination)
    save_contaminated_samples(cleaned_contamination_samples, path_contaminated_samples)


if __name__ == "__main__":
    main()
