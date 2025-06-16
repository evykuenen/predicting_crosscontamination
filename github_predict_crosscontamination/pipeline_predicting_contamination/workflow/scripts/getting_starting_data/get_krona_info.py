#!/usr/bin/env python
"""
Author: Evy Kuenen
Date: 26-2-2025
Functionality: This script gets from every sample the host and viruses.

path_krona_dir = "..\\..\\test_data\\data\\clc_output_with_krona_virus_taxon_contigs_output\\*"

Usage:
    Required directory structure:
                    get_krona_info.py needs to be in scripts/getting_starting_data directory
                    The krona files with host, virus and contigs needs to be retrieved from the
                    clc_output_with_krona_virus_taxon_contigs_output what is the output of the read_krona_files.py
                    script.
    Required files: files in clc_output_with_krona_virus_taxon_contigs_output
    Required environment: conda install xlrd
    Calling script: "python 3 get_krona_info.py"
"""
import glob
import mimetypes
import os
import re
import csv
import string
import sys
import pandas as pd
import warnings


def make_output_dir(dir_name):
    """
    Makes output directory where output files of KRONA analysis get stored.
    """
    try:
        os.mkdir(dir_name)
        print(f"Directory '{dir_name}' created successfully.")
    except FileExistsError:
        print(f"Directory '{dir_name}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{dir_name}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


def get_all_krona_files(path_krona_dir):
    """
    Gets all KRONA file names from given directory.
    param path_krona_dir: path to directory with all KRONA files.
    return krona_files: list with paths to KRONA txt files.
    """
    krona_dirs = glob.glob(path_krona_dir)
    krona_files = []
    
    for directory in krona_dirs:
        files = glob.glob(os.path.join(directory, "*", "*.txt"))
        print(files)
        for item in files:
             if os.path.isfile(item):
                  krona_files.append(item)

    print(krona_files)
    return krona_files


def get_host_name(filename):
    """
    Get the names of the host from the filename and remove parts that are not part of the organism name and remove
    organism names that are not species.
    param filename: list with names of KRONA files.
    :return: name of organism or 'Unknown'
    """
    removed_word = {"batch", "blad", "zaden", "wortels", "zaad", "rnaseq", "rna",
                    "dnaseq", "placenta", "vrucht", "van", "txt", "va", "seeds",
                    "leaf", "toetsplant", "baby", "fruit", "healty", "healthy", "wortel",
                    'bent', 'qui', "rca", "ppc", "pc", "sp.", "spp", "var", 'extract', 'subsp', 'sp', 'dna', 'vir',
                    'pcc', 'm', 'nieuw', 'kroontjes', "gevriesdroogd", "bla", "mengmonster", "rko", "peul", 'back',
                    'batch.txt', '(rca)', 'bkd', '(vrucht)', '(zaad)', '(blad)', 'ca', 'cgn', '(in', 'vivo)', 'extra',
                    'data', 'nakt', 'p', 'ins', "(jong", "blad)", 'c', 'n', '(syst)', 's', 'calyx', 'sampled',
                    '(kroonslipje)', '(kroonslipjes)', 'sub', 'v', '(mashua)', '(kroonslip)', '(apocynaceae)',
                    '(wortel)', 'mix', 'buitenkant', 'vruchtjes', '(calyx)', 'paul', 'wild',
                    }
    pattern = r"""  
        _([A-Za-z]+(?:[_ ]+[A-Za-z]+)*)\s |  
        _([A-Za-z]\.?[A-Za-z]+(?:[_ ]+[A-Za-z]+)*)\b |   
        _(?:\d+_)?([A-Za-z][A-Za-z_. ]+\([A-Za-z_ ]+\)|[A-Za-z][A-Za-z_ ]+) |  
        _(?:\d+[-_]?\d+[-_]?\d+_)?([A-Z]?[a-z]\.?\s?[a-z]+(?:\sva\s[A-Za-z]+)?(?:\s[A-Za-z]+)*)   |
        \d+_\d+\s([A-Z][a-z]+(?:\s[a-z]+)?) |
        \d+_\d+\s([A-Za-z]+\.[A-Za-z]+(?:\s[a-z]+)?(?:\sva\s[a-z]+)?) 
    """
    compiled_pattern = re.compile(pattern, re.VERBOSE)  # Compile regex with re.VERBOSE to ignore spaces
    match = re.search(compiled_pattern, filename)

    if match and 'ppc' not in filename.lower():
        for group in match.groups():
            if group:
                name = group.replace("_", " ")
                words = name.split()
                filtered_words = [word for word in words if word.lower() not in removed_word]
                clean_name = " ".join(filtered_words).replace('\"', '')
                return clean_name if clean_name else "unknown"
    if 'ppc' not in filename.lower():
        return "unknown"


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


def match_lims_to_host(excel_file):
    """
    Retrieve limsnumber from excel file
    :param excel_file: list with content of excel file
    :return lims_row: list with list with every line in excel with lims number.
    """
    removed_word = {"batch", "blad", "zaden", "wortels", "zaad", "rnaseq", "rna", 'wortel', 'subs', 'moeten', '26)',
                    "dnaseq", "placenta", "vrucht", "van", "txt", "va", "seeds", '(seeds)', '(sepals)', '(blad)\nrna,',
                    "leaf", "toetsplant", "baby", "fruit", "healty", "healthy", 'kroonslipjes\nrna', 'dieven\nrna',
                    "rca", "ppc", "pc", "sp.", "spp", "var", 'extract', 'subsp', 'sp', 'dna', 'vir', 'zaden\nrna', 'de',
                    'pcc', 'm', 'nieuw', 'kroontjes', "gevriesdroogd", "bla", "mengmonster", "rko", "peul", 'back',
                    'batch.txt', '(rca)', 'bkd', '(vrucht)', '(zaad)', '(blad)', 'ca', 'cgn', '(in', 'vivo)', 'extra',
                    'data', 'nakt', 'p', 'ins', "(jong", "blad)", 'c', 'n', '(syst)', 's', 'calyx', 'sampled', 'wb,',
                    '(kroonslipje)', '(kroonslipjes)', 'sub', 'v', '(mashua)', '(kroonslip)', '(apocynaceae)', 'pcr',
                    '(wortel)', 'mix', 'buitenkant', 'vruchtjes', '(calyx)', 'paul', 'wild', 'spp.', '(zaad)\nrna',
                    'restant', 'welke', 'opgestuurd', 'was', 'voor', 'hts,', 'hts', 'bcf', '105447-073', '(vir', 'wk',
                    '49,', '(rna,', '-80)', '(van', 'plant)', '(van', 'nvwa', 'terrein)', '(kroontjes)', '(bladnerven)',
                    'glut', 'rna,', 'geisoleerd', 'door', 'naktuinbouw,', 'bakje', '-80', 'wb', 'tracering', '2023)',
                    'geisoleerd', 'door', 'naktuinbouw', 'bakje', 'komk', 'en', 'extract.', 'Ligt', 'in', '(stond',
                    'geïsoleerd', 'ligt', '(blad)\nrna', 'al', '+', '(rna', 'weeklijst', '18)', 'samples', 'toetsing',
                    '(knol)', '(dna', '28)', 'dna)', '(va', 'bladnerven)', 'blad)\nlet', 'op:', 'naktuinbouw)', 'zaad,',
                    'gepooled', 'worden', '8]', 'week', '16', '(va', 'bent)', '(mengmonster)', 'backup', 'of23)', '(2',
                    '(leaf?)', '-81', 'blad,', 'blad\n', '(chilipeper)', 'opslag', "springs'", "baby'", 'leaf\nrna',
                    'blad\nrna', 'var.', 'nwva', 'tobrfv', 'lijst', '(kroontjes)', "'lady", "rosetta'", 'tobamovirus',
                    'qui,', '[va', 'vrucht]', 'zaad\nrna', '-80.', '(seeds)\nrna,', '(zaad)\nrna,', 'knol', 'dat',
                    '(blad).', 'subsp.', 'malmeanum', 'andigena', '(pcr', '14)', 'leaf\nRNA', '(pstvd)', 'pap2', '->',
                    'lade', 'vrucht.', '(blad),', '\n(wine', '3)', '2)', '1)', 'al,', 'wk43', 'zie', 'dieven', '(h28))',
                    '42)', '[blad]', '(pooled', 'leaf)', '(leaf)', '(peul)', 'uit', '(extract', 'loof)', '(h25))',
                    'v.a', '11)', '[vrucht]', '-80.\nstond', 'op', '2023', 'enkele', 'plant', 'p1)', 'kroonslipjes',
                    'blad\ntracering', 'naktuinbouw\n(stond', 'op', '52)', '(blad)\ngewas', 'is', 'onzeker', '50)',
                    '(blad)\ndna', 'samples', 'toetsing', 'lade', 'blad??\nrna', '2021-108', 't/m', '2021-112', '-20?',
                    '2021-098', '2021-107', '2021-097', '2021-088', '2021-087', '2021-078', 'pap', 'tom', '(back-up',
                    'doos', '22)', '(vrucht)\nrna', '104326-080', 'mag', 'gebruikt', '49)', 'geisoleerd:', 'pic', '8,',
                    'bcf104326-154', '(originele', 'dna)\ndna', '(blue', 'bird)', '-', '48)\n', 'moet', 'nog', '\n(2',
                    'worden)\n', '(met', 'symptomen)', '(vrucht)\nrna', 'kan', 'gebruikt', 'seeds\nrna', 'wortel', 'ct',
                    'leaf\nrna', '41)', '40)', "lady'", "bird'", "bride'", "angel'", '(burgundy', 'glow)', '-p1', 'p1,',
                    'lokaal', '(al', 'lokale', 'bent', 'symptomen]', '(monster', 'nog', 'opzoeken', 'glow"', ',', '31',
                    '-80,', 'vrucht\nrna', 'nivip,', '37\n', '(blad+nerf)', 'p1\nRNA', '\'', '(atropurpurea)', 'loof',
                    'rneasy\nf-mol-132-002', 'm&w', 'vivo', 'collectie)\nrna', 'geisoleerd.', 'stond', '32.', 'vir,',
                    'met', 'cucurbit', 'chlorotic', 'yellows', 'virus', 'ccyv', '(wk', '30-2-tobrfv)', '30-tobrfv)',
                    '26-tobrfv)', '(uit', 'collectie)', 'worden)\n', 'p1\nrna', '&', '\n', 'naktuinbouw.', 'I-MOL-061',
                    'calyx\nrna', 'ingezonden', 'bkd)', 'geisoleerd\nstond', '23)', 'bks.', 'zakjes)', '32)', 'zaad)',
                    'liggen', '1,', 'witte', 'doosje', 'zaden.', 'nakt.', '2,', '(rna)', '(floeem)', '(bladnerf)', '?',
                    'blad.', 'ages.', 'nrc', '11d', 'blad.', 'naktuinbouw.tobrfv', 'bakje)', '9b', 'hybride', '(?)',
                    'hongarije', '48)', '(zie', '44)', 'ook', 'deze', 'rna\nrna', '25x', 'blad.', 'tpo', 'inzendingen',
                    'p1.', 'geisoleerd,', 'floeem.', 'kroonblaadjes', '-20', '2021', '104326-058.', 'real-time', 'ToCV',
                    '2020', 'v.a.', '(blad', 'geeaxtraheerd)', '(wortels)', 'zaad\ntracering', 'PAC)', 'f-mol-068-009',
                    '(zaad),', 'geisoleerd,', 'rna-extract', '13-3', '(vruchten)', 'worden)\n', 'bcf-104326-040',
                    'vaatbundels.', 'herkomst', 'china,', 'israel,', 'molbio,', "'paul", "wild'", 'isolatie', 'seeds.',
                    'uitgevoerd-', 'kelkblad', '(28-9-2020)', 'zaad\n\nopslag', '7a', '"vir', 'niet', 'weggooien",',
                    '(10', 'vruchten)', 'schubben', 'qui', 'getoetst', 'molbio', 'ingestuurd', 'geïsoleerd,', '28.84',
                    'opsturen', 'naar', 'genomescan,', 'dus', 'onder', 'bcf-104326-028-008', 'f-mol-071-005\npospi2',
                    'molbio\n20210506_', '', 'project', '104326-023-015', 'bent/p1', '-20,', 'f-mol-071-004\npospi2',
                    '.', 'als', 'hts_20210310,', '104326-006', 'batchnummer', '104326-023)', 'hesp', '(10', 'vruchten)',
                    'stengel', 'kroontje', 'bakje,', 'kelkblad\nrna', 'kroontjes\nrna', 'zaden,', 'molbio012', '(10x',
                    'molbio013', 'molbio014', 'molbio015', 'molbio016', 'v10_vaccin', '(zaad)\nsub', 'geïsoleerd)',
                    '100ul', '(rna-isolate', 'maalzakjes', 'verdunning)',
                    }
    lims_row = []
    lims_with_host = {}
    for file in excel_file:
        for line in file:
            if str(line[2]).isdigit():
                if int(line[2]) > 10:
                    lims_row.append(line)
                    lims_nmbr = line[2]
                    host_item = line[4]
                    host_list = re.split(r'\s+', host_item)
                    host = [word for word in host_list if word.lower().strip() not in removed_word
                            and not word.strip().isdigit()]
                    host = " ".join(host).replace('\"', "").strip()
                    if len(host) > 0 and not host.isdigit():
                        host = ''.join(host)
                        lims_with_host[str(lims_nmbr)] = host
                    if len(host) == 0 or host.isdigit():
                        host = 'unknown'
                        lims_with_host[str(lims_nmbr)] = host

    return lims_with_host


def read_krona_files(krona_files, sample_host, sample_host_logboeken):
    """
    Read content of all KRONA files and get file with samplename, hostname, virusname, taxid of virus and level of
    taxonomy.
    param krona_files: list with paths to KRONA txt files.
    return krona_data: 2d list with samplename, hostname, virusname, taxid of virus and level of taxonomy.
    """
    krona_data = []
    print(krona_files)
    for file in krona_files:
        filename = str(os.path.basename(file).strip())
        id_pattern = r"(\d{6,}-\d{3,}-\d{3,})"  # Matches only sample name
        id_match = re.search(id_pattern, filename)

        if id_match:
            krona_line = batch_match(id_match, filename, sample_host, sample_host_logboeken, file)
            krona_data.extend(krona_line)

        if not id_match:
            try:
                with open(file, "r", encoding='utf-8') as filecontent:
                    lines = filecontent.readlines()
                    for line in lines[1:2]:
                        sample = re.search(id_pattern, line).group(1)
                        if sample in sample_host:
                            organism_name = sample_host[sample]
                            for line in lines[1:]:
                                krona_line = [sample, organism_name] + line.strip().split("\t")
                                krona_data.append(krona_line)
                        if sample in sample_host_logboeken and sample not in sample_host:
                            organism_name = sample_host_logboeken[sample]
                            for line in lines[1:]:
                                krona_line = [sample, organism_name] + line.strip().split("\t")
                                krona_data.append(krona_line)

                        if sample not in sample_host and sample not in sample_host_logboeken:
                            organism_name = 'unknown'
                            for line in lines[1:]:
                                krona_line = [sample, organism_name] + line.strip().split("\t")
                                krona_data.append(krona_line)
            except:
                print('did not find sample name or host name in:', file)
    return krona_data


def batch_match(id_match, filename, sample_host, sample_host_logboeken, file):
    """
    if whole sample is found but not in front of back of sample name get host name.
    :param id_match: sample id str
    :param filename: name of file str
    :param sample_host: dictionary with sample and host from bcf files
    :param sample_host_logboeken: dictionary with sample and host from logboeken virology
    :param file: krona file from sample
    :return krona_line: list with sample name and host name.
    """
    samplename = id_match.group(1)
    organism_name = get_host_name(filename)

    if organism_name is None or organism_name == 'unknown':
        samplename = samplename.strip()
        if samplename in sample_host:
            organism_name = sample_host[samplename]
        elif samplename in sample_host_logboeken:
            organism_name = sample_host_logboeken[samplename]

    krona_lines = []
    try:
        with open(file, "r", encoding='utf-8') as filecontent:
            next(filecontent)  # Skip header
            for line in filecontent:
                if organism_name is not None:
                    krona_line = [samplename, organism_name] + line.strip().split("\t")
                    krona_lines.append(krona_line)
    except Exception as e:
        print(f'Error processing file {file}: {e}')
    return krona_lines


def match_lims_on_sample(lims_with_host, path_lims_to_sample):
    """
     This function reads a CSV file containing sample IDs and LIMS numbers. It uses a dictionary
    (`lims_with_host`) that maps LIMS numbers to host names, and builds a new dictionary that maps
    sample IDs to their respective host names.
    :param lims_with_host: A dictionary where the keys are LIMS numbers (as strings) and the values are host names.
    :param path_lims_to_sample: The path to a CSV file where each row contains a sample ID and its corresponding
    LIMS number.
    :return: A dictionary where the keys are sample IDs and the values are host names matched using the LIMS number.

    """
    sample_host = {}

    try:
        with open(path_lims_to_sample, "r", encoding='utf-8') as filecontent:
            next(filecontent)
            for line in filecontent:
                splitted_line = line.split(',')
                if len(splitted_line) >= 2:
                    lims_number = splitted_line[1].replace('\n', '').strip().replace('{', '').replace('}', '')
                    sample_id = splitted_line[0].strip()
                    if str(lims_number) in lims_with_host:
                        host = lims_with_host[lims_number]
                        sample_host[sample_id] = host
        return sample_host
    except:
        print('did not find sample name or host name')


def get_all_excel_files(path_excel_dir):
    """
    Gets all KRONA file names from given directory.
    :param path_excel_dir: path to directory with all KRONA files.
    :return: krona_files: list with paths to KRONA txt files.
    """
    file_paths = []
    for root, _, files in os.walk(path_excel_dir):
        for file in files:
            if not root.endswith('2013') and not root.endswith('2014') and not root.endswith('2015'):
                file_paths.append(os.path.join(root, file))
    return file_paths


def filter_krona_known_host_species(krona_data):
    """
    Filter on level of taxonomy virus needs to be classified.
    :param krona_data: 2d list with samplename, hostname, virus name, tax id of virus and level of taxonomy.
    :return: 2d list with samplename, hostname, virus name, tax id of virus and
    level of taxonomy that is species.
    """
    filtered_krona_data_species = []
    for row in krona_data:
        if row is not None and len(row) > 4 and row[4] == 'species':
            filtered_krona_data_species.append(row)

    return filtered_krona_data_species


def find_host_if_unknown(path_lims_to_sample, path_extra_hosts):
    removed_word = {"nan", "bevroren", "gh+", "vreemde", '?', "2022-053", "ingezonden", 'extract', "diamantino", "aan",
                    "rna", "uit", "zaden", "cultivar", "onbekend", "cv", "div", "ve/ibad/p", "seviocard", "gewas",
                    "'blue", "Wish'", "zie", "foto.", "betreft", "C.", "spp", "niet", "27-10", "rivolo", "let", "op,",
                    "inzend/TPO", "formulier", "correct", "herkomst", "spanje,", "drie", "vruchten.", "zijn", "op",
                    "volledig", "doorgekleurd.", "losse", "oeganda", "vruchten,", "enkele", "lijken", "groene/gele",
                    "cv", "bronze", "beauty", "\n", "va", "p1", "pap", "bianca", "zaad", "hi-power", "zwarte", "plek",
                    "courgette", 'sp', "ruby", "glow", "happy", "returns", "dit", "is", "gecheckt", "bij", "roy", "gh",
                    "smit", "van", "de", "naktuinbouw", "‘marina’", "‘speciosus’", "ins-19-29037", "ins-19-29039", "en",
                    "gleaming", "gold", "swahili", "indian", "chief", "joanieke/marleen", "zit", "hier", "de", "hts",
                    "tymo", "in", "of", "een", "ander", "virus?", "zie", "ook", "5998981", "contaminatie", "heel",
                    "mogelijk", "contaminatie", "de", "kas)", "va", "tpo", "baby", "-", "van", "veen", "gele", "groene",
                    "met", "symptomen", "zonder", "alkemade", "pier/robert", "waarschijnlijk", "zonale", "pva", "pvs",
                    "cip314909.279", "gar", "11-02", "her-45.1", "red", "wing", "gewas", "blue", "bird", "p1", "blad",
                    "k2024-px1374", "robert", "bcf", "105976-023", "[op", "stond", "ipomoea", "naktuinbouw", "groot",
                    "[grote", "pot", "virus", "beeld/", "klein", "potje", "gezonde", "stek", "courgettes", "uiteinde",
                    "kleine", "stekjes", "bruinverkleuring", "rondom", "nerven.", "echt", "virus", "verdacht", "maar",
                    "wel", "leuk", "om", "eens", "te", "kijken", "wat", "er", "deze", "combi", "komt", "kenia", "ent",
                    "vlekken", "op", "vruchten].", "kleine", "rode", "pepertjes", "verschillende", "vlekjes.", "goed",
                    "doorgekleurd", "licht", "i", "rna", "grote", "beeld", "nerven", "israel", "zaden", "nextstrain",
                    "info", "opgeslagen", "t:\\pd\\nivip\\virologie\\qs", "nl\\2019_tobrfv_tomaat\\toetsing", "voor",
                    "ct", "china", "rossen", "seeds", "cata28", "csp1325", "sub", "naktuinbouw", "vlekjes", "sommige",
                    "spanje", "cgn", "4049929959651", "vruchten", "lichtgroene", "chlorotische", "vrucht", "hevig",
                    "andere", "twee", "chlorotische", "vlekkerigheid", "vruchten", "die", "slecht", "kringen", "die",
                    "deels", "ingezonken", "info", "over", "telefonisch", "contact", "monster", "costa", "rica", "necr",
                    "gedeeld", "myc", "samengesteld", "blad", "deelbladeren", "chl", "tot", "planten/bladeren",
                    "necr", "vlekjes", "super", "extract", "buffer", "nederland", "steeltje", "onregelmatig", "vlekjes",
                    "pleksgewijs", "bruine", "stippen/vlekken", "vooral", "kringachtig", "oppervlakkig", "virologisch",
                    "verspreid", "twijfel", "zekerheid",
                    }
    warnings.simplefilter("ignore")

    sample_host_logboeken = {}
    sample_to_lims = read_sample_to_lims(path_lims_to_sample)

    files = glob.glob(os.path.join(path_extra_hosts, "*"))
    file_content = make_2d_list(files)

    for line in sample_to_lims:
        sample = line.split(',')[0].strip()
        lims_number = sample_to_lims[sample]
        for excel_line in file_content:
            for item in excel_line:
                if lims_number in str(item):
                    if len(item) > 8:
                        if '2019' in str(item[0]) or '2018' in str(item[0]) or '2017' in str(item[0]) \
                                or '2016' in str(item[0]):
                            host = str(str(item[6]) + " " + str(item[7]))
                        if '2020' in str(item[0]):
                            host = str(str(item[7]))
                        if '2024' in str(item[0]):
                            host = str(str(item[9]) + " " + str(item[10]))
                        if '2021' in str(item[0]) or '2022' in str(item[0]) or '2023' in str(item[0]) \
                                or '2025' in str(item[0]):
                            host = str(str(item[8]) + " " + str(item[9]))

                        cleaned_host_words = [
                            word.strip(string.punctuation).replace('’', '').replace(',', '') for word in host.strip()
                            .replace('\n', ' ').replace('\xa0', ' ').replace(',', ' ').replace('\"', "").split(' ')
                            if word.lower().strip(string.punctuation) not in removed_word
                               and not is_number(word.lower().strip(string.punctuation))
                               and 'rko' not in word.lower()
                               and '2023' not in word.lower()
                               and '2025' not in word.lower()
                        ]
                        if cleaned_host_words:
                            if len(cleaned_host_words) > 0:
                                cleaned_host = " ".join(cleaned_host_words)
                                if cleaned_host.strip() is not None and len(cleaned_host) > 3:
                                    sample_host_logboeken[sample] = cleaned_host
    return sample_host_logboeken


def read_sample_to_lims(path_lims_to_sample):
    """
    read file with samples and lims number translation
    :param path_lims_to_sample: path to file
    :return sample_to_lims: dictionary with sample and lims translation
    """
    sample_to_lims = {}
    with open(path_lims_to_sample, "r", encoding='utf-8') as filecontent:
        for item in filecontent:
            splitted_line = item.strip().split(',')
            if len(splitted_line) >= 2:
                sample_id = splitted_line[0].strip()
                lims_number = splitted_line[1].strip().replace('{', '').replace('}', '')
                sample_to_lims[sample_id] = lims_number
    return sample_to_lims


def is_number(word):
    """
    check if word is a number
    :param word: string of word that needs to be checked
    :return: true if string is number or false if string is not number
    """
    try:
        float(word)
        return True
    except ValueError:
        return False


def filter_phages_and_non_plant_virus(filtered_krona_data_species):
    """
    Filter on phages and corona viruses out of dataset.
    :param filtered_krona_data_species: 2d list with samplename, hostname, virus name, tax id of virus and
    level of taxonomy that is species.
    :return: 2d list with samplename, hostname, virus name, tax id of virus and
    level of taxonomy that is species.
    """
    filtered_phages_plant_virus = []
    for line in filtered_krona_data_species:
        if 'phage' not in line[2].lower() and 'corona' not in line[2].lower():
            filtered_phages_plant_virus.append(line)
    return filtered_phages_plant_virus


def filter_synonyms_host(filtered_phages_plant_virus):
    """
    Normalize host names so there is only one name for every host and save file with all used hosts.
    param filtered_krona_data: 2d list with samplename, hostname, virus name and
    level of taxonomy that is species.
    return normalized_data:  filtered_krona_data with one name per host.
    """
    name_mapping = {
        "tomaat": "solanum lycopersicum",
        "s lycopersicum": "solanum lycopersicum",
        "s.lycoperiscum": "solanum lycopersicum",
        "s. betaceum": "solanum betaceum",
        "s.lycopersicum": "solanum lycopersicum",
        "s. lycopersicum": "solanum lycopersicum",
        "solanum_lycopersicum": "solanum lycopersicum",
        "solanum  lycopersicum": "solanum lycopersicum",
        "s. lycoperosicum": "solanum lycopersicum",
        "n. benthamiana": "nicotiana benthamiana",
        "s. aethiopicum": "solanum aethiopicum",
        "lycopersicum": "solanum lycopersicum",
        "s.lycoperiscum p1": "solanum lycopersicum",
        "s.lycopersicum n.benthamiana": "solanum lycopersicum",
        "solanum lycpersicum": "solanum lycopersicum",
        "lycopersicum tofbv": "solanum lycopersicum",
        "lycopersicum   tofbv": "solanum lycopersicum",
        "solanum psuedocapsicum": "solanum pseudocapsicum",
        "solanum macrocarpum": "solanum macrocarpon",
        "solanum_tuberosum": "solanum tuberosum",
        "solanum  tuberosum": "solanum tuberosum",
        "solanum  tuberosum ramos": "solanum tuberosum",
        "solanum tubersom p1": "solanum tuberosum",
        "tuberosum": "solanum tuberosum",
        "alium tuberosum" : "allium tuberosum",
        "tulipa p1": "tulipa",
        "tulipa  royal ten": "tulipa",
        "tulipa_royal_ten": "tulipa",
        "n.occidentalis": "nicotiana occidentalis",
        "hibiscus rosa-sinensis": "hibiscus rosa sinensis",
        "cucumis  melo": "cucumis melo",
        "capsicum  annuum": "capsicum annuum",
        "cucumis  sativus": "cucumis sativus",
        "solanum tuberosum p1": "solanum tuberosum",
        "solanum tubersum p1": "solanum tuberosum",
        "s. tuberosum": "solanum tuberosum",
        "s. tuberosum p1": "solanum tuberosum",
        "s. tuberosum chenopodium quinoa": "solanum tuberosum",
        "s. tuberosum nicotiana occidentalis": "solanum tuberosum",
        "aardappel": "solanum tuberosum",
        "solanum lycoperisum": "solanum lycopersicum",
        "solanum lycopersiucm": "solanum lycopersicum",
        "s tuberosum": "solanum tuberosum",
        "solanum tubersom": "solanum tuberosum",
        "beaucarnia": "beaucarnea",
        "cucurbita  pepo": "cucurbita pepo",
        "ipomaea batata": "ipomaea batatas",
        "ipomeao batatas": "ipomaea batatas",
        "ipomae batatas": "ipomaea batatas",
        "s betaceum": "solanum betaceum",
        "s.betaceum": "solanum betaceum",
        "c. annuum": "capsicum annuum",
        "paprika": "capsicum annuum",
        "hibiscus": "hibiscus rosa sinensis",
        "chenopodiumam album": "chenopodium album",
        "curcurbita pepo": "cucurbita pepo",
        "cucurbita pepo cucumis sativus": "cucurbita pepo",
        "hoya suntice": "hibiscus syriacus",
        "ipomoea batatas": "ipomaea batatas",
        "ipomoae batatas": "ipomaea batatas",
        "hibiscus rosa chinensis": "hibiscus rosa sinensis",
        "physalis peruviana wb": "physalis peruviana",
        "prunus serrulata prunus avium": "prunus serrulata",
        "solanum melonga": "solanum melongena",
        "tropaeolum tuberosia": "tropaeolum tuberosum",
        "ajuga_reptans": "ajuga reptans",
        "lelie": "lillium",
        "c. pepo": "cucurbita pepo",
        "unknown": 'unknown',
        "lilllium": "lillium",
        "vaccium": "vaccinium",
        "asperge": "asparagus officinalis",
        "calibrachoa petchoa en petunia": "calibrachoa petchoa petunia",
        "c.annuum": 'capsicum annuum',
        "c.sativus": "cucumis sativus",
        "beta": "beta vulgaris",
        "capsicum anuum": "capsicum annuum",
        "pastinaak": "pastinaca sativa",
        "roos": "rosa",
        "lactuca  sativa": "lactuca sativa",
        "pier": "pieris japonica",
        "callibrachoa, petunia": "calibrachoa petunia",
        "calibrachoa, petchoa petunia": "calibrachoa, petunia",
        "komkommer": "cucumis sativus",
        "lolium_perenne": "lolium perenne",
        "lonicera japonica aureoretuculata": "lonicera japonica aureoreticulata",
        "p1": "nicotiana occidentalis",
        "physalis peruviana qui tomaat": "physalis peruviana",
        "solanum lycopersicum_zaad": "solanum lycopersicum",
        "calibrachoa, petunia": "calibrachoa petunia",
        "cucurbita pepo p1": "cucurbita pepo",
        "euphorbia palcherima": "euphorbia pulcherrima",
        "hybiscus syriacus": "hibiscus syriacus",
        "physalis peruviana tomaat": "physalis peruviana",
        "ipomoea batatas bellevue": "ipomoae batatas bellevue"
    }
    normalized_data = []
    for line in filtered_phages_plant_virus:
        hostname = line[1]
        lower_hostname = hostname.lower()
        if lower_hostname in name_mapping:
            line[1] = name_mapping.get(lower_hostname, lower_hostname)
            normalized_data.append(line)
        if lower_hostname not in name_mapping:
            if not lower_hostname.isdigit() and len(lower_hostname) > 3:
                line[1] = lower_hostname.strip().replace('\'', '').replace(',', '')
                normalized_data.append(line)

    hostsset = set()
    for line in normalized_data:
        hostsset.add(line[1])
    hostsetlist = sorted(list(hostsset))
    return hostsetlist, normalized_data


def save_hosts_in_krona_files(hostsetlist, path_hosts):
    """
    Save set with host names in csv file.
    :param hostsetlist: set with host names.
    :param path_hosts: path to file where hosts are saved
    """
    with open(path_hosts, "w", encoding="utf-8") as f:
        for item in hostsetlist:
            f.write(item + '\n')


def save_normalized_hosts_viruses(normalized_data, path_normalized_data):
    """
    Save samplename, host name, virus name and contigs in csv file.
    :param normalized_data: 2d list with samplename, host name, virus name and contigs.
    :param path_normalized_data: path where normalized data is saved
    """
    with open(path_normalized_data, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(normalized_data)


def main():
    # input
    path_krona_dir = sys.argv[1]
    path_lims_to_sample = sys.argv[2]
    path_excel_dir = sys.argv[3]
    path_extra_hosts = sys.argv[4]
    
    # output
    dir_name = sys.argv[5]
    path_normalized_data = sys.argv[6]
    path_hosts = sys.argv[7]

    make_output_dir(dir_name)

    file_paths = get_all_excel_files(path_excel_dir)
    excel_file = make_2d_list(file_paths)

    lims_with_host = match_lims_to_host(excel_file)
    sample_host = match_lims_on_sample(lims_with_host, path_lims_to_sample)
    sample_host_logboeken = find_host_if_unknown(path_lims_to_sample, path_extra_hosts)

    krona_files = get_all_krona_files(path_krona_dir)
    krona_data = read_krona_files(krona_files, sample_host, sample_host_logboeken)

    filtered_krona_data_species = filter_krona_known_host_species(krona_data)
    filtered_phages_plant_virus = filter_phages_and_non_plant_virus(filtered_krona_data_species)
    hostsetlist, normalized_data = filter_synonyms_host(filtered_phages_plant_virus)

    # save files
    print(normalized_data)
    save_hosts_in_krona_files(hostsetlist, path_hosts)
    save_normalized_hosts_viruses(normalized_data, path_normalized_data)


if __name__ == "__main__":
    main()
