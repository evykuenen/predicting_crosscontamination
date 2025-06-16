# !/usr/bin/env python
"""
Author: Evy Kuenen
Date: 23-4-2025
Functionality: This script subsamples files and puts them in other files to get files with contamination.

Usage:
    Required directory structure:
                    make_contamination_files.py needs to be in scripts/in_silico_contamination directory
                    The raw reads fastq files needs to be retrieved from the
                    X:\1.RawData\4.Sequencing\Illumina\External\GenomeScan path on the x server of the NIVIP.
    Required environment: seqtk needs to be installed in this environment
    Required files: files in X:\1.RawData\4.Sequencing\Illumina\External\GenomeScan
    Calling script: "python 3 make_contamination_files.py"
"""
import os
import re
import subprocess
from glob import glob
from collections import defaultdict
from pathlib import Path
import gzip


def make_output_dir(output_dir, contamination_dir):
    """
    Makes output directory where output files of insilico analysis get stored.
    :param output_dir: path to output directory
    :param contamination_dir: path tot contamination directory
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(contamination_dir, exist_ok=True)
        print(f"Directory '{output_dir}' created successfully.")
    except FileExistsError:
        print(f"Directory '{output_dir}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{output_dir}'.")
    except Exception as e:
        print(f"An error occurred: {e}")


def extract_sample_id(filename):
    """
    extract sample id in filename
    :param filename: str name of file
    :return: str of sample id or none
    """
    match = re.search(r"(\d{6,}-\d{3,}-\d{3,})", filename)
    return match.group(1) if match else None


def get_sample_info(raw_data_dir):
    """
    group files if they have the same sample id.
    :param raw_data_dir: raw data fastq files
    :return selected_samples: list with samples with and their fastq data
    :return sample_files: dict with sample id and fastq data
    """
    sample_files = defaultdict(list)

    for fq in glob(os.path.join(raw_data_dir, "*.fq.gz")) + glob(os.path.join(raw_data_dir, "*.fastq.gz")):
        sample_id = extract_sample_id(Path(fq).name)
        print(sample_id)
        if sample_id:
            sample_files[sample_id].append(fq)
        print(fq)

    selected_samples = list(sample_files.keys())
    return selected_samples, sample_files


def make_contamination(selected_samples, output_dir, sample_files, seqtk_path):
    """
    merge and subsample samples to get subsamples of every sample that is 10 percent of the total
    reads.
    :param selected_samples: list with samples with and their fastq data
    :param output_dir: path to output
    :param sample_files: dict with sample id and fastq data
    :param seqtk_path: path to seqtk conda environment
    :return merged_files: dictionary with sample id and information of merged files
    :return subsampled_files: dictionary with sample id and information of subsampled_files
    """
    subsampled_files = {}
    merged_files = {}
    for sample in selected_samples:
        merged_file = os.path.join(output_dir, f"{sample}_merged.fastq")
        print(merged_file)
        subsampled_file = os.path.join(output_dir, f"{sample}_subsampled.fastq.gz")
        print(subsampled_file)

        # Merge files
        print('merging files')
        merges_same_sample_files(merged_file, sample_files, sample)
        merged_files[sample] = merged_file
        # Subsampling
        print('subsampling')
        subsample_files(subsampled_file, merged_file, seqtk_path)
        subsampled_files[sample] = subsampled_file
    return merged_files, subsampled_files


def subsample_files(subsampled_file, merged_file, seqtk_path):
    """
    subsample files so only 10 percent of reads is in file.
    :param subsampled_file: path to subsampled file
    :param merged_file: path to merged file
    :param seqtk_path: path to seqtk
    """
    with open(subsampled_file, 'wb') as out_f:
        lines_file = subprocess.run(["wc", "-l", merged_file], capture_output=True, text=True)
        line_count = int(lines_file.stdout.strip().split()[0])
        reads = line_count // 4  # only one every 4 rows is a read in a fastq
        print(f"Total reads: {reads}")
        subsample_reads = reads // 10  # 10% of reads
        subprocess.run([seqtk_path, "sample", merged_file, str(subsample_reads)], stdout=out_f)


def merges_same_sample_files(merged_file, sample_files, sample):
    """
    merge files of same sample.
    :param merged_file: path to merged_file
    :param sample_files: files for every sample in dict
    :param sample: name of sample
    """
    with open(merged_file, 'wb') as out_f:
        for f in sample_files[sample]:
            if f.endswith(".gz"):
                with gzip.open(f, 'rb') as in_f:
                    out_f.write(in_f.read())
            else:
                with open(f, 'rb') as in_f:
                    out_f.write(in_f.read())  # if files in onedrive this does not work


def merge_sample_with_dif_subsample(contamination_dir, selected_samples, merged_files, subsampled_files):
    """
    merge merged sample fastq with subsample of other sample and save in directory.
    :param contamination_dir:  path to directory contamination needs to be saved
    :param selected_samples: list with samples with and their fastq data
    :param merged_files: dictionary with sample id and information of merged files
    :param subsampled_files: dictionary with sample id and information of subsampled_files
    """
    # sample a is full sample, sample b is subsample
    for sample_a in selected_samples:
        for sample_b in selected_samples:
            if sample_a == sample_b:
                continue  # skip self-sample combination

            file_a = merged_files[sample_a]
            file_b = subsampled_files[sample_b]

            combined_filename = os.path.join(
                contamination_dir, f"{sample_a}_with_{sample_b}.fastq"
            )

            # Merge 2 files
            with open(combined_filename, 'wb') as out_f:
                for f in [file_a, file_b]:
                    with open(f, 'rb') as in_f:
                        out_f.write(in_f.read())

            print(f" {combined_filename} is merged")

            checked_filename = os.path.join(
                contamination_dir, f"{combined_filename}_checked.fastq"
            )
            subprocess.run([
                "repair.sh",
                f"in={combined_filename}",
                f"out={checked_filename}",
                "overwrite=t",
                "-Xmx48000m"
            ])
            print(f" {combined_filename} is repaired")

            splitted_filename_r1 = os.path.join(
                contamination_dir, f"{combined_filename}_splitted_r1.fastq"
            )
            splitted_filename_r2 = os.path.join(
                contamination_dir, f"{combined_filename}_splitted_r2.fastq"
            )
            subprocess.run([
                "reformat.sh",
                f"in={checked_filename}",
                f"out={splitted_filename_r1}",
                f"out2={splitted_filename_r2}",
                "overwrite=t",
                "-Xmx48000m"
            ])
            print("splitted!")


def main():
    raw_data_dir = r"/MolbioStorage/6.Temp/Temp_Evy/RawData_small/"
    output_dir = r"subsampled"
    seqtk_path = "seqtk"  # seqtk needs to be in path (conda install -c bioconda seqtk)
    contamination_dir = os.path.join(output_dir, "in_silico_contaminations")

    make_output_dir(output_dir, contamination_dir)
    selected_samples, sample_files = get_sample_info(raw_data_dir)
    merged_files, subsampled_files = make_contamination(selected_samples, output_dir, sample_files, seqtk_path, )
    merge_sample_with_dif_subsample(contamination_dir, selected_samples, merged_files, subsampled_files)


main()
