# -*- coding: utf-8 -*-

import yaml # read config file
import sys
import os
import pathlib
import subprocess

import prepare_and_check

def bam2pbc(bam, pbc):

    cmd = 'bedtools genomecov -ibam "{}" -d > "{}"'.format(bam, pbc)
    try:
        out = subprocess.run([cmd], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n{}'.format(cmd))


def mapping(mapping_dir, fasta_path, read_paths, insert_size, bam_path):

    if type(read_paths) != list:
        read_paths = [read_paths]

    cmd_build = 'bowtie2-build "{}" {}assembly'.format(fasta_path, mapping_dir)

    if len(read_paths) == 2:
        cmd_mapping =  'bowtie2 --sensitive -I {} -X {} -a -p 16 -x {}assembly -1 {} -2 {} -S {}mapping.sam'.format(str(insert_size-200), str(insert_size+200), mapping_dir, read_paths[0], read_paths[1], mapping_dir)
    elif len(read_paths) == 1:
        cmd_mapping =  'bowtie2 --sensitive -a -p 16 -x {}assembly -U {} -S {}mapping.sam'.format(mapping_dir, read_paths[0], mapping_dir)
    else:
        print('Error: unexpected number of read paths detected for coverage set (max allowed: 2; given: {}). Please recheck your input.'.format(str(len(read_paths))))

    cmd_view = 'samtools view -b {}mapping.sam | samtools sort -o {}'.format(mapping_dir, bam_path)

    try:
        out = subprocess.run([cmd_build], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n{}'.format(cmd_build))
    try:
        out = subprocess.run([cmd_mapping], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n{}'.format(cmd_mapping))


def process_coverage(cfg):

    if cfg.compute_coverage:
        for cov_set, status in cfg.cov_set_exists.items():
            if status == 'pbc_paths':
                continue

            elif status == 'bam_paths':
                bam_path = cfg.bam_paths.get(cov_set)
                pbc_path = cfg.pbc_paths.get(cov_set)
                bam2pbc(bam_path, pbc_path)

            elif status == 'read_paths':
                mapping_dir = cfg.output_path+'mapping_files_{}/'.format(str(cov_set))
                pathlib.Path(mapping_dir).mkdir(parents=True, exist_ok=True)

                read_paths_set = cfg.read_paths.get(cov_set)
                insert_size_set = int(cfg.insert_size.get(cov_set))
                bam_path = cfg.bam_paths.get(cov_set)
                pbc_path = cfg.pbc_paths.get(cov_set)


                mapping(mapping_dir, fasta_path, read_paths_set,
                        insert_size_set, bam_path)
                bam2pbc(bam_path, pbc_path)

            else:
                print("Error: status of coverage set {} unidentified".format(cov_set))


def main():

    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(config_path)

    process_coverage(cfg)


if __name__ == '__main__':
    main()
