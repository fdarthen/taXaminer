#!/usr/bin/env python

"""Prepare coverage

Prepares coverage data for further processing. Raw reads are mappend
and mapping files (BAM) are converted to per base coverage (PBC) files

Expects prepared config file with information regarding available
coverage information and corresponding paths to find the files
(preparation of config by prepare_and_check.py)
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

import logging
import sys
import pathlib
import subprocess

# package specific modules
import prepare_and_check


def bam2pbc(bam, pbc):
    """Convert BAM file to PBC.

    Uses bedtools genomecov with option -d to convert mapping
    information in BAM format to per base coverage file (PBC)
    Option -d reports read depth at each position

    Args:
      bam(str): path to BAM file
      pbc(str): path to PBC file (output)
    """

    # FIXME: check this here (add specific output filtering)
    # --> probably error because samtools view command was skipped in mapping()
    cmd = 'bedtools genomecov -ibam "{}" -d > "{}"'.format(bam, pbc)
    out = subprocess.run([cmd], shell=True, capture_output=True)
    if out.returncode != 0:
        logging.error('converting BAM to PBC failed:\n{}'.format(cmd))
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()


def mapping(mapping_dir, fasta_path, read_paths, min_insert_set, max_insert_set, read_orientation_set, bam_path, threads):
    """Mapping reads to assembly using Bowtie2.

    Uses Bowtie2 to map reads to reference assembly. Options used:
      -a: every alignment is reported (multiple alignments enabled)
      --sensitive: sensitivity mode
      -I/-X: insert size (lower and upper limit)

    Aligns paired and single end reads sets based on how many paths
    there are in read_paths

    Args:
      mapping_dir(str): path to output directory of mapping files
      fasta_path(str): path to FASTA file of assembly
      read_paths(str,list): (list of) path(s) to read file(s)
      insert_size(int): insert size for paired end reads, else None
      bam_path(str): path to BAM file
    """

    # TODO: check if SAM file exists at default location
    if type(read_paths) != list:
        read_paths = [read_paths]

    ## prepare commands
    cmd_build = 'bowtie2-build "{}" {}assembly'.format(fasta_path, mapping_dir)
    # depending on number of paths in list read_paths different cmd
    # for mapping of paired end or single end data is used
    if len(read_paths) == 2:
        cmd_mapping =  'bowtie2 --sensitive -I {} -X {} --{} -a -p {} -x {}assembly -1 {} -2 {} -S {}mapping.sam'.format(min_insert_set, max_insert_set, read_orientation_set, threads, mapping_dir, read_paths[0], read_paths[1], mapping_dir)
    elif len(read_paths) == 1:
        cmd_mapping =  'bowtie2 --sensitive -a -p {} -x {}assembly -U {} -S {}mapping.sam'.format(threads, mapping_dir, read_paths[0], mapping_dir)
    else:
        print('Error: unexpected number of read paths detected for coverage set (max allowed: 2; given: {}). Please recheck your input.'.format(str(len(read_paths))))
    cmd_view = 'samtools view -b {}mapping.sam | samtools sort -o {}'.format(mapping_dir, bam_path)

    ## run commands
    # build index
    out = subprocess.run([cmd_build], shell=True, capture_output=True, check=True)
    if out.returncode != 0:
        logging.error('failed building index for mapping:\n{}'.format(cmd_build))
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()
    # mapping
    out = subprocess.run([cmd_mapping], shell=True, capture_output=True, check=True)
    if out.returncode != 0:
        logging.error('failed aligning reads to assembly:\n{}'.format(cmd_mapping))
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()
    # convert SAM to BAM
    out = subprocess.run([cmd_view], shell=True, capture_output=True, check=True)
    if out.returncode != 0:
        logging.error('sorting and converting SAM to BAM failed:\n{}'.format(cmd_view))
        logging.error('Error message:\n'+out.stderr.decode())
        sys.exit()

def process_coverage(cfg):
    """Prepare coverage information for further analysis

    Uses information in config object to perform computations for
    further use of coverage information. Calls
      mapping()
      bam2pbc()

    Args:
      cfg: Config object with config parameters
    """

    if cfg.compute_coverage:
        # cfg.cov_set_exists: {index: status}
        # for every set of coverage data there is
        for cov_set, status in cfg.cov_set_exists.items():
            # check the status
            if status == 'pbc_paths':
                # coverage data already in PBC format
                continue

            elif status == 'bam_paths':
                # coveage data in BAM format -> convert to PBC
                bam_path = cfg.bam_paths.get(cov_set)
                pbc_path = cfg.pbc_paths.get(cov_set)
                bam2pbc(bam_path, pbc_path)

            elif status == 'read_paths':
                # coverage data as read files -> map + convert
                # create directory to store mapping files
                mapping_dir = cfg.output_path+'mapping_files_{}/'.format(str(cov_set))
                pathlib.Path(mapping_dir).mkdir(parents=True, exist_ok=True)

                read_paths_set = cfg.read_paths.get(cov_set)
                min_insert_set = int(cfg.min_insert.get(cov_set))
                max_insert_set = cfg.max_insert.get(cov_set)
                read_orientation_set = cfg.read_orientation.get(cov_set)
                bam_path = cfg.bam_paths.get(cov_set)
                pbc_path = cfg.pbc_paths.get(cov_set)

                num_threads = cfg.threads if cfg.threads != 'auto' else 1


                mapping(mapping_dir, cfg.fasta_path, read_paths_set,
                        min_insert_set, max_insert_set, read_orientation_set, bam_path, num_threads)
                bam2pbc(bam_path, pbc_path)

            else:
                print("Error: status of coverage set {} unidentified".format(cov_set))


def main():
    """Call module directly with preprocessed config file"""

    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(config_path)

    process_coverage(cfg)


if __name__ == '__main__':
    main()
