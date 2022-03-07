#!/usr/bin/env python

"""taXaminer main script

Runs the pipeline by subsequently calling the submodules. Reads
path to raw config file, extracts and prepares information on this
basis and calls necessary processing steps.

Expects path to raw config file
"""

__author__ = "Freya Arthen"
__version__ = "0.6.0"

# package specific modules
import prepare_and_check
import prepare_coverage
import produce_gene_info
import taxonomic_assignment
import modify_html
import extract_prot_seq

import os
import sys
import yaml
import pathlib
import subprocess
import shutil
import logging
import time
import argparse



# custom logging format per level
class MyFormatter(logging.Formatter):
    """Object for logging."""
    def format(self, record):
        """Set different logging styles for different log levels."""
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        else:
            self._style._fmt = "[%(asctime)s] %(levelname)s [%(funcName)s:%(lineno)d] %(message)s"
        self.datefmt = "%Y-%m-%d %H:%M:%S"
        return super().format(record)

def main():
    """Run the program and call submodules."""

    # location of script on system to run from outside of directory
    SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

    start_time = time.time()

    # TODO: add more options to pass config parameters via command line
    # TODO: add levels of cleaning (temporarily) created data
    # process logging level information from command line input

    parser = argparse.ArgumentParser()
    parser.add_argument("config_path", help="path to config file")
    parser.add_argument("-q","--quiet", help="print only warnings and error messages",
                    action="store_true", default=False)
    parser.add_argument("-v","--verbose", help="print detailed log",
                    action="store_true", default=False)
    args = parser.parse_args()

    log_level = logging.INFO # default log level
    if args.quiet:
        log_level = logging.WARNING
    elif args.verbose:
        log_level = logging.DEBUG
    # initiate the logger
    logger = logging.getLogger()
    handler = logging.StreamHandler()
    handler.setFormatter(MyFormatter())
    logger.setLevel(log_level)
    logger.addHandler(handler)

    # read config file and retrieve output path
    if pathlib.Path(args.config_path).is_file():
        user_config_dict = yaml.safe_load(open(args.config_path, 'r'))
    else:
        logging.error('config file not found:   {}'.format(args.config_path))
        sys.exit()
    output_path = user_config_dict.get('output_path') # location for results
    if not output_path.endswith('/'): # sanity check path to dir
        output_path = output_path + '/'
    # create dictionary for output
    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path+'tmp/').mkdir(parents=True, exist_ok=True)

    ## create config file to be used by subsequent modules
    # includes default for parameters not given by user and checks
    # file existence of required files and files to be computed in pipeline
    prepare_and_check.process_config(args.config_path, SCRIPT_DIR)
    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(output_path+'tmp/tmp.cfg.yml')

    if not cfg.update_plots:

        pre_time = time.time()
        logging.info('>>> creating genomic FASTA index')
        cmd_faidx_g = 'samtools faidx "{}" -o "{}tmp/tmp.taxaminer.fasta.fai"'.format(
                            cfg.fasta_path, cfg.output_path)
        out_faidx_g = subprocess.run([cmd_faidx_g], shell=True, capture_output=True)
        if out_faidx_g.returncode != 0:
            logging.error('creation of genomic FASTA index failed:\n{}'.format(cmd_faidx_g))
            logging.error('Error message:\n'+out_faidx_g.stderr.decode())
            sys.exit()
        logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))

        ## make data integrity checks
        # check GFF and assembly FASTA for ID compatibility
        prepare_and_check.check_assembly_ids(cfg)

        if cfg.compute_coverage:
            pre_time = time.time()
            logging.info('>>> calculating coverage information')
            prepare_coverage.process_coverage(cfg)
            logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))

        if cfg.include_coverage:
            prepare_and_check.check_pbc_ids(cfg)

        pre_time = time.time()
        logging.info('>>> computing gene descriptors')
        produce_gene_info.process_gene_info(cfg)
        logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))

        pre_time = time.time()
        logging.info('>>> executing PCA and clustering ')
        cmd_r_pca = 'Rscript {}/perform_PCA_and_clustering.R "{}"'.format(SCRIPT_DIR, cfg.cfg_path)
        out_r_pca = subprocess.run([cmd_r_pca], shell=True, capture_output=True)
        if out_r_pca.returncode != 0:
            logging.error('running PCA failed:\n{}'.format(cmd_r_pca))
            if "gene_info/imputed_gene_table.txt': No such file or directory\n" in out_r_pca.stderr.decode():
                # file from previous step not available
                logging.error('Missing gene info file: "gene_info/imputed_gene_table.txt".\
                               Please check for previous errors and try restarting.')
            else:
                logging.error('PCA Error message:\n'+out_r_pca.stderr.decode())
            sys.exit()
        logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))

        if cfg.extract_proteins:
            pre_time = time.time()
            logging.info('>>> extracting protein sequences')
            extract_prot_seq.generate_fasta(cfg)
            logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))
        else:
            # create symlink if proteins.faa file does not exist in taxaminer
            # report (for GUI)
            if not pathlib.Path(cfg.output_path+'proteins.faa').exists():
                pathlib.Path(cfg.output_path+'proteins.faa').symlink_to(pathlib.Path(cfg.proteins_path).resolve())

    pre_time = time.time()
    logging.info('>>> creating protein FASTA index')
    cmd_faidx_p = 'samtools faidx "{}" -o "{}tmp/tmp.proteins.fa.fai"'.format(
                        cfg.proteins_path, cfg.output_path)
    out_faidx_p = subprocess.run([cmd_faidx_p], shell=True, capture_output=True)
    if out_faidx_p.returncode != 0:
        logging.error('creation of protein FASTA index failed:\n{}'.format(cmd_faidx_p))
        logging.error('Error message:\n'+out_faidx_p.stderr.decode())
        sys.exit()
    logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))

    pre_time = time.time()
    logging.info('>>> running taxonomic assignment')
    taxonomic_assignment.run_assignment(cfg)
    logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))

    pre_time = time.time()
    logging.info('>>> creating plots')
    cmd_r_plot = 'Rscript {}/plotting.R "{}"'.format(SCRIPT_DIR, cfg.cfg_path)
    out_r_plot = subprocess.run([cmd_r_plot], shell=True, capture_output=True)
    if out_r_plot.returncode != 0:
        logging.error('creation of plots failed:\n{}'.format(cmd_r_plot))
        if "taxonomic_assignment/gene_table_taxon_assignment.csv': No such file or directory\n" in out_r_plot.stderr.decode():
            # file from previous step not available
            logging.error('Missing file: "taxonomic_assignment/gene_table_taxon_assignment.csv".\
                           Please check for previous errors and try restarting.')
        else:
            logging.error('Error message:\n'+out_r_plot.stderr.decode())
        sys.exit()

    # make HTML file self-contained, text selectable and change title
    modify_html.perform_adjustments(cfg)

    # create static PDF versions of 3D plot
    # (only works on workstation with an graphical interface)
    cmds_orca = []
    if cfg.output_pdf:
        cmds_orca.append('orca graph {}tmp/*.json -f "pdf" -d "{}taxonomic_assignment/"'.format(cfg.output_path, cfg.output_path))
    if cfg.output_png:
        cmds_orca.append('orca graph {}tmp/*.json -f "png" -d "{}taxonomic_assignment/"'.format(cfg.output_path, cfg.output_path))
    out_orca = ''

    for cmd_orca in cmds_orca:
        logging.debug(cmd_orca)
        out_orca = subprocess.run([cmd_orca], shell=True, capture_output=True)
        out_orca_e = out_orca.stderr.decode()
        for line in out_orca_e.split('\n'):
            # parsing regular orca errors
            if 'error' in line.lower():
                if '[.DisplayCompositor]GL ERROR :GL_INVALID_OPERATION : glBufferData' in line:
                    pass
            elif 'Gtk-WARNING' in line and 'cannot open display' in line:
                logging.warning('Creation of static versions of 3D plot unavailable on cluster environment. \
                                Rerun with setting option "update_plots" to TRUE on a computer with a graphics unit to create plots.')
            elif 'done with code' in line:
                    pass
            else:
                print(line)
    logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))



    try:
        pre_time = time.time()
        logging.info('>>> deleting temporary files')
        shutil.rmtree("{}tmp/".format(cfg.output_path))
        logging.debug('finished [{}s]\n'.format(int(time.time()-pre_time)))
    except OSError as e:
        print("Error: %s : %s" % ("{}tmp/".format(cfg.output_path), e.strerror))


    logging.info('Total runtime: {}s\n'.format(int(time.time()-start_time)))


if __name__ == '__main__':
    main()
