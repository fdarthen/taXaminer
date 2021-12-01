# -*- coding: utf-8 -*-

import os
import sys
import yaml
import pathlib
import subprocess
import shutil
import logging
from logging.handlers import RotatingFileHandler

import prepare_and_check
import prepare_coverage
import produce_gene_info
import taxonomic_assignment
import modify_html
import extract_prot_seq


def main():


    SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

    user_config_dict = yaml.safe_load(open(sys.argv[1], 'r'))
    output_path = user_config_dict.get('output_path')
    log_level = user_config_dict.get('log_level')
    if not log_level:
        if len(sys.argv) > 2:
            if sys.argv[2] == '--quiet':
                log_level = logging.WARNING
            elif sys.argv[2] == '--verbose':
                log_level = logging.DEBUG
    if not log_level:
        log_level = logging.INFO

    logger = logging.basicConfig(level=log_level,
            handlers=[RotatingFileHandler(output_path+'milts.log', maxBytes=100000, backupCount=10)],
            format="[%(asctime)s] %(levelname)s [%(funcName)s:%(lineno)d] %(message)s",
            datefmt='%Y-%m-%dT%H:%M:%S')

    if not output_path.endswith('/'):
        output_path = output_path + '/'

    pathlib.Path(output_path).mkdir(parents=True, exist_ok=True)
    pathlib.Path(output_path+'tmp/').mkdir(parents=True, exist_ok=True)

    # set default for parameters not given by user
    # check file existence of required files and files to be computed in pipeline
    prepare_and_check.process_config(sys.argv[1], SCRIPT_DIR)

    # create class object with configuration parameters
    cfg = prepare_and_check.cfg2obj(output_path+'tmp/tmp.cfg.yml')

    if not cfg.update_plots:


        cmd_faidx_g = 'samtools faidx "{}" -o "{}tmp/tmp.MILTS.fasta.fai"'.format(
                            cfg.fasta_path, cfg.output_path)
        try:
            logging.info('creating genomic FASTA file index')
            out_faidx_g = subprocess.run([cmd_faidx_g], shell=True, capture_output=True, check=True)

        except:
            logging.error('Creation of genomic FASTA file failed')
            sys.exit('Error running\n{}'.format(cmd_faidx_g))


        prepare_and_check.make_checks(cfg)


        if cfg.compute_coverage:
            logging.info('Coverage information is calculated')
            prepare_coverage.process_coverage(cfg)

        logging.info('computing gene descriptors')
        produce_gene_info.process_gene_info(cfg)


        cmd_r_pca = 'Rscript {}/perform_PCA_and_clustering.R "{}"'.format(SCRIPT_DIR, cfg.cfg_path)
        try:
            out_r_pca = subprocess.run([cmd_r_pca], shell=True, capture_output=True, check=True)
        except:
            sys.exit('Error running\n{}'.format(cmd_r_pca))


        if cfg.extract_proteins:
            #TODO: update extraction to gff source rule
            cmd_faa = ' gffread -S --table "@id" -y {} -g {} {}'.format(
                            cfg.proteins_path+"gffread.faa", cfg.fasta_path, cfg.gff_path)
            try:
                out_faa = subprocess.run([cmd_faa], shell=True, capture_output=True, check=True)
            except:
                sys.exit('Error running\n{}'.format(cmd_faa))

            extract_prot_seq.generate_fasta(cfg)


    cmd_faidx_p = 'samtools faidx "{}" -o "{}tmp/tmp.proteins.fa.fai"'.format(
                        cfg.proteins_path, cfg.output_path)
    try:
        out_faidx_p = subprocess.run([cmd_faidx_p], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n{}'.format(cmd_faidx_p))


    taxonomic_assignment.run_assignment(cfg)


    cmd_r_plot = 'Rscript {}/plotting.R "{}"'.format(SCRIPT_DIR, cfg.cfg_path)
    try:
        out_r_plot = subprocess.run([cmd_r_plot], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Error running\n{}'.format(cmd_r_plot))


    modify_html.perform_adjustments(cfg)



    cmds_orca = []
    if cfg.output_pdf:
        cmds_orca.append('orca graph {}tmp/*.json -f "pdf" -d "{}taxonomic_assignment/"'.format(cfg.output_path, cfg.output_path))
    if cfg.output_png:
        cmds_orca.append('orca graph {}tmp/*.json -f "png" -d "{}taxonomic_assignment/"'.format(cfg.output_path, cfg.output_path))
    out_orca = ''

    #TODO: try fails because of regular errors
    for cmd_orca in cmds_orca:
        print(cmd_orca)
        try:
            out_orca = out_orca + subprocess.run([cmd_orca], shell=True, capture_output=True, check=True)
        except:
            sys.exit('Error running\n{}'.format(cmd_orca))

    try:
        shutil.rmtree("{}tmp/".format(cfg.output_path))
    except OSError as e:
        print("Error: %s : %s" % ("{}tmp/".format(cfg.output_path), e.strerror))


if __name__ == '__main__':
    main()
