# -*- coding: utf-8 -*-

import os
import sys
import yaml
import pathlib
import subprocess
import shutil
import logging

import prepare_and_check
import prepare_coverage
import produce_gene_info
import taxonomic_assignment
import modify_html
import extract_prot_seq

# custom logging format per level
class MyFormatter(logging.Formatter):
    def format(self, record):
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        else:
            self._style._fmt = "[%(asctime)s] %(levelname)s [%(funcName)s:%(lineno)d] %(message)s"
        self.datefmt = "%Y-%m-%d %H:%M:%S"
        return super().format(record)

def main():

    SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

    log_level = logging.INFO
    if len(sys.argv) > 2:
        if sys.argv[2] == '--quiet':
            log_level = logging.WARNING
        elif sys.argv[2] == '--verbose':
            log_level = logging.DEBUG

    # initiate the logger
    logger = logging.getLogger()
    handler = logging.StreamHandler()
    handler.setFormatter(MyFormatter())
    logger.setLevel(log_level)
    logger.addHandler(handler)

    try:
        user_config_dict = yaml.safe_load(open(sys.argv[1], 'r'))
    except:
        logging.error('config file not found:   {}'.format(sys.argv[1]))
        sys.exit()

    output_path = user_config_dict.get('output_path')

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
        logging.info('>>> creating genomic FASTA file index')
        try:
            out_faidx_g = subprocess.run([cmd_faidx_g], shell=True, capture_output=True, check=True)
        except:
            logging.error('creation of genomic FASTA file failed')
            sys.exit('Error running\n{}'.format(cmd_faidx_g))


        prepare_and_check.make_checks(cfg)


        if cfg.compute_coverage:
            logging.info('>>> calculating coverage information')
            prepare_coverage.process_coverage(cfg)

        logging.info('>>> computing gene descriptors')
        try:
            produce_gene_info.process_gene_info(cfg)
        except:
            logging.error('computing gene descriptors failed')
            sys.exit()

        logging.info('>>> executing PCA and clustering ')
        cmd_r_pca = 'Rscript {}/perform_PCA_and_clustering.R "{}"'.format(SCRIPT_DIR, cfg.cfg_path)
        try:
            out_r_pca = subprocess.run([cmd_r_pca], shell=True, capture_output=True, check=True)
        except:
            sys.exit('Error running\n{}'.format(cmd_r_pca))


        if cfg.extract_proteins:
            #TODO: update extraction to gff source rule
            cmd_faa = ' gffread -S --table "@id" -y {} -g {} {}'.format(
                            cfg.proteins_path+"gffread.faa", cfg.fasta_path, cfg.gff_path)
            logging.info('>>> extracting protein sequences')
            try:
                out_faa = subprocess.run([cmd_faa], shell=True, capture_output=True, check=True)
            except:
                logging.error('extraction of protein sequences into FASTA file failed')
                sys.exit('Error running\n{}'.format(cmd_faa))

            extract_prot_seq.generate_fasta(cfg)


    cmd_faidx_p = 'samtools faidx "{}" -o "{}tmp/tmp.proteins.fa.fai"'.format(
                        cfg.proteins_path, cfg.output_path)
    logging.info('>>> creating protein FASTA file index')
    try:
        out_faidx_p = subprocess.run([cmd_faidx_p], shell=True, capture_output=True, check=True)
    except:
        logging.error('creation of protein FASTA file failed')
        sys.exit('Error running\n{}'.format(cmd_faidx_p))

    logging.info('>>> running taxonomic assignment')
    taxonomic_assignment.run_assignment(cfg)


    cmd_r_plot = 'Rscript {}/plotting.R "{}"'.format(SCRIPT_DIR, cfg.cfg_path)
    logging.info('>>> creating plots')
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

    #FIXME: try fails because of regular errors
    for cmd_orca in cmds_orca:
        print(cmd_orca)
        out_orca = subprocess.run([cmd_orca], shell=True, capture_output=True)
        out_orca_e = out_orca.stderr.decode()
        for line in out_orca_e.split('\n'):
            if 'error' in line.lower():
                if '[.DisplayCompositor]GL ERROR :GL_INVALID_OPERATION : glBufferData' in line:
                    pass
                elif 'Gtk-WARNING' in line and 'cannot open display' in line:
                    logging.warning('Creation of static versions of 3D plot unavailable on cluster environment. \
                                    Rerun with setting option "update_plots" to TRUE on a workstation.')
                elif 'done with code' in line:
                    pass
            else:
                print(line)


    # delete temporary files
    try:
        logging.info('>>> deleting temporary files')
        shutil.rmtree("{}tmp/".format(cfg.output_path))
    except OSError as e:
        print("Error: %s : %s" % ("{}tmp/".format(cfg.output_path), e.strerror))


if __name__ == '__main__':
    main()
