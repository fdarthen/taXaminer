#!/usr/bin/env python

"""Add taxonomic assignment of genes to GFF

Adds the taxonomic assignment of each genes as attribute to GFF

Expects prepared config file
(preparation of config by prepare_and_check.py)
"""
__author__ = "Freya Arthen"
__version__ = "0.6.0"

import logging
import sys
import pathlib
import subprocess
import csv
import argparse

# package specific modules
import prepare_and_check

def remove_prefix(text, prefix):
    """

    Args:
      text:
      prefix:

    Returns:

    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def strip_ID(id):
    """
    Remove GFFs prefixes from IDs

    Args:
      id:

    Returns:

    """
    for prefix in ['gene:','gene-','transcript:','transcript-','rna:','rna-','cds:','cds-']:
        id = remove_prefix(id,prefix)
    return id

def write2gff(gff_in_path, gff_out_path, genes_taxon_dict):
    """

    Args:
      cfg: Config object with config parameters
    """

    with open(gff_in_path, "r") as gff_in: # open GFF file
        # open output file to write to
        with open(gff_out_path, "w") as gff_out:

            # read every line of the original GFF file
            for line in gff_in:
                # remove newline from line
                line = line.strip()

                # if anntations have been reached:
                if not line.startswith("#"):
                    # split GFF line to access the individual column entries
                    gff_array = line.split('\t')
                    # get gene name from attributes entry in 9th column
                    id = strip_ID(gff_array[8].split('ID=')[1].split(';')[0])

                    # if there is taXaminer info on this gene
                    if id in genes_taxon_dict.keys():
                        # add a taXaminer attribute to the list of attributes in the 9th column of the GFF file
                        line = '{};assigned_taxon={}'.format(line,genes_taxon_dict.get(id))
                # whether the line has been modified or not, write it to the output GFF
                gff_out.write(line + "\n")


def read_assignments(taxaminer_csv_path):

    genes_taxon_dict = {}

    # read taXaminer file to retrieve for all genes taXaminer information is available on
    with open(taxaminer_csv_path, "r") as taxaminer:
        taxaminer_dict = csv.DictReader(taxaminer, delimiter=',')
        for row intaxaminer_dict:
            genes_taxon_dict[row.get('g_name')] = row.get('taxon_assignment')

    return genes_taxon_dict


def add2gff(taxaminer_csv_path, gff_in, gff_out):

    genes_taxon_dict = read_assignments(taxaminer_csv_path)
    write2gff(gff_in, gff_out, genes_taxon_dict)


def main():
    """Call module directly with preprocessed config file"""

    parser = argparse.ArgumentParser(description='Add taxonomic assignment to GFF')

    parser.add_argument('--config-file', '-c', type=str, action='store',
                       help='path to config file')
    parser.add_argument('--gff-in', '-g', action='store',
                        help='path of GFF input file')
    parser.add_argument('--gff-out', '-o', action='store',
                        help='path to save modified GFF to')
    parser.add_argument('--taxaminer', '-m', action='store',
                        help='path to either taXaminer report or taXaminer taxonomic assignment tabluar output')
    args = parser.parse_args()

    if args.config_file:
        config_path = sys.argv[1]
        # create class object with configuration parameters
        cfg = prepare_and_check.cfg2obj(args.config_file)
        add2gff(cfg.output_path+'taxonomic_assignment/gene_table_taxon_assignment.csv',
                cfg.gff_path, cfg.gff_ta_path)
    else:
        if args.taxaminer.endswith('.csv'):
            add2gff(args.taxaminer,args.gff_in,args.gff_out)
        else:
            add2gff(args.taxaminer+'taxonomic_assignment/gene_table_taxon_assignment.csv',args.gff_in,args.gff_out)


if __name__ == '__main__':
    main()
