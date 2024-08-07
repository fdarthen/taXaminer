"""Regroup the labels for the taxonomic assignments in PCA plot.

"""

__author__ = "Freya Arthen"

from . import compTaxonomicAssignment


import pandas as pd
import taxopy
import numpy as np
import sys
import argparse
import os
import yaml


def main():
    """."""

    parser = argparse.ArgumentParser(description=__doc__)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-i", "--input",
                          help="Path to taxaminer report directory or config file")
    optional.add_argument("-o", "--output",
                          help="path where output should be written if no overwrite",
                          default=None, type=str)
    optional.add_argument("-t", "--taxid",
                          help="NCBI taxonomy ID of target organism",
                          default=1, type=int)
    optional.add_argument("-m", "--merger",
                          help="works like merging_labels",
                          default=None)
    optional.add_argument("-n", "--num",
                          help="number of labels to show",
                          default=None, type=int)

    args = parser.parse_args()

    if os.path.isfile(args.input):
        cfg = yaml.safe_load(open(args.input, 'r'))
        gene_table_path = f"{cfg.get('output_path')}/taxonomic_assignment/gene_table_taxon_assignment.csv"
    elif os.path.isdir(args.input):
        gene_table_path = f"{args.input}/taxonomic_assignment/gene_table_taxon_assignment.csv"
    else:
        sys.exit('Input does not exist.')

    if args.num:
        num_groups_plot = args.num
        print(f"Number of labels limited to {num_groups_plot}")
    else:
        num_groups_plot = 25
        print(f"Defaulting number of labels to {num_groups_plot}")

    if args.output:
        output_path = args.output
    else:
        output_path = gene_table_path


    #TODO: check when correct taxonomy ID of target is required


    SCRIPT_DIR = os.path.abspath(
        os.path.dirname(os.path.realpath(__file__)))
    database_dir = eval(open(f"{SCRIPT_DIR}/pathconfig.txt",
                             'r').readline().strip()).get('data_path')

    TAX_DB = taxopy.TaxDb(nodes_dmp=database_dir + "/nodes.dmp",
                          names_dmp=database_dir + "/names.dmp",
                          merged_dmp=database_dir + "/merged.dmp",
                          keep_files=True)

    missing_ids = set()

    target_taxon = compTaxonomicAssignment.init_taxopyon(args.taxid,
                                                         missing_ids, TAX_DB)

    print(target_taxon)
    if not target_taxon:
        sys.exit("NCBI Taxonomy ID of target not recognized")

    # get gene table
    genes_df = pd.read_csv(gene_table_path, sep=',', index_col=None,
                           dtype={'start': int, 'end': int,
                                  'taxon_assignmentID': 'str', 'best_hitID': 'str',
                                  'refinedLCAID': 'str', 'lcaID': 'str'},
                           keep_default_na=False)

    genes_df = genes_df.replace(r'^\s*$', None, regex=True)

    genes_df[['plot_label','plot_labelID']] = None

    genes_df = compTaxonomicAssignment.merge_labels(genes_df, target_taxon, args.merger,
                                         missing_ids, TAX_DB)
    if num_groups_plot:
        num_labels = len(genes_df['plot_labelID'].unique())
        genes_df = compTaxonomicAssignment.iterative_merging_top(genes_df,
                                                      num_groups_plot - num_labels,
                                                      missing_ids, TAX_DB)
    genes_df = compTaxonomicAssignment.set_unassigned_labels(genes_df)

    print(genes_df['plot_label'].value_counts())

    #genes_df = genes_df.astype({'plot_labelID': 'Int64'})
    genes_df.to_csv(output_path, index=False)


if __name__ == '__main__':
    main()