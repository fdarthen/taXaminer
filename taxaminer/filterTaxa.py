#!/usr/bin/env python


"""setup external dependencies for taXaminer

"""

__author__ = "Freya Arthen"
__version__ = "0.6.0"

from . import checkInput
from . import parseGFF

import argparse
import os
import taxopy
import pandas as pd
import sys
import yaml
from Bio import Seq

class Config:
    def __init__(self, cfg_dict):
        self.gff_path = cfg_dict.get('gff_path')
        self.include_pseudogenes = cfg_dict.get('include_pseudogenes')

def init_gff_df(gff_path, filter_seq_ids):


    cfg_dict = {'gff_path': gff_path,
           'include_pseudogenes': True}
    cfg = Config(cfg_dict)

    gff_df = parseGFF.parse_genes(cfg)
    gff_df = parseGFF.strip_ids(gff_df)
    gff_df = parseGFF.add_missing_ids(gff_df)
    gff_df = parseGFF.remove_duplicate_ids(gff_df)
    gff_df = parseGFF.set_ids2index(gff_df)

    filtered_gff_df = gff_df.loc[(gff_df['type'] == 'gene') & (gff_df.index.isin(filter_seq_ids))]

    return filtered_gff_df

def df2gene_dict(gff_df):
    gff_df['pos'] = list(zip(zip(gff_df['start'],gff_df['end'])))
    dict= {}
    for scaffold in gff_df['scaffold'].unique().tolist():
        dict[scaffold] = gff_df.loc[gff_df['scaffold']==scaffold, ['pos', 'strand']].to_dict('index')


    return dict


def df2cds_dict(gff_df):
    dict = {}
    gff_df.rename(columns={'coding_coordinates': 'pos'}, inplace=True)
    for scaffold in gff_df['scaffold'].unique().tolist():
        dict[scaffold] = gff_df.loc[
            gff_df['scaffold'] == scaffold, ['pos', 'strand']].to_dict('index')


    return dict



def get_gene_coordinates(gene_table_path, gff_path, gene_ids):
    all_gene_loc_df = pd.read_csv(gene_table_path,
                                   usecols=["g_name", "c_name", "start","end"])
    gene_loc_df = all_gene_loc_df.loc[all_gene_loc_df['g_name'].isin(gene_ids)]
    gene_loc_df['pos'] = list(zip(gene_loc_df['start'], gene_loc_df['end']))

    gene_locs = {}
    for contig_id in gene_loc_df['c_name'].unique().tolist():
        contig_genes = gene_loc_df.loc[gene_loc_df['c_name']==contig_id]
        gene_locs[contig_id] = {row['g_name']: {'pos': (row['pos'],) } for index, row in contig_genes.iterrows()}

    with open(gff_path, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                spline = line.strip().split('\t')
                if spline[2] == 'gene':
                    gene_id = spline[-1].split('ID=')[1].split(';')[0]
                    if gene_id in gene_ids:
                        gene_locs[spline[0]][gene_id]['strand'] = spline[6]

    return gene_locs


def filter_genes(fasta_path, out_path, seq_locs):
    save_seq=False
    with open(fasta_path, 'r') as in_file, \
            open(out_path, 'w') as out_file:
        for line in in_file:
            if line.startswith('>'):
                if save_seq:
                    for subseq_id, info in seq_locs.get(header).items():
                        out_file.write(f">{subseq_id}\n")
                        subseq = ''
                        for coord in info.get('pos'):
                            subseq = subseq + seq[coord[0]-1:coord[1]]
                        if info.get('strand') == '-':
                            subseq = str(Seq.Seq(subseq).reverse_complement())
                        out_file.write(f"{subseq}\n")

                save_seq = False
                header = line[1:].split()[0]
                if header in seq_locs.keys():
                    # header is in list and header in list are to be written
                    save_seq=True
                    seq = ''
                else:
                    save_seq = False
            else:
                if save_seq:
                    seq = seq + line.strip()

        if save_seq:
            for subseq_id, info in seq_locs.get(header).items():
                out_file.write(f">{subseq_id}\n")
                subseq = ''
                for coord in info.get('pos'):
                    subseq = subseq + seq[coord[0] - 1:coord[1]]
                if info.get('strand') == '-':
                    subseq = str(Seq.Seq(subseq).reverse_complement())
                out_file.write(subseq)



def filter_fasta(fasta_path, out_path, keep_ids):
    write = False
    with open(fasta_path, 'r') as in_file, \
        open(out_path, 'w') as out_file:
        for line in in_file:
            if line.startswith('>'):
                header = line[1:].split()[0]
                if header in keep_ids:
                    # header is in list and header in list are to be written
                    write = True
                    out_file.write(line)
                else:
                    write = False
            else:
                if write:
                    out_file.write(line)


def has_assignment(TAX_DB, assignment_ids, filter_taxon):
    # find taxonomic assignments, that have filter_taxon in their lineage

    valid_ids = []

    for id in assignment_ids:
        try:
            taxon = taxopy.Taxon(id, TAX_DB)
        except:
            continue
        if filter_taxon.taxid in taxon.taxid_lineage:
            valid_ids.append(taxon.taxid)

    return valid_ids


def ident_matched_genes(TAX_DB, gene_table_path, filter_taxon):
    gene_assignments = pd.read_csv(gene_table_path,
                                   usecols=["g_name", "c_name", "fasta_header",
                                            "best_hit", "best_hitID",
                                            "bh_pident", "bh_evalue", "bh_bitscore",
                                            "lca", "lcaID",
                                            "refined_lca", "refined_lcaID",
                                            "taxon_assignment", "taxon_assignmentID",
                                            "is_target", "start","end",
                                            "upstream_gene", "downstream_gene"])

    selected_taxon_assignments = has_assignment(
        TAX_DB,
        gene_assignments['taxon_assignmentID'].unique(),
        filter_taxon)

    selected_genes = gene_assignments[
        gene_assignments['taxon_assignmentID'].isin(selected_taxon_assignments)]

    return selected_genes['g_name']



def ident_matched_contigs(TAX_DB, contig_table_path, filter_taxon):
    contig_assignments = pd.read_csv(contig_table_path)

    selected_taxon_assignments = has_assignment(
        TAX_DB,
        contig_assignments['lcaID'].unique(),
        filter_taxon)

    selected_contigs = contig_assignments[
        contig_assignments['lcaID'].isin(selected_taxon_assignments)]

    return selected_contigs['c_name']

def taxon_based_filtering(args, TAX_DB, gene_table_path,contig_table_path,
                          id_mapper, fasta_path, gff_path):
    filter_taxon = taxopy.Taxon(args.taxon, TAX_DB)

    print(f"Returning {args.out_seq_type} sequences with {args.assignment_level} "
          f"assignment to the taxon {filter_taxon.name} (ID: {filter_taxon.taxid})")
    if args.assignment_level == 'gene':
        keep_seq_ids = ident_matched_genes(TAX_DB, gene_table_path,
                                           filter_taxon)
    else:
        keep_seq_ids = ident_matched_contigs(TAX_DB, contig_table_path,
                                             filter_taxon)

    if args.out_seq_type == "contig":
        if args.assignment_level == 'contig':
            filter_seq_ids = keep_seq_ids
        else:
            filter_seq_ids = \
                id_mapper.loc[id_mapper['g_name'].isin(keep_seq_ids)][
                    'c_name'].unique()
        print(f"Returning {len(filter_seq_ids)} contigs.")
        filter_fasta(fasta_path, args.output, filter_seq_ids)
    elif args.out_seq_type == "gene":
        if args.assignment_level == 'contig':
            # return all gene sequences of considered contigs
            filter_seq_ids = \
                id_mapper.loc[id_mapper['c_name'].isin(keep_seq_ids)][
                    'g_name'].unique()
        else:
            filter_seq_ids = keep_seq_ids.values
        print(f"Returning {len(filter_seq_ids)} genes.")
        gff_df = init_gff_df(gff_path, filter_seq_ids)
        gene_locs = df2gene_dict(gff_df)
        filter_genes(fasta_path, args.output, gene_locs)
    elif args.out_seq_type == "cds":
        if args.assignment_level == 'contig':
            # return all gene sequences of considered contigs
            filter_seq_ids = \
                id_mapper.loc[id_mapper['c_name'].isin(keep_seq_ids)][
                    'g_name'].unique()
        else:
            filter_seq_ids = keep_seq_ids.values
        print(f"Returning CDS sequences of {len(filter_seq_ids)} genes.")
        gff_df = init_gff_df(gff_path, filter_seq_ids)
        cds_locs = df2cds_dict(gff_df)
        filter_genes(fasta_path, args.output, cds_locs)
    elif args.out_seq_type == "protein":
        if args.assignment_level == 'contig':
            filter_seq_ids = \
                id_mapper.loc[id_mapper['c_name'].isin(keep_seq_ids)][
                    'fasta_header'].unique()
        else:
            filter_seq_ids = \
                id_mapper.loc[id_mapper['g_name'].isin(keep_seq_ids)][
                    'fasta_header'].unique()
        print(f"Returning {len(filter_seq_ids)} proteins.")
        filter_fasta(fasta_path, args.output, filter_seq_ids)
    else:
        sys.exit(
            'Returning reads for a certain taxon is not supported yet.')


def seqid_based_filtering(args, seq_ids, id_mapper, gene_table_path, fasta_path, gff_path):
    print(f"Returning {args.out_seq_type} sequences with the specified sequence"
          f" identifiers")

    if args.in_seqid_type:
        if args.out_seq_type == "contig":
            if args.in_seqid_type == "contig":
                filter_seq_ids = seq_ids
            elif args.in_seqid_type == "gene":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['g_name'].isin(seq_ids)][
                        'c_name'].unique()
            elif args.in_seqid_type == "protein":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['fasta_header'].isin(seq_ids)][
                        'c_name'].unique()
            else:
                sys.exit("Read filtering not supported yet.")

            filter_fasta(fasta_path, args.output, filter_seq_ids)

        elif args.out_seq_type == "gene":
            if args.in_seqid_type == "contig":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['c_name'].isin(seq_ids)][
                        'g_name']
            elif args.in_seqid_type == "gene":
                filter_seq_ids = seq_ids
            elif args.in_seqid_type == "protein":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['fasta_header'].isin(seq_ids)][
                        'g_name'].unique()
            else:
                sys.exit("Read filtering not supported yet.")
            gff_df = init_gff_df(gff_path, filter_seq_ids)
            gene_locs = df2gene_dict(gff_df)
            filter_genes(fasta_path, args.output, gene_locs)

        elif args.out_seq_type == "cds":
            if args.in_seqid_type == "contig":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['c_name'].isin(seq_ids)][
                        'g_name']
            elif args.in_seqid_type == "gene":
                filter_seq_ids = seq_ids
            elif args.in_seqid_type == "protein":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['fasta_header'].isin(seq_ids)][
                        'g_name'].unique()
            else:
                sys.exit("Read filtering not supported yet.")
            gff_df = init_gff_df(gff_path, filter_seq_ids)
            cds_locs = df2cds_dict(gff_df)
            filter_genes(fasta_path, args.output, cds_locs)


        elif args.out_seq_type == "protein":
            if args.in_seqid_type == "contig":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['c_name'].isin(seq_ids)][
                        'fasta_header'].unique()
            elif args.in_seqid_type == "gene":
                filter_seq_ids = \
                    id_mapper.loc[id_mapper['g_name'].isin(seq_ids)][
                        'fasta_header'].unique()
                print(filter_seq_ids)
            elif args.in_seqid_type == "protein":
                filter_seq_ids = seq_ids
            else:
                sys.exit("Read filtering not supported yet.")
            filter_fasta(fasta_path, args.output, filter_seq_ids)

        else:
            sys.exit("Read filtering not supported yet.")


    else:
        # in type == out type
        filter_seq_ids = seq_ids
        if args.out_seq_type == "gene":
            gff_df = init_gff_df(gff_path, filter_seq_ids)
            gene_locs = df2gene_dict(gff_df)
            filter_genes(fasta_path, args.output, gene_locs)
        elif args.out_seq_type == "cds":
            gff_df = init_gff_df(gff_path, filter_seq_ids)
            cds_locs = df2cds_dict(gff_df)
            filter_genes(fasta_path, args.output, cds_locs)
        else: # contig and protein return whole sequence of fasta
            filter_fasta(fasta_path, args.output, filter_seq_ids)




def main():
    """  """
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-i", "--input",
                          help="Path to taXaminer report or config file. If "
                               "report folder provided, assembly FASTA and GFF "
                               "may be required.",
                          action="store")
    required.add_argument("-t", "--taxon",
                          help="NCBI Taxonomy ID to filter FASTA by.",
                          action="store", type=int)
    required.add_argument("-ids", "--seq-ids",
                          help="Sequence IDs to retain.",
                          action="store")
    optional.add_argument("--in-seqid-type",
                          help="Type of input sequence IDs. If not specified,"
                               "the in-seq-type is assumed to match type of out-seq-type.",
                          choices=['gene', 'contig', 'protein', 'read'])
    required.add_argument("-o", "--output",
                          help="Path to filtered output file.",
                          action="store")
    optional.add_argument("-s", "--out-seq-type",
                          help="Sequence type to return.",
                          choices=['gene', 'contig', 'protein', 'read', 'cds'],
                          default='contig')
    optional.add_argument('-l', '--assignment-level',
                          help="Assignment level to base the taxon filtering on.",
                          choices=['gene', 'contig'], default='gene')
    optional.add_argument('-f', '--fasta',
                          help="Path to FASTA file that should be filtered. "
                               "If protein sequences should be reported, FASTA "
                               "file in taXaminer output directory is used.",
                          action="store")
    optional.add_argument('-gff', '--annotation',
                          help="Path to GFF file of gene annotation. Is required "
                               "if gene or cds sequences should be reported",
                          action="store")

    args = parser.parse_args()

    gff_path = None
    if os.path.isfile(args.input):
        #TODO: add official run log to txmnr output
        input_type = "config"
        cfg = yaml.safe_load(open(args.input,'r'))
        gene_table_path = f"{cfg.get('output_path')}/taxonomic_assignment/gene_table_taxon_assignment.csv"
        contig_table_path = f"{cfg.get('output_path')}/taxonomic_assignment/contig_assignments.csv"
        if args.out_seq_type == "contig":
            fasta_path = cfg.get('fasta_path')
        elif args.out_seq_type == "protein":
            fasta_path = f"{cfg.get('output_path')}/proteins.faa"
        elif args.out_seq_type == "gene":
            fasta_path = cfg.get('fasta_path')
            gff_path = cfg.get('gff_path')
        elif args.out_seq_type == "cds":
            fasta_path = cfg.get('fasta_path')
            gff_path = cfg.get('gff_path')
        else:
            sys.exit('Returning reads for a certain taxon is not supported yet.')
        print(f"Using FASTA file: {fasta_path} (based on config file)")
    elif os.path.isdir(args.input):
        input_type = "dir"
        gene_table_path = f"{args.input}/taxonomic_assignment/gene_table_taxon_assignment.csv"
        contig_table_path = f"{args.input}/taxonomic_assignment/contig_assignments.csv"
        if args.out_seq_type == "contig":
            if args.fasta:
                fasta_path = args.fasta
                print(f"Using FASTA file: {fasta_path} (based on user input)")
            else:
                sys.exit("Please provide FASTA file of assembly")
        elif args.out_seq_type == "protein":
            fasta_path = f"{args.input}/proteins.faa"
            print(f"Using FASTA file: {fasta_path} (based on taXaminer report)")
        elif args.out_seq_type == "gene":
            fasta_path = args.fasta
            gff_path = args.annotation
        elif args.out_seq_type == "cds":
            fasta_path = args.fasta
            gff_path = args.annotation
        else:
            sys.exit(
                'Returning reads for a certain taxon is not supported yet.')
    else:
        sys.exit('Input does not exist.')

    SCRIPT_DIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    database_dir = eval(open(f"{SCRIPT_DIR}/pathconfig.txt",
                             'r').readline().strip()).get('data_path')

    TAX_DB = taxopy.TaxDb(nodes_dmp=database_dir + "/nodes.dmp",
                          names_dmp=database_dir + "/names.dmp",
                          merged_dmp=database_dir + "/merged.dmp",
                          keep_files=True)

    id_mapper = pd.read_csv(gene_table_path,
                            usecols=["g_name", "c_name", "fasta_header"])

    if args.taxon:
        taxon_based_filtering(args, TAX_DB, gene_table_path, contig_table_path,
                              id_mapper, fasta_path, gff_path)

    elif args.seq_ids:
        if os.path.exists(args.seq_ids):
            print("read sequence IDs from file")
            with open(args.seq_ids, 'r') as ids_file:
                seq_ids = [id.strip() for id in ids_file.readlines()]
        else:
            seq_ids = args.seq_ids.split(',')

        seqid_based_filtering(args, seq_ids, id_mapper, gene_table_path, fasta_path,
                              gff_path)

    else:
        sys.exit(
            'Please either specify taxon or sequence identifiers to use for'
            'filtering.')


if __name__ == "__main__":
    main()