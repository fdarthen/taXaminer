#!/usr/bin/env python

"""Prepare data for running taXaminer

Parses the GFF file

Prepares coverage data for further processing. Raw reads are mappend
and mapping files (BAM) are converted to per base coverage (PBC) files

Compute protein FASTA by parsing gene sequences baed on coordinates from
GFF and FASTA with scaffold sequences (assembly FASTA) and converting
it to amino acid sequence

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
import time
import timeit

from Bio import Seq
import pandas as pd

from . import checkInput

def assess_parent_path(features, gene, transcript, exon, cds):
    transcript_type = 'gene'  # feature type of which to use in the protein FASTA
    if transcript:
        transcript_avail = True
        transcript_type = 'mRNA'
    if exon:
        exon_avail = True
        if transcript_type == 'gene':
            transcript_type = 'exon'
    if cds:
        cds_avail = True
        coding_feature = cds
        if transcript_type == 'gene' or transcript_type == 'exon':
            transcript_type = 'CDS'
    else:
        if exon:
            coding_feature = exon
        else:
            if transcript:
                coding_feature = transcript
            else:
                coding_feature = gene

    coding_type = coding_feature.get('type')
    parent_path = [coding_type]
    parent = coding_feature.get('id')
    while not features.get(parent).get('type') == 'gene':
        parent = features.get(parent).get('parent_id')
        parent_path.append(features.get(parent).get('type'))

    return coding_type, transcript_type, parent_path

def get_trancsript_type(cfg, features, gene):
    # if proteins will not be extracted by taxaminer
    # check of which feature type the ID is used in the protein FASTA
    transcript_type = None
    with open(cfg.proteins_path, 'r') as p_file:
        for line in p_file:
            if line.startswith('>'):
                protein_id = line.strip().split('|')[0][1:]
                # assess the feature whose ID is given in the first position
                # of the header in the fasta file
                for id, feat in features.items():
                    if id == protein_id:
                        transcript_type = feat.get('type')
                        break
                    elif feat.get('add_attrs').get('ID') == protein_id:
                        transcript_type = feat.get('type')
                        break
            else:
                if transcript_type:
                    break
                continue

    return transcript_type


def check_gff(cfg):

    gene, transcript, exon, cds = None, None, None, None
    features = {}
    var_list = []
    with open(cfg.gff_path, 'r') as gff_file:
        for line in gff_file:
            # skip comments
            if line.startswith('#'):
                continue
            spline = line.strip().split('\t')
            feature_dict = spline2dict(spline, cfg.include_pseudogenes)
            if feature_dict.get('type') == 'gene':
                # read only first gene
                if gene:
                    coding_type, transcript_type, parent_path = assess_parent_path(
                        features, gene, transcript, exon, cds)
                    if not cfg.extract_proteins:
                        transcript_type = get_trancsript_type(
                            cfg, features, gene)
                    if var_list:
                        if (var_list == [coding_type, transcript_type, parent_path])\
                                and (not cfg.extract_proteins and transcript_type):
                            break
                    var_list = [coding_type, transcript_type, parent_path]
                else:
                    transcript, exon, cds = None, None, None
                    features = {}
                gene = feature_dict
                features[gene.get('id')] = gene
            elif feature_dict.get('type') in ['tRNA', 'rRNA', 'lnc_RNA', 'lncRNA', 'ncRNA']:
                gene, transcript, exon, cds = None, None, None, None
                features = {}
            elif feature_dict.get('type') == 'mRNA':
                transcript = feature_dict
                features[transcript.get('id')] = transcript
            elif feature_dict.get('type') == 'exon':
                exon = feature_dict
                features[exon.get('id')] = exon
            elif feature_dict.get('type') == 'CDS':
                cds = feature_dict
                features[cds.get('id')] = cds
            else:
                continue

    return coding_type, transcript_type, parent_path


def match_fasta2id(cfg, gff_df):
    """ Match the header of a precomputed protein FASTA file to gene ID

    :param geneid_in_fasta:
    :return:
    """

    header_dict = {}
    if not cfg.extract_proteins:
        # if proteins will not be extracted by taxaminer
        # check of which feature type the ID is used in the protein FASTA
        with open(cfg.proteins_path, 'r') as p_file:
            for line in p_file:
                if line.startswith('>'):
                    protein_id = line.strip().split('|')[0][1:]
                    # diamond uses everything until first whitespace as ID
                    header_dict[protein_id] = line[1:].strip().split()[0]
                else:
                    continue
        # if protein was passed by user, 'transcript_id' resembles the id
        # of the feature that denotes the header in the fasta file
        gff_df['fasta_header'] = gff_df['transcript_id'].map(header_dict)
        missing_header = gff_df.loc[gff_df['type'] == 'gene', 'fasta_header'].isna()
        if missing_header.any(): # and geneid_in_fasta:
            if gff_df.loc[gff_df['type'] == 'gene'][missing_header].shape[0] == gff_df.loc[gff_df['type'] == 'gene'].shape[0]:
                # if not any gene was mapped to a fasta header
                logging.error(
                    f"Unable to map the gene IDs to the header in the protein "
                    f"FASTA file. Header must start with ID of either the gene,"
                    f"the transcript or the coding sequence in the GFF.")
                sys.exit(1)
            else:
                if gff_df.loc[gff_df['type'] == 'gene'][~missing_header].shape[0] == len(header_dict):
                    # if every fasta file header was mapped
                    logging.info(f"All headers of the fasta file was mapped to "
                                 f"a gene ID. Additional {gff_df.loc[gff_df['type'] == 'gene'][missing_header].shape[0]} "
                                 f"genes remain without a protein sequence.")
                else:
                    logging.info(
                        f"Unable to map the following "
                        f"{gff_df.loc[gff_df['type'] == 'gene'][missing_header].shape[0]} "
                        f"gene IDs to the headers in the protein FASTA file:\n"
                        f"{gff_df.loc[gff_df['type'] == 'gene'][missing_header].index.to_list()}")
    else:
        gff_df['fasta_header'] = gff_df['transcript_id']


def remove_prefix(text, prefix):
    """
    Removes given prefix from text, returns stripped text.

    Args:
      text(str): text to remove prefix from
      prefix(str): prefix to remove from text

    Returns:
      (str) stripped text
    """
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def strip_ID(id):
    """
    Removes attribute-linked prefixes from IDs, returns stripped ID.

    Args:
      id(str): ID with prefix specifying type of attribute

    Returns:
      (str) stripped ID

    """
    if not id:
        return None

    for prefix in ['gene:', 'gene-', 'transcript:', 'transcript-',
                   'rna:', 'rna-', 'cds:', 'cds-']:
        id = remove_prefix(id, prefix)
    return id


def attrs2dict(attr_list):
    """ Parse the GFF additional attributes list into a dict

    :param attr_list:
    :return:
    """

    return {attr.split('=')[0]: attr.split('=')[1] for attr in
            attr_list.split(';')}


def spline2dict(spline, include_pseudogenes):

    spline_dict = {}
    # direct inference
    spline_dict['scaffold'] = spline[0]
    spline_dict['source'] = spline[1]
    if include_pseudogenes and spline[2] == 'pseudogene':
        spline_dict['type'] = 'gene'
    else:
        spline_dict['type'] = spline[2]
    spline_dict['start'] = int(spline[3])
    spline_dict['end'] = int(spline[4])
    spline_dict['score'] = spline[5]
    spline_dict['strand'] = spline[6]
    spline_dict['phase'] = spline[7]
    spline_dict['add_attrs'] = attrs2dict(spline[8])
    spline_dict['add_attrs']['ID'] = strip_ID(spline_dict['add_attrs'].get('ID'))
    spline_dict['add_attrs']['Parent'] = strip_ID(spline_dict['add_attrs'].get('Parent'))
    spline_dict['id'] = spline_dict['add_attrs'].get('ID')
    spline_dict['parent_id'] = strip_ID(spline_dict['add_attrs'].get('Parent'))
    transl_table = spline_dict['add_attrs'].get('transl_table')
    spline_dict['transl_table'] = transl_table if transl_table else '1'

    # indirect inference
    spline_dict['gene_id'] = spline_dict['id'] if spline_dict['type'] == 'gene' else None
    spline_dict['length'] = spline_dict['end'] - spline_dict['start'] + 1
    if not spline_dict.get('id'):
        # no ID information only acceptable for features without subfeatures
        spline_dict['id'] = f"{spline_dict.get('parent_id')}-{spline_dict['type']}-{spline_dict['start']}"

    return spline_dict


def replace_longest_transcript(gene, transcript, transcript_type, transcript_cds_features):

    cds_ids = [cds.get('id') for cds in transcript_cds_features]
    # check if CDS features have unique identifiers
    # else enumerate them
    if len(set(cds_ids)) != len(cds_ids):
        cds_ids = []
        for i, cds in enumerate(transcript_cds_features):
            cds['id'] = f"{cds.get('id')}-{i}"
            cds_ids.append(cds.get('id'))
    cds_coordindates = [(cds.get('start'), cds.get('end')) for cds in
                        transcript_cds_features]
    transcript['coding_features'] = cds_ids
    gene['coding_features'] = cds_ids
    transcript['coding_coordindates'] = cds_coordindates
    gene['coding_coordindates'] = cds_coordindates
    if transcript_type:
        if transcript_type == transcript.get('type'):
            gene['transcript_id'] = transcript.get('id')
        else:
            gene['transcript_id'] = transcript_cds_features[0].get('add_attrs').get('ID')
    else:
        gene['transcript_id'] = transcript.get('id')
    # use translational table of CDS feature is gene does not have info
    if gene.get('transl_table') == '1' and transcript_cds_features[0].get('transl_table') != '1':
        # if coding feature has different transl table info than default use this
        gene['transl_table'] = transcript_cds_features[0].get('transl_table')
    gene_cds_features = [transcript]
    gene_cds_features += transcript_cds_features

    return gene_cds_features


def save_gene_features(pandas_row_list, transcript, gene,
                       transcript_cds_features, transcript_cds_length,
                       gene_cds_features, max_cds_len, transcript_type):
    # process previous gene
    if transcript:
        if transcript_cds_length >= max_cds_len:
            gene_cds_features = replace_longest_transcript(
                gene, transcript, transcript_type, transcript_cds_features)
    if gene:
        if not transcript:
            cds_ids = [cds.get('id') for cds in transcript_cds_features]
            # check if CDS features have unique identifiers
            # else enumerate them
            if len(set(cds_ids)) != len(cds_ids):
                cds_ids = []
                for i, cds in enumerate(transcript_cds_features):
                    cds['id'] = f"{cds.get('id')}-{i}"
                    cds_ids.append(cds.get('id'))
            cds_coordindates = [(cds.get('start'), cds.get('end')) for cds in
                                transcript_cds_features]
            gene['coding_features'] = cds_ids
            gene['coding_coordindates'] = cds_coordindates
            gene['transcript_id'] = gene.get('id')
            gene_cds_features = transcript_cds_features
        if gene.get('coding_features'):
            pandas_row_list.append(gene)
            pandas_row_list += gene_cds_features
        else:
            # gene has no subfeatures -> gene = transcript = CDS
            gene['transcript_id'] = gene.get('id')
            gene['coding_features'] = [gene.get('id')]
            pandas_row_list.append(gene)



def parse_file(cfg, coding_type, transcript_type, parent_path):
    """

    :param coding_type:
    :param parent_path:
    :param cfg:
    :return:
    """

    pandas_row_list = []

    gene, transcript, exon, cds = None, None, None, None
    transcript_cds_features, transcript_cds_length, gene_cds_features = None, None, None
    max_cds_len = None

    with open(cfg.gff_path, 'r') as gff_file:
        for line in gff_file:

            # skip comments
            if line.startswith('#'):
                # stop parsing when FASTA block is reached
                if line.startswith('##FASTA'):
                    break
                else:
                    continue

            spline = line.strip().split('\t')

            # only relevant lines will be further processed
            ## type in parent_path or a pseudogene
            if (not spline[2] in parent_path) and (spline[2] != 'pseudogene'):
                if spline[2] in ['tRNA', 'rRNA', 'lnc_RNA', 'lncRNA', 'ncRNA']:
                    # only process protein coding genes
                    if gene:
                        if (gene['id'] == feature_dict['gene_id']):
                            # last parsed feature was assigned to active gene
                            # i.e.: there was a coding feature assigned to the gene
                            # catches: saving genes of coding features that are
                            # followed by non-coding RNAs without a gene feature
                            save_gene_features(pandas_row_list, transcript,
                                               gene, transcript_cds_features,
                                               transcript_cds_length,
                                               gene_cds_features, max_cds_len,
                                               transcript_type)
                    gene = None
                continue
            ## no further processing if pseudogenes are not ought to be included
            elif not cfg.include_pseudogenes and spline[2] == 'pseudogene':
                if gene:
                    # process previous gene
                    save_gene_features(pandas_row_list, transcript, gene,
                                       transcript_cds_features,
                                       transcript_cds_length, gene_cds_features,
                                       max_cds_len, transcript_type)
                gene = None
                continue

            # pseudogenes are returned with type "gene" by this function
            feature_dict = spline2dict(spline, cfg.include_pseudogenes)

            # feature == gene
            if (feature_dict.get('type') == parent_path[-1]):
                if gene:
                    # process previous gene
                    save_gene_features(pandas_row_list, transcript, gene,
                                       transcript_cds_features,
                                       transcript_cds_length, gene_cds_features,
                                       max_cds_len, transcript_type)

                # init new gene
                gene = feature_dict
                max_cds_len = 0
                gene_cds_features = []
                transcript, cds = None, None
                transcript_cds_features = []
                transcript_cds_length = 0

            # feature == coding feature
            elif feature_dict.get('type') == parent_path[0]:
                if not gene:
                    continue
                cds = feature_dict
                cds['gene_id'] = gene.get('id')
                """if cds.get('gene_id') == cds.get('parent_id'):
                    # no transcript features
                    if cds.get('length') > max_cds_len:
                        max_cds_len = cds.get('length')
                        gene_cds_features = [cds]
                        # update gene with info from cds
                        if transcript_type:
                            if transcript_type == gene.get('type'):
                                gene['transcript_id'] = gene.get('id')
                            else:
                                gene['transcript_id'] = cds.get('id')
                        else:
                            gene['transcript_id'] = gene.get('id')
                        gene['coding_features'] = [cds.get('id')]
                        gene['coding_coordindates'] = [(cds.get('start'), cds.get('end'))]
                        # use translational table of CDS feature is gene has default info
                        if gene.get('transl_table') == '1' and cds.get('transl_table') != '1':
                            gene['transl_table'] = cds.get('transl_table')
                else:"""
                if feature_dict.get('type') == 'mRNA':
                    # mRNA features, i.e. transcripts, are the coding features
                    # in the GFF; choose the longest mRNA per gene
                    if transcript:
                        if transcript.get('length') >= feature_dict.get('length'):
                            continue
                    transcript = feature_dict
                    transcript['gene_id'] = gene.get('id')
                    transcript['coding_features'] = [transcript.get('id')]
                    gene['coding_features'] = transcript.get('id')
                    transcript['coding_coordindates'] = [(transcript.get('start'), transcript.get('end'))]
                    gene['coding_coordindates'] = [(transcript.get('start'), transcript.get('end'))]
                else:
                    # add to list of CDS for transcript
                    transcript_cds_features.append(cds)
                    transcript_cds_length += cds.get('length')

             # feature between gene and CDS
            else:
                if not gene:
                    continue
                if feature_dict.get('type') == 'mRNA':
                    # process previous transcript
                    if transcript:
                        # update the rows to be added to the gff data frame
                        # if new transcript is longer;
                        # also update gene's transcript id
                        if transcript_cds_length >= max_cds_len:
                            max_cds_len = transcript_cds_length
                            gene_cds_features = replace_longest_transcript(
                                gene, transcript, transcript_type,
                                transcript_cds_features)
                    # init new transcript
                    transcript = feature_dict
                    transcript['gene_id'] = gene.get('id')
                    transcript_cds_features = []
                    transcript_cds_length = 0
                else:
                    # feature only required for matching to gene ID, which
                    # is already fulfilled by order of the GFF file
                    continue

    if gene:
        save_gene_features(pandas_row_list, transcript, gene,
                           transcript_cds_features, transcript_cds_length,
                           gene_cds_features, max_cds_len, transcript_type)

    gff_df = pd.DataFrame(pandas_row_list)
    # if features don't have ID attribute, concatenate type and line number
    no_id = gff_df['id'].isnull()
    gff_df.loc[no_id, 'id'] = gff_df['type'] + '-' + gff_df.index.astype(str)
    # check for duplicate gene IDs
    if gff_df.loc[gff_df['type'] == 'gene'].duplicated(subset='id', keep=False).any():
        logging.error('Duplicate gene IDs detected. Please correct your input GFF.')
        sys.exit(1)
    # append line number to duplicate ids
    duplicate_index = gff_df.duplicated(subset='id', keep=False)
    gff_df.loc[(duplicate_index) & (gff_df['type'] != 'gene'), 'id'] = gff_df['id'] + '-' + gff_df.index.astype(str)

    gff_df.set_index('id', inplace=True)

    return gff_df

def set_gene_neighbours(gff_df):

    cur_gene, cur_contig = None, None
    upstream_gene = None
    for row in gff_df.loc[gff_df['type'] == 'gene'].itertuples():
        if cur_contig == row.scaffold:
            gff_df.at[cur_gene, 'upstream_gene'] = upstream_gene
            gff_df.at[cur_gene, 'downstream_gene'] = row.Index
            upstream_gene = cur_gene
            cur_gene = row.Index
        else:
            if cur_gene:
                gff_df.at[cur_gene, 'upstream_gene'] = upstream_gene
                gff_df.at[cur_gene, 'downstream_gene'] = None

            upstream_gene = None
            cur_contig = row.scaffold
            cur_gene = row.Index


def process(cfg):
    coding_type, transcript_type, parent_path = check_gff(cfg)
    gff_df = parse_file(cfg, coding_type, transcript_type, parent_path)
    set_gene_neighbours(gff_df)
    match_fasta2id(cfg, gff_df)

    return gff_df


def main():
    """Call module directly with preprocessed config file"""
    config_path = sys.argv[1]
    # create class object with configuration parameters
    cfg = checkInput.cfg2obj(config_path)

    process(cfg)


if __name__ == '__main__':
    main()
