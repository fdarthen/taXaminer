# -*- coding: utf-8 -*-

import yaml # read config file
import taxopy
import time
import sys
import numpy as np
import csv
import os
import pathlib
import subprocess


script_dir = sys.argv[2]

TAX_DB = taxopy.TaxDb(nodes_dmp=script_dir+"/nodes.dmp",
                        names_dmp=script_dir+"/names.dmp", keep_files=True)


missing_taxids = set()

class Gene:
    def __init__(self, line_array, header_index):
        self.g_name = line_array[header_index.get('g_name')]
        self.c_name = line_array[header_index.get('c_name')]
        self.c_num_of_genes = line_array[header_index.get('c_num_of_genes')]
        self.c_len = line_array[header_index.get('c_len')]
        self.c_pct_assemby_len = line_array[header_index.get('c_pct_assemby_len')]
        self.c_genelenm = line_array[header_index.get('c_genelenm')]
        self.c_genelensd = line_array[header_index.get('c_genelensd')]

        self.c_cov = tuple([line_array[index] for index in
                            header_index.get('c_cov')])
        self.c_covsd = tuple([line_array[index] for index in
                            header_index.get('c_covsd')])
        self.c_covdev = tuple([line_array[index] for index in
                            header_index.get('c_covdev')])
        self.c_genecovm = tuple([line_array[index] for index in
                            header_index.get('c_genecovm')])
        self.c_genecovsd = tuple([line_array[index] for index in
                            header_index.get('c_genecovsd')])

        self.c_pearson_r = line_array[header_index.get('c_pearson_r')]
        self.c_pearson_p = line_array[header_index.get('c_pearson_p')]
        self.c_gc_cont = line_array[header_index.get('c_gc_cont')]
        self.c_gcdev = line_array[header_index.get('c_gcdev')]

        self.g_len = line_array[header_index.get('g_len')]
        self.g_lendev_c = line_array[header_index.get('g_lendev_c')]
        self.g_lendev_o = line_array[header_index.get('g_lendev_o')]
        self.g_abspos = line_array[header_index.get('g_abspos')]
        self.g_terminal = line_array[header_index.get('g_terminal')]
        self.g_single = line_array[header_index.get('g_single')]

        self.g_cov = tuple([line_array[index] for index in
                            header_index.get('g_cov')])
        self.g_covsd = tuple([line_array[index] for index in
                            header_index.get('g_covsd')])
        self.g_covdev_c = tuple([line_array[index] for index in
                            header_index.get('g_covdev_c')])
        self.g_covdev_o = tuple([line_array[index] for index in
                            header_index.get('g_covdev_o')])

        self.g_pearson_r_o = line_array[header_index.get('g_pearson_r_o')]
        self.g_pearson_p_o = line_array[header_index.get('g_pearson_p_o')]
        self.g_pearson_r_c = line_array[header_index.get('g_pearson_r_c')]
        self.g_pearson_p_c = line_array[header_index.get('g_pearson_p_c')]
        self.g_gc_cont = line_array[header_index.get('g_gc_cont')]
        self.g_gcdev_c = line_array[header_index.get('g_gcdev_c')]
        self.g_gcdev_o = line_array[header_index.get('g_gcdev_o')]

        self.coords = tuple([line_array[i] for i in header_index.get('Dim.')])

        self.protID = None # ID of transcript with longest CDS
        self.lencds = 0 # lenght of self.protID

        self.tax_assignments = [] # all taxonomic hits for the gene
        self.best_hit = None
        self.best_hitID = None
        self.bh_evalue = None
        self.lca = None
        self.lcaID = None
        self.corrected_lca = None
        self.corrected_lcaID = None
        self.taxon_assignment = None
        self.taxon_assignmentID = None

        self.plot_label = None
        self.plot_labelID = None

###############################################################################
############################ HELPER FUNCTIONS #################################
###############################################################################

def remove_prefix(text, prefix):
    """Removes given prefix from text -> return text w/o prefix."""
    if text.startswith(prefix):
        return text[len(prefix):]
    return text

def strip_ID(id):
    """Removes GFFs prefixes from ID -> return striped ID."""
    for prefix in ['gene:','gene-','transcript:','transcript-',
                    'rna:','rna-','cds:','cds-']:
        id = remove_prefix(id,prefix)
    return id


def get_id_for_rank_of_species(queryID, rank):
    """Returns ID of taxon at the given rank for given taxon ID.

    Args:
      queryID (int): NCBI ID of taxon for wich to find ID at given rank
      rank (str): rank of which to get the NCBI ID of given taxon with queryID

    Returns:
      int: NCBI ID of taxon at the given rank for given taxon ID
    """
    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_ranks = query_taxon.rank_name_dictionary
    rank_of_interest = query_ranks.get(rank)
    # get ID for the specific class (function taxopy.taxid_from_name() might get multiple IDs because of homonyms)
    if rank_of_interest:
        rank_index = query_taxon.name_lineage.index(rank_of_interest)
        rank_id = query_taxon.taxid_lineage[rank_index]
        return rank_id
    else:
        return None

def get_children_ids(parent):
    """return list of taxon IDs which are children of given parent taxon"""
    # taxid2parent: dict with key = child, value = parent
    return [key for key in TAX_DB.taxid2parent.keys() if TAX_DB.taxid2parent[key] == parent]

def get_gff_attribute(attr_list, attr):
    """return value for GFF attribute of given list of attributes"""
    tmp = attr_list.split(attr+'=')[1]
    if len(tmp.split(';')) > 1:
        value = strip_ID(tmp.split(';')[0])
    else:
        value = strip_ID(tmp)
    return value


###############################################################################
######################### PLOT LABEL COMPUTATION ##############################
###############################################################################


def init_label_taxa_dict(genes):
    """ initializes dictionary which stores which taxonomic assignments are summarized in each label
    every taxonomic assignment is assigned to 'root' in the beginning
    # also set initial plot label and ID (to root) for each gene """

    labelID_taxa = {1:set()} # taxa that have this label; label : set(taxa)
    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue
        labelID_taxa[1].add(gene.taxon_assignmentID)
        # gene.plot_label = 'root'
        # gene.plot_labelID = 1
    return labelID_taxa

def get_children_for_label(taxa, labelID):
    """ gets all descendants of label which are in lineage of one of the taxa the label represents
    also counts how many of these descendants are already in list of current labels """
    children = {}
    for taxID in taxa:
        taxon = taxopy.Taxon(taxID, TAX_DB)
        child = taxon.taxid_lineage[max(taxon.taxid_lineage.index(labelID)-1,0)] # max 0 if taxon is end of lineage
        if child in children:
            children[child].add(taxID)
        else:
            children[child] = set([taxID])

    return children

def merge_assignments_at_rank(genes, rank):
    """Merge taxonomic assignments at given rank.
    Args:
      genes (dict): dictionary of genes {g_name: gene object}
      rank (str): rank at which to merge all taxonomic assignments at
    """
    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue
        taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
        ranks = taxon.rank_name_dictionary
        merger_rank = ranks.get(rank)
        if merger_rank:
            gene.plot_label = merger_rank
            gene.plot_labelID = get_id_for_rank_of_species(gene.taxon_assignmentID, rank)
        else:
            gene.plot_label = gene.taxon_assignment
            gene.plot_labelID = gene.taxon_assignmentID

def merge_assignments_at_id(genes, ids):
    """Merge taxonomic assignments at given NCBI taxon IDs.

    Args:
      genes (dict): dictionary of genes {g_name: gene object}
      ids (int,list): list of taxon IDs at which to merge taxonomic assignments
    """

    if type(ids) != list:
        ids = [ids]

    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            gene.plot_label = 'Unassigned'
            gene.plot_labelID = 'NA'
            continue

        taxon = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB)
        lineage = taxon.taxid_lineage
        for id in ids:
            if id in lineage:
                id_merger_taxon = taxopy.Taxon(id, TAX_DB)
                gene.plot_label = id_merger_taxon.name
                gene.plot_labelID = id_merger_taxon.taxid
                lineage = id_merger_taxon.taxid_lineage

def iterative_merging_top(genes, num_groups_plot):
    """ identify at which rank in query lineage to merge best """

    labelID_taxa = init_label_taxa_dict(genes) # taxa that have this label; label : set(taxa)

    update_dict = {'nonempty-dummy'}
    while update_dict:
        update_dict = {}
        del_dict_keys = set()
        remaining_labels_num = num_groups_plot - len(labelID_taxa)
        for labelID, taxa in labelID_taxa.items():
            if len(taxa) >= 2:
                lca_candidates = []
                for taxon in taxa:
                    lca_candidates.append(taxopy.Taxon(taxon, TAX_DB))
                lca = taxopy.find_lca(lca_candidates, TAX_DB)
                if lca.taxid != labelID: # if it is possible to find a more specific label than the current one with the lca
                    del_dict_keys.add(labelID)
                    update_dict[lca.taxid] = taxa
                else:
                    # get all children for label that are in list taxa
                    children = get_children_for_label(taxa, labelID)
                    existing_labels = len(children.keys())-len(list(set(set(children.keys())-labelID_taxa.keys())))

                    # check if number of new labels (= num of children - labels that are already in lablels)
                    # fits in remainung number of possible new labels
                    if len(children)-existing_labels <= remaining_labels_num:
                        update_dict.update(children)
                        del_dict_keys.add(labelID)
                        remaining_labels_num -= (len(children)-existing_labels)
            else:
                pass
        for key in del_dict_keys:
            del labelID_taxa[key]
        labelID_taxa.update(update_dict)

    taxid2label = {} # {taxonomic assignment : the label it gets}
    for labelID, taxa in labelID_taxa.items():
        if len(taxa) == 1:
            taxid2label[min(taxa)] = min(taxa)
        else:
            for taxon in taxa:
                taxid2label[taxon] = labelID


    for g_name, gene in genes.items():
        if gene.taxon_assignmentID == 'NA':
            continue
        else:
            gene.plot_label = TAX_DB.taxid2name[taxid2label.get(gene.taxon_assignmentID)]
            gene.plot_labelID = taxid2label.get(gene.taxon_assignmentID)


def filter_labeled_genes(genes):
    """Filter genes which already have a label and return number of labels.

    Args:
      genes (dict): dictionary of genes {g_name: gene object}

    Returns:
      dict: dictionary of genes without label {g_name: gene object}
      int: number of distinct labels in the set of labeled genes
    """
    unlabeled_genes = {}
    labels = set()

    for g_name, gene in genes.items():
        if not gene.plot_label:
            unlabeled_genes[g_name] = gene
        else:
            labels.add(gene.plot_labelID)


    num_labels = len(labels)
    return unlabeled_genes, num_labels


def merge_labels(genes, queryID, merging_labels):
    """Merge plot labels based on user input.

    Args:
      genes (dict): dictionary of genes {g_name: gene object}
      queryID (int): NCBI taxon ID for query species
      merging_labels (list/int): user input for merging labels
        1. taxonomic rank with suffix '-all'
            -> merge all taxonomic assignment at that rank
        2. taxonomic rank
            -> merge query lineage at that rank
        3. list of NCBI taxon IDs
            -> merge at each of the taxa provided
    """

    merger_ids = []
    all_merger_rank = None
    query_merger_id = None

    if type(merging_labels) != list:
        merging_labels = [merging_labels]

    for merger_label in merging_labels:
        if merger_label.isalpha():
            # option 1: string end with '-all' -> every taxonomic assignment is merged at that rank
            # option 2: string does not end with '-all' -> only taxonomic assignments with same rank <merger_label> as the query species are merged
            if merger_label.endswith('-all'):
                all_merger_rank = merger_label[:-len('-all')]
                break
            else:
                # get ncbiID for rank <merger_label> of query
                query_merger_id = get_id_for_rank_of_species(queryID, merger_label)
                break
        else:
            if type(merger_label) == int or merger_label.isdigit():
                merger_ids.append(int(merger_label)) # ncbiID of <merger_label> to merge at




    # if there is an -all merger-label; this overwrites everything
    # else; merge at query_merger_id and merger_ids
    if all_merger_rank:
        merge_assignments_at_rank(genes, all_merger_rank)
    elif query_merger_id:
        merge_assignments_at_id(genes, query_merger_id)
    elif merger_ids:
        merge_assignments_at_id(genes, merger_ids)


def set_unassigned_labels(genes):
    """Set missing plot labels to taxonomic assignment or 'Unassigned'."""
    for g_name, gene in genes.items():
        if not gene.plot_label:
            if gene.taxon_assignmentID == 'NA':
                gene.plot_label = 'Unassigned'
                gene.plot_labelID = 'NA'
            else:
                gene.plot_label = gene.taxon_assignment
                gene.plot_labelID = gene.taxon_assignmentID


def ident_query_label(genes, queryID):
    """Identify which label represents the query species."""


    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    label_candidate = (None, len(query_lineage))

    for g_name, gene in genes.items():
        if gene.plot_labelID in query_lineage:
            if query_lineage.index(gene.plot_labelID) < label_candidate[1]:
                label_candidate = (gene.plot_label, query_lineage.index(gene.plot_labelID))

    return label_candidate[0]


###############################################################################
########################## TAXONOMIC ASSIGNMENT ###############################
###############################################################################

def get_subset_id(queryID, subset_marker):
    """Return input if input.isdigit() or return ID for rank of query."""
    if type(subset_marker) == int:
        return subset_marker
    elif subset_marker.isdigit():
        return int(subset_marker)
    else:
        return get_id_for_rank_of_species(queryID, subset_marker)


def perform_quick_search_1(perform_diamond, diamond_cmd, tax_assignment_path_1,
                            quick_mode_search_rank, taxon_exclude, queryID):
    """Perform first DIAMOND run for quick assignment mode."""
    print("Quick mode for taxonomic assignment selected")

    if perform_diamond:
        # quick_mode_search_rank can be either taxon ID or rank
        quick_mode_search_id = get_subset_id(queryID, quick_mode_search_rank)

        if taxon_exclude:
            q1_exclude = taxon_exclude.rstrip('"') + ',' + ','.join(diamond_inclusion_by_exclusion(quick_mode_search_id))
        else:
            q1_exclude = ' --taxon-exclude "{}"'.format(','.join(diamond_inclusion_by_exclusion(quick_mode_search_id)))
        diamond_cmd_1 = diamond_cmd + ' -o "{}" {}'.format(tax_assignment_path_1, q1_exclude)
        ###taxonlist###  diamond_cmd_1 = diamond_cmd + diamond_o1 + ' --taxonlist "' + str(quick_mode_search_id) + '" --taxon-exclude "' + str(queryID) + '"'

        print("Diamond command for quick search round 1: ")
        print(diamond_cmd_1)
        try:
            dmdnOut = subprocess.run([diamond_cmd_1], shell=True, capture_output=True, check=True)
        except:
            sys.exit('Error running\n%s' % diamond_cmd_1)

def perform_quick_search_2(perform_diamond, diamond_cmd, tax_assignment_path_2,
                        quick_mode_match_rank, genes, tmp_prot_path,
                        taxon_exclude, queryID):
    """Perform second DIAMOND run for quick assignment mode."""
    print("Quick mode for taxonomic assignment selected")

    if perform_diamond:
        quick_mode_match_id = get_subset_id(queryID, quick_mode_match_rank)

        filter_prots_quick_hits(genes, quick_mode_match_id, tmp_prot_path, tmp_prot_path)

        diamond_cmd_2 = diamond_cmd + ' -o "{}" '.format(tax_assignment_path_2) + taxon_exclude
        print("Diamond command for quick search round 2: ")
        print(diamond_cmd_2)
        try:
            dmdnOut = subprocess.run([diamond_cmd_2], shell=True, capture_output=True, check=True)
        except:
            sys.exit('Error running\n%s' % diamond_cmd_2)


def perform_exhaustive_search(perform_diamond, diamond_cmd,
                                tax_assignment_path, taxon_exclude):
    """Perform  DIAMOND run for exhaustive assignment mode."""
    print("Exhaustive mode for taxonomic assignment selected")
    if perform_diamond:
        diamond_cmd = diamond_cmd + ' -o "{}" '.format(tax_assignment_path) + taxon_exclude
        print("Diamond command for exhaustive search: ")
        print(diamond_cmd)
        try:
            dmdnOut = subprocess.run([diamond_cmd], shell=True, capture_output=True, check=True)
        except:
            sys.exit('Error running\n%s' % diamond_cmd)


def assess_best_of_multi_hits(taxIDs, queryID):
    """for a list of taxa IDs, get the one where LCA with queryID is closest to query"""

    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    bestID = (None, len(query_lineage)+1) #(id, index in lineage)
    for taxID in taxIDs:
        if taxID in TAX_DB.taxid2name.keys(): # avoid crashes because of ncbiIDs not being represented in TAX_DB
            taxon = taxopy.Taxon(taxID, TAX_DB)
            lca = taxopy.find_lca([query_taxon, taxon], TAX_DB)
            lca_index = query_lineage.index(lca.taxid)
            if bestID[1] > lca_index: # check if new best hit is closer to query in lineage as current best hit
                bestID = (taxon, lca_index)
        else:
            missing_taxids.add(taxID)


    return bestID[0]

def set_taxon_assignment(genes):
    """Assess taxonomic assignment for each gene.

    If there is a corrected_lca; taxonomic assignment = corrected_lca
    else; taxonomic assignment = lca
    """
    for g_name, gene in genes.items():
        if gene.corrected_lca != 'NA':
            gene.taxon_assignment = gene.corrected_lca
            gene.taxon_assignmentID = gene.corrected_lcaID
        else:
            gene.taxon_assignment = gene.lca
            gene.taxon_assignmentID = gene.lcaID


def calc_corrected_lca(genes, queryID):
    """Calculate the corrected LCA for each gene based on best hit and query."""
    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    for g_name, gene in genes.items():
        if gene.lcaID in query_lineage and gene.best_hitID != 'NA':
            if gene.best_hitID in TAX_DB.taxid2name.keys(): # avoid crashes because of ncbiIDs not being represented in TAX_DB
                bh_taxon = taxopy.Taxon(gene.best_hitID, TAX_DB)
                corrected_lca = taxopy.find_lca([query_taxon, bh_taxon], TAX_DB)
                gene.corrected_lcaID = corrected_lca.taxid
                gene.corrected_lca = corrected_lca.name
            else:
                missing_taxids.add(gene.best_hitID)
                gene.corrected_lcaID = 'NA'
                gene.corrected_lca = 'NA'
        else:
            gene.corrected_lcaID = 'NA'
            gene.corrected_lca = 'NA'


def calc_lca(genes):
    """Calculate the LCA for each gene based on taxonomic hits."""

    for g_name, gene in genes.items():
        lca_query = []
        for tax_assignment in gene.tax_assignments:
            if tax_assignment[-2]:
                tax_id = int(tax_assignment[-2].split(';')[0])
                if tax_id in TAX_DB.taxid2name.keys(): # avoid crashes because of ncbiIDs not being represented in TAX_DB
                    taxon = taxopy.Taxon(tax_id, TAX_DB)
                    lca_query.append(taxon)
                else:
                    missing_taxids.add(tax_id)
        if len(lca_query) > 1:
            lca = taxopy.find_lca(lca_query, TAX_DB)
            gene.lca = lca.name
            gene.lcaID = int(lca.taxid)
        elif not lca_query:
            gene.lca = 'NA'
            gene.lcaID = 'NA'
            pass
        else: # if only one hit for taxon
            gene.lca = lca_query[0].name
            gene.lcaID = int(lca_query[0].taxid)


def read_tax_assignments(tax_assignment_path, prots, queryID):
    """Read DIAMOND output file and assign hits to genes."""

    with open(tax_assignment_path, 'r') as tax_assignment:
        # cols: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore taxid taxname
        for line in tax_assignment:
            # first is "best hit" in terms of alignment score
            spline = line.strip().split('\t')
            gene = prots.get(spline[0])
            closest_hit = None
            if gene:
                if not gene.best_hit or gene.best_hit == 'NA': # if first best hit had no matching to NCBI ID
                    if spline[-2]:
                        if ';' in spline[-2]:
                            hit_ids = [int(id) for id in spline[-2].split(';')]
                            closest_hit = assess_best_of_multi_hits(hit_ids, queryID)
                            if closest_hit:
                                gene.best_hit = closest_hit.name
                                gene.best_hitID = closest_hit.taxid
                            else:
                                gene.best_hit = 'NA'
                                gene.best_hitID = 'NA'
                        else:
                            gene.best_hit = spline[-1]
                            gene.best_hitID = int(spline[-2])
                        gene.bh_evalue = spline[10]
                    else: # ncbiID for match could not be found
                        gene.best_hit = 'NA'
                        gene.best_hitID = 'NA'
                        gene.bh_evalue = 'NA'

                if closest_hit:
                    gene.tax_assignments.append(spline[1:-2]+[str(closest_hit.taxid),closest_hit.name])
                elif ';' in spline[-2]:
                    hit_ids = [int(id) for id in spline[-2].split(';')]
                    closest_hit = assess_best_of_multi_hits(hit_ids, queryID)
                    if closest_hit:
                        gene.tax_assignments.append(spline[1:-2]+[str(closest_hit.taxid),closest_hit.name])
                    else:
                        gene.tax_assignments.append(spline[1:])
                else:
                    gene.tax_assignments.append(spline[1:])

def taxonomic_assignment(tax_assignment_path, genes, prots, queryID):
    """Run functions to save taxonomic assignment to gene objects"""
    read_tax_assignments(tax_assignment_path, prots, queryID)
    calc_lca(genes)
    calc_corrected_lca(genes, queryID)
    set_taxon_assignment(genes)

def reset_tax_assignment(gene):
    """Delete all taxonomc assignment info for given gene"""
    gene.tax_assignments = []
    gene.best_hit = None
    gene.best_hitID = None
    gene.lca = None
    gene.lcaID = None
    gene.corrected_lca = None
    gene.corrected_lcaID = None
    gene.taxon_assignment = None
    gene.taxon_assignmentID = None



###############################################################################
############################# DATA PREPARATION ################################
###############################################################################


def get_child_parent(attrs, child_parent_dict):
    """Use GFF attributes and retrieve child-parent pairs, given certain attributes."""
    child_attr = child_parent_dict.get('child')
    parent_attr = child_parent_dict.get('parent')
    if child_attr+"=" in attrs and parent_attr+"=" in attrs:
        childID = get_gff_attribute(attrs, child_attr)
        parentID = get_gff_attribute(attrs, parent_attr)
    else:
        return None, None
    return childID, parentID


def prot_gene_matching(output_path, gff_path, genes, gff_rule):
    """Match headers in protein FASTA to gene ID.

    Parses the GFF using the given gff parsing rule
    """

    prots = {}

    child_parent_dict = {} # key = child, value = parent
    with open(gff_path, 'r') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                spline = line.strip().split('\t')
                # read lines where fasta headers are contained and/or child-parent-relations are to be found
                if gff_rule.get('source') == "default" or spline[1] == gff_rule.get('source'):
                    if spline[2] in gff_rule.get('fasta_header_type'):
                        if gff_rule.get('gene_connection') == "parent/child":
                            # when attribute for child-parent-matching unequals attribute where header can be found,
                            # match attribute of fasta header and child of parent/child matching together
                            if gff_rule.get('parent_child_attr').get('child') not in gff_rule.get('fasta_header_attr'):
                                header_child_dict = {'child': gff_rule.get('fasta_header_attr'), 'parent': gff_rule.get('parent_child_attr').get('child')}
                                childID, parentID = get_child_parent(spline[8], header_child_dict)
                                if childID:
                                    if not child_parent_dict.get(strip_ID(childID)):
                                        child_parent_dict[strip_ID(childID)] = strip_ID(parentID)
                                    elif child_parent_dict.get(strip_ID(childID)) == strip_ID(childID):
                                        child_parent_dict[strip_ID(childID)] = strip_ID(parentID)

                            # add the parent child matching if fasta header line (also) contains p/c matching info
                            if spline[2] in gff_rule.get('parent_child_types'):
                                childID, parentID = get_child_parent(spline[8], gff_rule.get('parent_child_attr'))
                                if childID:
                                    if not child_parent_dict.get(strip_ID(childID)):
                                        child_parent_dict[strip_ID(childID)] = strip_ID(parentID)
                                    elif child_parent_dict.get(strip_ID(childID)) == strip_ID(childID):
                                        child_parent_dict[strip_ID(childID)] = strip_ID(parentID)
                        elif gff_rule.get('gene_connection') == "inline":
                            headerID = get_gff_attribute(spline[8], gff_rule.get('fasta_header_attr'))
                            geneID = get_gff_attribute(spline[8], gff_rule.get('gene_attr'))
                            if not child_parent_dict.get(strip_ID(headerID)):
                                child_parent_dict[strip_ID(headerID)] = strip_ID(geneID)
                            elif child_parent_dict.get(strip_ID(headerID)) == strip_ID(headerID):
                                child_parent_dict[strip_ID(headerID)] = strip_ID(geneID)
                    elif spline[2] in gff_rule.get('parent_child_types'):
                        childID, parentID = get_child_parent(spline[8], gff_rule.get('parent_child_attr'))
                        if childID:
                            if not child_parent_dict.get(strip_ID(childID)):
                                child_parent_dict[strip_ID(childID)] = strip_ID(parentID)
                            elif child_parent_dict.get(strip_ID(childID)) == strip_ID(childID):
                                child_parent_dict[strip_ID(childID)] = strip_ID(parentID)

            elif "#FASTA" in line: # if FASTA block has been reached
                break


    prot_index_path = output_path + 'tmp/tmp.proteins.fa.fai'

    unmatched = [] #list of proteins where gene could not be matched

    with open(prot_index_path, 'r') as prot_index:
        for line in prot_index:
            id, length = line.split()[0], ((int(line.split()[1])*3)+3)
            parent = child_parent_dict.get(strip_ID(id))

            if parent: # id in proteins fasta index matches the IDs in the GFF
                while parent in child_parent_dict.keys():
                    if parent == child_parent_dict.get(parent):
                        break
                    else:
                        parent = child_parent_dict.get(parent)
            elif child_parent_dict.get(strip_ID(id.split('|')[0])): # header might contain pipes
                parent = child_parent_dict.get(id.split('|')[0])
                while parent in child_parent_dict.keys():
                    if parent == child_parent_dict.get(parent):
                        break
                    else:
                        parent = child_parent_dict.get(parent)
            elif strip_ID(id) in genes.keys():
                # fasta header is gene ID
                parent = strip_ID(id)
            elif strip_ID(id).split('|')[0] in genes.keys():
                parent = strip_ID(id)
            else: # if not, look if any child ID is contained in the id
                for child in child_parent_dict.keys():
                    if child in id:
                        parent = child_parent_dict.get(child)
                        while parent in child_parent_dict.keys():
                            if parent == child_parent_dict.get(parent):
                                break
                            else:
                                parent = child_parent_dict.get(parent)
                if not parent: # still no success try to see if you can find gene ID in id
                    for g_name in genes.keys():
                        if g_name in id:
                            parent = g_name
            if not parent:
                unmatched.append(id)


            else:
                gene = genes.get(parent)
                # put the ID as the protID, as this will be what is used in DIAMOND results
                if gene:
                    if gene.protID:
                        if length > gene.lencds:
                            gene.protID = id
                            gene.lencds = length
                            prots[gene.protID] = gene
                    else:
                        gene.protID = id
                        gene.lencds = length
                        prots[gene.protID] = gene


    if len(prots) == 0:
        print("Proteins were unable to be matched with gene IDs via GFF. Please check that FASTA headers fit with IDs in GFF (string before first whitespace or pipe must match ID of mRNA or CDS of corresponding gene)")
    else:
        if unmatched:
            print("Protein(s) with FASTA following header could not be matched to a gene ID:" )
            for id in unmatched:
                print(id)
    return prots


def parse_header(header):
    """Returns dictionary mapping column names to index in the header of gene_table_coords."""

    indexing = {}
    for index, item in enumerate(header.split(',')):
        raw_item = item.strip().rstrip('_0123456789')
        if raw_item in indexing.keys(): # for coverage variables
            indexing[raw_item] = indexing.get(raw_item) + (index,)
        elif raw_item == 'Dim.' or ('cov' in raw_item and 'bool' not in raw_item):
            indexing[raw_item] = (index,)
        else:
            indexing[raw_item] = index
    return indexing


def read_genes_coords(output_path):
    """Read gene_table_coords.csv and save information to Gene objects."""

    input_path = output_path + 'PCA_and_clustering/gene_table_coords.csv'

    with open(input_path, 'r') as input_file:
        header = next(input_file)
        header_index = parse_header(header)

        genes = {}
        for line in input_file:
            gene = Gene(line.strip().split(','), header_index)

            # filter genes without PCA coordinates
            if gene.coords[0] != "NA":
                genes[gene.g_name] = gene

    return genes, header


###############################################################################
########################### PROTEIN FASTA FILE ################################
###############################################################################


def subset_protein_fasta(proteins_path, prot_list, path_out, mode):
    """use the list of IDs to write to the tmp prot fasta only those
    proteins whose ID appears (not) in the list
    """

    # mode: include or exclude list
    if mode == 'include':
        mode = True
    else:
        mode = False

    with open(proteins_path, 'r') as file_prot:
        lines = file_prot.readlines()
    with open(path_out, 'w') as file_out:
        write_bool = False
        for line in lines:
            if line.startswith('>'):
                transcript = line.split()[0].lstrip(">")
                #gene = line.split()[0].split('|')[1].lstrip(">")

                if ((transcript in prot_list) and mode) or ((transcript not in prot_list) and not mode):
                    write_bool = True
                    file_out.write(line)
                else:
                    write_bool = False
            else:
                if write_bool:
                    file_out.write(line)


def subset_prots_longest_cds(genes,proteins_path, path_out):
    """ the gene-protein ID matching table, generated with gffread is used for this table
    (it contains the length of the CDS ) to infer the transcript with the longest cds for each gene
    this transcript is the written to a tmp fasta file"""

    # when matching prot ID to gene ID it is already checked for the one with the longest CDS
    longest_transcripts = [gene.protID for gene in genes.values() if gene.protID]

    print(str(len(longest_transcripts)) + " proteins written to fasta file for taxonomic assignment (longest cds subsetting)")
    subset_protein_fasta(proteins_path, longest_transcripts, path_out, "include")


def filter_prots_quick_hits(genes, quick_mode_match_id, prot_path_in, prot_path_out):
    """write proteins to tmp protein FASTA which were not assigned in first
    quick search with DIAMOND
    """

    no_match = []

    for g_name, gene in genes.items():
        if gene.protID:
            if gene.taxon_assignmentID != 'NA':

                tax_lineage = taxopy.Taxon(gene.taxon_assignmentID, TAX_DB).taxid_lineage
                #if quick_mode_match_id not in tax_lineage or quick_mode_match_id == tax_lineage[0]: # <- option 2: exclusive
                if quick_mode_match_id not in tax_lineage:                                            # <- option 1: inclusive
                    # inclusive: if match id equals class, assignment must be at least same class as query
                    # exclusive: if match id equals class, assignment must be at least one rank closer to query
                    no_match.append(gene.protID)
                    reset_tax_assignment(gene)
            else:
                no_match.append(gene.protID)
                reset_tax_assignment(gene)

    print(str(len(no_match)) + " proteins written to fasta file for taxonomic assignment (quick search subsetting)")
    subset_protein_fasta(prot_path_in, no_match, prot_path_out, "include")


###############################################################################
################################### MISC ######################################
###############################################################################


def diamond_inclusion_by_exclusion(include_id):
    """ diamond can not exclude and include taxa simultaneously
    thus, this script identifies all taxIDs to exclude to simulate inclusion
    """

    exclude_list = []

    include_taxon = taxopy.Taxon(include_id, TAX_DB)
    include_lineage = include_taxon.taxid_lineage
    for parent in include_lineage[1:]:
        childs = get_children_ids(parent)
        exclude_list = exclude_list + [str(child) for child in childs if child not in include_lineage]

    return exclude_list


def write_output(output_path, genes, header):
    """Write information in gene objects to gene_table_taxon_assignment.csv."""

    out_path = output_path + 'taxonomic_assignment/gene_table_taxon_assignment.csv'

    with open(out_path, 'w') as out_file:

        csv_columns = header.strip().split(',') + ['protID', 'lcaID', 'lca',  'best_hitID', 'best_hit', 'bh_evalue', 'corrected_lca', 'taxon_assignment', 'plot_label']
        writer = csv.DictWriter(out_file, fieldnames=csv_columns, extrasaction='ignore')
        writer.writeheader()

        for gene in genes.values():
            gene_dict = gene.__dict__
            to_delete = []
            for attr, value in list(gene_dict.items()):
                if type(value) == tuple:
                    if attr == "coords":
                        name = "Dim."
                        add = 1
                    else:
                        name = attr + "_"
                        add = 0
                    for index, sub_val in enumerate(value):
                        gene_dict[name + str(index+add)] = sub_val
                    to_delete.append(attr)

            for key in to_delete:
                del gene_dict[key]

            writer.writerow(gene_dict)


def parse_gff_source_rule(gff_source):
    """ parsing which features in GFF to read (as genes) """
    if gff_source == "default":
        rule_dict = {}

    elif gff_source == "maker":
        rule_dict = {'source': 'maker',
            'gene_tag': 'gene',
            'fasta_header_type': 'mRNA',
            'fasta_header_attr': 'ID',
            'gene_connection': 'parent/child',
            'parent_child_types': 'mRNA',
            'parent_child_attr': {'parent': 'Parent', 'child': 'ID'}}
    elif gff_source == "augustus_masked":
        rule_dict = {'source': 'augustus_masked',
            'gene_tag': 'match',
            'fasta_header_type': 'match',
            'fasta_header_attr': 'Name',
            'gene_connection': 'inline',
            'gene_attr': 'ID'}

    else:
        rule_file_path = pathlib.Path(gff_source)
        if rule_file_path.is_file():
            rule_dict = {}
            with open(rule_file_path, 'r') as rule_file:
                for line in rule_file:
                    rule_dict[line.split(":")[0].strip()] = ":".join(line.split(":")[1:]).strip()
            if "parent_child_attr" in rule_dict.keys():
                pc_dict = {}
                for kv_pair in rule_dict.get("parent_child_attr").strip("{}").split(","):
                    pc_dict[kv_pair.split(":")[0].strip()] = kv_pair.split(":")[1].strip()
                rule_dict["parent_child_attr"] = pc_dict
        else:
            rule_dict = {}
            print("ERROR: source type for GFF could not be interpreted. Please check your input. Computations will continue with default setting")

    if not rule_dict: #if not been set -> set to default
        rule_dict = {'source': 'default',
            'gene_tag': 'gene',
            'fasta_header_type': 'mRNA,CDS',
            'fasta_header_attr': 'ID',
            'gene_connection': 'parent/child',
            'parent_child_types': 'mRNA,CDS',
            'parent_child_attr': {'parent': 'Parent', 'child': 'ID'}}


    return rule_dict

def write_tmp_query_info(output_path, genes, queryID):
    """Write temporary file which has information about query taxon."""
    query_label = ident_query_label(genes, queryID)
    query_name = TAX_DB.taxid2name[queryID]
    with open(output_path+'tmp/tmp.query_label', 'w') as tmp_file:
        tmp_file.write(query_label+'\n')
        tmp_file.write(query_name+'\n')



###############################################################################
###############################################################################

def main():

    config_path = sys.argv[1]

    # read parameters from config file
    config_obj=yaml.safe_load(open(config_path,'r'))
    output_path=config_obj.get('output_path') # complete output path (ENDING ON A SLASH!)
    nr_db_path=config_obj.get('database_path') # path to database
    proteins_path=config_obj.get('proteins_path', None) # path to FASTA w/ protein seqs
    gff_path=config_obj.get('gff_path') # path to GFF
    taxon_exclude = config_obj.get('taxon_exclude') # bool to exclude query taxon from sequence alignment
    compute_tax_assignment = config_obj.get('compute_tax_assignment') # whether taxonomic assigned should be computed
    only_plotting = config_obj.get('update_plots') # whether only the plots should be update
    # -> recompute the taxonomic assignment based on exisiting alignment hits as settings might be changed
    assignment_mode = config_obj.get('assignment_mode')
    quick_mode_search_rank = config_obj.get('quick_mode_search_rank', None)
    quick_mode_match_rank = config_obj.get('quick_mode_match_rank', None)
    if assignment_mode == "quick":
        tax_assignment_path_1, tax_assignment_path_2 = tax_assignment_path = config_obj.get('tax_assignment_path')
    else:
        tax_assignment_path = config_obj.get('tax_assignment_path')
    queryID = int(config_obj.get('taxon_id'))
    merging_labels = config_obj.get('merging_labels')
    num_groups_plot = config_obj.get('num_groups_plot')
    gff_source = config_obj.get('gff_source', 'default')

    gff_rule = parse_gff_source_rule(gff_source)

    tmp_prot_path = output_path+"tmp/tmp.subset.protein.fasta"

    diamond_cmd = 'diamond blastp -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames -b2.0 --tmpdir /dev/shm --sensitive -c1 --top 10 -q "{}" -d "{}"'.format(tmp_prot_path, nr_db_path)

    if taxon_exclude:
        taxon_exclude = ' --taxon-exclude "{}"'.format(queryID)
    else:
        taxon_exclude = ''
    if type(num_groups_plot) == int:
        pass
    elif num_groups_plot.isdigit():
        num_groups_plot = int(num_groups_plot)
    elif num_groups_plot == 'all':
        # no merging of taxonomic assignments desired
        num_groups_plot = False
    else:
        print('No valid option for "num_groups_plot"')
        return

    # read file
    genes, header = read_genes_coords(output_path)

    prots = prot_gene_matching(output_path, gff_path, genes, gff_rule)

    perform_diamond = compute_tax_assignment and not only_plotting

    if perform_diamond:
        subset_prots_longest_cds(genes, proteins_path, tmp_prot_path)

    if assignment_mode == 'quick':
        perform_quick_search_1(perform_diamond, diamond_cmd,
                                tax_assignment_path_1, quick_mode_search_rank,
                                taxon_exclude, queryID)
        taxonomic_assignment(tax_assignment_path_1, genes, prots, queryID)
        perform_quick_search_2(perform_diamond, diamond_cmd,
                                tax_assignment_path_2, quick_mode_match_rank,
                                genes, tmp_prot_path, taxon_exclude, queryID)
        taxonomic_assignment(tax_assignment_path_2, genes, prots, queryID)
    elif assignment_mode == 'exhaustive':
        perform_exhaustive_search(perform_diamond, diamond_cmd,
                                    tax_assignment_path, taxon_exclude)
        taxonomic_assignment(tax_assignment_path, genes, prots, queryID)
    else:
        print('Assignment mode not one of quick or exhaustive')

    merge_labels(genes, queryID, merging_labels)
    if num_groups_plot:
        unlabeled_genes, num_labels = filter_labeled_genes(genes)
        iterative_merging_top(unlabeled_genes, num_groups_plot-num_labels)
    else:
        set_unassigned_labels(genes)

    write_tmp_query_info(output_path, genes, queryID)

    if missing_taxids:
        print("The following Taxon ID(s) could not be found in the NCBI: ")
        print(missing_taxids)
        print("Skipped for taxonomic assignment.")

    pathlib.Path(output_path+'taxonomic_assignment/').mkdir(parents=True, exist_ok=True)

    write_output(output_path, genes, header)


if __name__ == '__main__':
    main()
