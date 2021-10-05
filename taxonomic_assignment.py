# -*- coding: utf-8 -*-

import yaml # read config file
import taxopy
import time
import sys
import numpy as np
import csv
import os
import pathlib


TAX_DB = taxopy.TaxDb(keep_files=True)
# TAX_DB = taxopy.TaxDb(nodes_dmp="./nodes.dmp", names_dmp="./names.dmp", keep_files=True)

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

        self.c_cov = tuple([line_array[index] for index in header_index.get('c_cov')])
        self.c_covsd = tuple([line_array[index] for index in header_index.get('c_covsd')])
        self.c_covdev = tuple([line_array[index] for index in header_index.get('c_covdev')])
        self.c_genecovm = tuple([line_array[index] for index in header_index.get('c_genecovm')])
        self.c_genecovsd = tuple([line_array[index] for index in header_index.get('c_genecovsd')])

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

        self.g_cov = tuple([line_array[index] for index in header_index.get('g_cov')])
        self.g_covsd = tuple([line_array[index] for index in header_index.get('g_covsd')])
        self.g_covdev_c = tuple([line_array[index] for index in header_index.get('g_covdev_c')])
        self.g_covdev_o = tuple([line_array[index] for index in header_index.get('g_covdev_o')])
        self.g_cov_zscore = tuple([line_array[index] for index in header_index.get('g_cov_zscore')])
        self.g_cov_z_bool = line_array[header_index.get('g_cov_z_bool')]

        self.g_pearson_r_o = line_array[header_index.get('g_pearson_r_o')]
        self.g_pearson_p_o = line_array[header_index.get('g_pearson_p_o')]
        self.g_pearson_r_c = line_array[header_index.get('g_pearson_r_c')]
        self.g_pearson_p_c = line_array[header_index.get('g_pearson_p_c')]
        self.g_gc_cont = line_array[header_index.get('g_gc_cont')]
        self.g_gcdev_c = line_array[header_index.get('g_gcdev_c')]
        self.g_gcdev_o = line_array[header_index.get('g_gcdev_o')]

        self.coords = tuple([line_array[index] for index in header_index.get('Dim.')])

        self.protID = None # ID of transcript with longest CDS
        self.lencds = 0 # lenght of self.protID

        self.tax_assignments = []
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



# ██   ██ ███████ ██      ██████  ███████ ██████      ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
# ██   ██ ██      ██      ██   ██ ██      ██   ██     ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
# ███████ █████   ██      ██████  █████   ██████      █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
# ██   ██ ██      ██      ██      ██      ██   ██     ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
# ██   ██ ███████ ███████ ██      ███████ ██   ██     ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████

def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text

def strip_ID(id):
    """ remove GFFs prefixes from IDs """
    for prefix in ['gene:','gene-','transcript:','transcript-','rna:','rna-','cds:','cds-']:
        id = remove_prefix(id,prefix)
    return id


def get_id_for_rank_of_species(queryID, rank):
    """ for a query taxon and a rank; get the ID of the taxon at the given rank for given taxon """
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
    # get children for parent
    # taxid2parent: dict with key = child, value = parent
    return [key for key in TAX_DB.taxid2parent.keys() if TAX_DB.taxid2parent[key] == parent]

def get_gff_attribute(attr_list, attr):
    """ stable version of: strip_ID(spline[8].split('ID=')[1].split(';')[0]) """
    tmp = attr_list.split(attr+'=')[1]
    if len(tmp.split(';')) > 1:
        value = strip_ID(tmp.split(';')[0])
    else:
        value = strip_ID(tmp)
    return value


# ██████  ██       ██████  ████████     ██       █████  ██████  ███████ ██      ███████
# ██   ██ ██      ██    ██    ██        ██      ██   ██ ██   ██ ██      ██      ██
# ██████  ██      ██    ██    ██        ██      ███████ ██████  █████   ██      ███████
# ██      ██      ██    ██    ██        ██      ██   ██ ██   ██ ██      ██           ██
# ██      ███████  ██████     ██        ███████ ██   ██ ██████  ███████ ███████ ███████


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
    """  """
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
    for g_name, gene in genes.items():
        if not gene.plot_label:
            if gene.taxon_assignmentID == 'NA':
                gene.plot_label = 'Unassigned'
                gene.plot_labelID = 'NA'
            else:
                gene.plot_label = gene.taxon_assignment
                gene.plot_labelID = gene.taxon_assignmentID


def ident_query_label(genes, queryID):
    """ identify which label represent the query species """

    query_taxon = taxopy.Taxon(queryID, TAX_DB)
    query_lineage = query_taxon.taxid_lineage

    label_candidate = (None, len(query_lineage))

    for g_name, gene in genes.items():
        if gene.plot_labelID in query_lineage:
            if query_lineage.index(gene.plot_labelID) < label_candidate[1]:
                label_candidate = (gene.plot_label, query_lineage.index(gene.plot_labelID))


    return label_candidate[0]


# ████████  █████  ██   ██  ██████  ███    ██  ██████  ███    ███ ██  ██████      █████  ███████ ███████ ██  ██████  ███    ██ ███    ███ ███████ ███    ██ ████████
#    ██    ██   ██  ██ ██  ██    ██ ████   ██ ██    ██ ████  ████ ██ ██          ██   ██ ██      ██      ██ ██       ████   ██ ████  ████ ██      ████   ██    ██
#    ██    ███████   ███   ██    ██ ██ ██  ██ ██    ██ ██ ████ ██ ██ ██          ███████ ███████ ███████ ██ ██   ███ ██ ██  ██ ██ ████ ██ █████   ██ ██  ██    ██
#    ██    ██   ██  ██ ██  ██    ██ ██  ██ ██ ██    ██ ██  ██  ██ ██ ██          ██   ██      ██      ██ ██ ██    ██ ██  ██ ██ ██  ██  ██ ██      ██  ██ ██    ██
#    ██    ██   ██ ██   ██  ██████  ██   ████  ██████  ██      ██ ██  ██████     ██   ██ ███████ ███████ ██  ██████  ██   ████ ██      ██ ███████ ██   ████    ██


def assess_best_of_multi_hits(taxIDs, queryID):
    """ for a list of taxa IDs, get the one where LCA with queryID is closest to query """

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
    # assess taxonomic assignment for each gene:
    # if there is a corrected_lca; taxonomic assignment = corrected_lca
    # else; taxonomic assignment = lca
    for g_name, gene in genes.items():
        if gene.corrected_lca != 'NA':
            gene.taxon_assignment = gene.corrected_lca
            gene.taxon_assignmentID = gene.corrected_lcaID
        else:
            gene.taxon_assignment = gene.lca
            gene.taxon_assignmentID = gene.lcaID


def calc_corrected_lca(genes, queryID):

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

    with open(tax_assignment_path, 'r') as tax_assignment:
        # cols: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore taxid taxname
        for line in tax_assignment:
            # first is "best hit" in terms of alignment score
            spline = line.strip().split('\t')
            gene = prots.get(spline[0])
            closest_hit = None
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
    read_tax_assignments(tax_assignment_path, prots, queryID)
    calc_lca(genes)
    calc_corrected_lca(genes, queryID)
    set_taxon_assignment(genes)

def reset_tax_assignment(gene):
    gene.tax_assignments = []
    gene.best_hit = None
    gene.best_hitID = None
    gene.lca = None
    gene.lcaID = None
    gene.corrected_lca = None
    gene.corrected_lcaID = None
    gene.taxon_assignment = None
    gene.taxon_assignmentID = None



# ██████   █████  ████████  █████      ██████  ██████  ███████ ██████   █████  ██████   █████  ████████ ██  ██████  ███    ██
# ██   ██ ██   ██    ██    ██   ██     ██   ██ ██   ██ ██      ██   ██ ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██
# ██   ██ ███████    ██    ███████     ██████  ██████  █████   ██████  ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██
# ██   ██ ██   ██    ██    ██   ██     ██      ██   ██ ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
# ██████  ██   ██    ██    ██   ██     ██      ██   ██ ███████ ██      ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████



def prot_gene_matching(output_path, gff_path, genes, gff_rule):

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
                            # match those attributes together
                            if gff_rule.get('parent_child_attr').get('child') != gff_rule.get('fasta_header_attr'):
                                child_attr = gff_rule.get('fasta_header_attr')
                                parent_attr = gff_rule.get('parent_child_attr').get('child')
                                if (spline[8].startswith(child_attr+"=") or ";"+child_attr+"=" in spline[8]) and \
                                   (spline[8].startswith(parent_attr+"=") or ";"+parent_attr+"=" in spline[8]):
                                    childID = get_gff_attribute(spline[8], child_attr)
                                    parentID = get_gff_attribute(spline[8], parent_attr)
                                    child_parent_dict[childID] = parentID
                            if gff_rule.get('fasta_header_type') in gff_rule.get('parent_child_types'):
                                child_attr = gff_rule.get('parent_child_attr').get('child')
                                parent_attr = gff_rule.get('parent_child_attr').get('parent')
                                if (spline[8].startswith(child_attr+"=") or ";"+child_attr+"=" in spline[8]) and \
                                   (spline[8].startswith(parent_attr+"=") or ";"+parent_attr+"=" in spline[8]):
                                    childID = get_gff_attribute(spline[8], child_attr)
                                    parentID = get_gff_attribute(spline[8], parent_attr)
                                    child_parent_dict[childID] = parentID
                        elif gff_rule.get('gene_connection') == "inline":
                            headerID = get_gff_attribute(spline[8], gff_rule.get('fasta_header_attr'))
                            geneID = get_gff_attribute(spline[8], gff_rule.get('gene_attr'))
                            child_parent_dict[headerID] = geneID
                    elif gff_rule.get('parent_child_types') and spline[2] in gff_rule.get('parent_child_types'):
                        child_attr = gff_rule.get('parent_child_attr').get('child')
                        parent_attr = gff_rule.get('parent_child_attr').get('parent')
                        if (spline[8].startswith(child_attr+"=") or ";"+child_attr+"=" in spline[8]) and \
                           (spline[8].startswith(parent_attr+"=") or ";"+parent_attr+"=" in spline[8]):
                            childID = get_gff_attribute(spline[8], child_attr)
                            parentID = get_gff_attribute(spline[8], parent_attr)
                            child_parent_dict[childID] = parentID

            elif "#FASTA" in line: # if FASTA block has been reached
                break

    prot_index_path = output_path + 'tmp/tmp.proteins.fa.fai'

    unmatched = [] #list of proteins where gene could not be matched

    with open(prot_index_path, 'r') as prot_index:
        for line in prot_index:
            id, length = line.split()[0], ((int(line.split()[1])*3)+3)
            parent = child_parent_dict.get(id)
            if parent: # id in proteins fasta index matches the IDs in the GFF
                while parent in child_parent_dict.keys():
                    if parent == child_parent_dict.get(parent):
                        break
                    else:
                        parent = child_parent_dict.get(parent)
            elif child_parent_dict.get(id.split('|')[0]): # header might contain pipes
                parent = child_parent_dict.get(id.split('|')[0])
                while parent in child_parent_dict.keys():
                    if parent == child_parent_dict.get(parent):
                        break
                    else:
                        parent = child_parent_dict.get(parent)
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



# ██████  ██████   ██████  ████████ ███████ ██ ███    ██     ███████  █████  ███████ ████████  █████
# ██   ██ ██   ██ ██    ██    ██    ██      ██ ████   ██     ██      ██   ██ ██         ██    ██   ██
# ██████  ██████  ██    ██    ██    █████   ██ ██ ██  ██     █████   ███████ ███████    ██    ███████
# ██      ██   ██ ██    ██    ██    ██      ██ ██  ██ ██     ██      ██   ██      ██    ██    ██   ██
# ██      ██   ██  ██████     ██    ███████ ██ ██   ████     ██      ██   ██ ███████    ██    ██   ██

def subset_protein_fasta(proteins_path, prot_list, path_out, mode):
    # use the list of IDs to write to the tmp prot fasta only those
    # proteins whose ID appears in (not) the list

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


# ███    ███ ██ ███████  ██████
# ████  ████ ██ ██      ██
# ██ ████ ██ ██ ███████ ██
# ██  ██  ██ ██      ██ ██
# ██      ██ ██ ███████  ██████


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



# ███    ███  █████  ██ ███    ██
# ████  ████ ██   ██ ██ ████   ██
# ██ ████ ██ ███████ ██ ██ ██  ██
# ██  ██  ██ ██   ██ ██ ██  ██ ██
# ██      ██ ██   ██ ██ ██   ████


def main():

    config_path = sys.argv[1]

    # read parameters from config file
    config_obj=yaml.safe_load(open(config_path,'r'))
    output_path=config_obj['output_path'] # complete output path (ENDING ON A SLASH!)
    nr_db_path=config_obj['database_path']
    proteins_path=config_obj['proteins_path'] # path to FASTA w/ protein seqs
    gff_path=config_obj['gff_path'] # path to GFF
    taxon_exclude = config_obj['taxon_exclude'] # bool to exclude query taxon from sequence alignment
    compute_tax_assignment = config_obj['compute_tax_assignment']
    only_plotting = config_obj['update_plots']
    assignment_mode = config_obj['assignment_mode']
    quick_mode_search_rank = config_obj['quick_mode_search_rank'] if 'quick_mode_search_rank' in config_obj.keys() else None
    quick_mode_match_rank = config_obj['quick_mode_match_rank'] if 'quick_mode_match_rank' in config_obj.keys() else None
    tax_assignment_path = config_obj['tax_assignment_path']
    queryID = int(config_obj['taxon_id'])
    merging_labels = config_obj['merging_labels']
    num_groups_plot = config_obj['num_groups_plot']
    gff_source = config_obj['gff_source'] if 'gff_source' in config_obj.keys() else "default"

    gff_rule = parse_gff_source_rule(gff_source)

    tmp_prot_path = output_path+"tmp/tmp.subset.protein.fasta"

    diamond_q = ' -q "' + tmp_prot_path + '"'
    raw_tax_assignment_path = '.'.join(tax_assignment_path.split('.')[:-1])
    diamond_o = ' -o "' + tax_assignment_path + '"'
    diamond_o1 = ' -o "' + raw_tax_assignment_path + ".quick_1.txt" + '"'
    diamond_o2 = ' -o "' + raw_tax_assignment_path + ".quick_2.txt" + '"'
    diamond_d = ' -d "' + nr_db_path + '"'
    diamond_cmd = 'diamond blastp -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames -b2.0 --tmpdir /dev/shm --sensitive -c1 --top 10' + diamond_q + diamond_d

    if taxon_exclude.upper() == "TRUE":
        taxon_exclude = ' --taxon-exclude "' + str(queryID) + '"'
    else:
        taxon_exclude = ''

    if compute_tax_assignment.upper() == "TRUE":
        compute_tax_assignment = True
    else:
        compute_tax_assignment = False

    if only_plotting.upper() == "TRUE":
        only_plotting = True
    else:
        only_plotting = False

    if num_groups_plot.isdigit() or type(num_groups_plot) == int:
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

    if compute_tax_assignment and not only_plotting:
        subset_prots_longest_cds(genes, proteins_path, tmp_prot_path)

    if assignment_mode == 'quick':
        print("Quick mode for taxonomic assignment selected")

        # TODO: check if file exists
        if compute_tax_assignment and not only_plotting:
            if quick_mode_search_rank.isdigit() or type(quick_mode_search_rank) == int:
                quick_mode_search_id = int(quick_mode_search_rank)
            else:
                quick_mode_search_id = get_id_for_rank_of_species(queryID, quick_mode_search_rank)
            exclude_ids_string = ','.join(diamond_inclusion_by_exclusion(quick_mode_search_id))
            diamond_cmd_1 = diamond_cmd + diamond_o1 + taxon_exclude.rstrip('"') + ',' + exclude_ids_string + '"'
            ###taxonlist###  diamond_cmd_1 = diamond_cmd + diamond_o1 + ' --taxonlist "' + str(quick_mode_search_id) + '" --taxon-exclude "' + str(queryID) + '"'


            print("Diamond command for quick search round 1: ")
            print(diamond_cmd_1)
            os.system(diamond_cmd_1)

        taxonomic_assignment(raw_tax_assignment_path + ".quick_1.txt", genes, prots, queryID)

        # write_output(output_path, genes, header)

        if compute_tax_assignment and not only_plotting:
            if quick_mode_match_rank.isdigit() or type(quick_mode_match_rank) == int:
                quick_mode_match_id = int(quick_mode_match_rank)
            else:
                quick_mode_match_id = get_id_for_rank_of_species(queryID, quick_mode_match_rank)
            filter_prots_quick_hits(genes, quick_mode_match_id, tmp_prot_path, tmp_prot_path)

            diamond_cmd_2 = diamond_cmd + diamond_o2 + taxon_exclude
            print("Diamond command for quick search round 2: ")
            print(diamond_cmd_2)
            os.system(diamond_cmd_2)

        taxonomic_assignment(raw_tax_assignment_path + ".quick_2.txt", genes, prots, queryID)


    elif assignment_mode == 'exhaustive':
        print("Exhaustive mode for taxonomic assignment selected")
        if compute_tax_assignment:
            diamond_cmd = diamond_cmd + diamond_o + taxon_exclude
            print("Diamond command for exhaustive search: ")
            print(diamond_cmd)
            os.system(diamond_cmd)

        taxonomic_assignment(tax_assignment_path, genes, prots, queryID)


    else:
        print('Assignment mode not one of quick or exhaustive')


    merge_labels(genes, queryID, merging_labels)
    if num_groups_plot:
        unlabeled_genes, num_labels = filter_labeled_genes(genes)
        iterative_merging_top(unlabeled_genes, num_groups_plot-num_labels)
    else:
        set_unassigned_labels(genes)

    query_label = ident_query_label(genes, queryID)
    with open(output_path+'tmp/tmp.query_label', 'w') as tmp_file:
        tmp_file.write(query_label+'\n')

    if missing_taxids:
        print("The following Taxon ID(s) could not be found in the NCBI: ")
        print(missing_taxids)
        print("Skipped for taxonomic assignment.")

    write_output(output_path, genes, header)


if __name__ == '__main__':
    main()
