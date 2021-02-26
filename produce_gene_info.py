# -*- coding: utf-8 -*-

# reads FASTA, GFF3 and per_base_coverage.txt file and produces gene_info directory
# author: Simonida Zehr
# date: 18 June 2019

from operator import itemgetter
import numpy as np
import pathlib # to create directories
import operator # for quick comparisons
import scipy.stats as stats # for Pearson's R
from itertools import product as itertools_product # to generate all possible oligonucleotides from base alphabet
from Bio.Seq import Seq as BioPython_Seq # to count oligonucleotides (also overlapping ones! not implemented in normal count)

# read parameters from config file ./cfg.json
# json version: config_obj=json.load(open('./cfg.json','r'))
import yaml
config_obj=yaml.safe_load(open('./config.yml','r'))
gff_path=config_obj['gff_path'] # GFF file path
pbc_paths=config_obj['pbc_paths'] # per base coverage (PBC) file path(s)
contig_tag=config_obj['contig_tag'] # contig label (e.g. "region", "contig", ...)
output_path=config_obj['output_path'] # complete output path (ENDING ON A SLASH!)
fasta_path=output_path+'tmp.MILTS.fasta' # path to tmp FASTA file
include_pseudogenes = config_obj['include_pseudogenes'] # boolean signifying whether pseudogenes should be included in the analysis
include_coverage = config_obj['include_coverage']

# read pseudogene decision:
# note: trying to catch all user inputs such as 'TRUE', 'tru', True', 't', 'T', ...
# --> read only the first letter
# .lower() transforms the given string to all lower case letters
# [0] treats the string as an array and takes the first item (= 1st letter of thsi word)
if include_pseudogenes.lower()[0] == 't': # if this 1st letter is 't'
    include_pseudogenes = True # set boolean in this script to True
elif include_pseudogenes.lower()[0] == 'f': # if this first letter is 'f'
    include_pseudogenes = False # set the boolean in this script to False
else: # if the user input is invalid,
    include_pseudogenes = False # set the boolean to false by default
print("include pseudogenes = " + str(include_pseudogenes))

if include_coverage.lower()[0] == 't': # if this 1st letter is 't'
    include_coverage = True # set boolean in this script to True
    num_pbc = len(pbc_paths) # number of pbc paths
elif include_coverage.lower()[0] == 'f': # if this first letter is 'f'
    include_coverage = False # set the boolean in this script to False
    num_pbc = 1 # number of pbc paths is set to 1 as placeholder (mock coverage is generated)
else: # if the user input is invalid,
    include_coverage = False # set the boolean to false by default
print("include coverage = " + str(include_coverage))


# ==============================================================
# =========================== CLASSES ==========================
# ==============================================================


#  ██████  ██████  ███    ██ ████████ ██  ██████
# ██      ██    ██ ████   ██    ██    ██ ██
# ██      ██    ██ ██ ██  ██    ██    ██ ██   ███
# ██      ██    ██ ██  ██ ██    ██    ██ ██    ██
#  ██████  ██████  ██   ████    ██    ██  ██████

class Contig:
    def __init__(self, name, length):
        self.name = name
        self.length = length

        self.centre_corrector = 0 # see method set_centre_corrector() for detailed explanation
        self.set_centre_corrector()


        self.per_base_coverages = [[] for pbc in range(num_pbc)] # array that holds the coverage per base in this contig
                                # note that the element at index i is the coverage at base i+1
                                # (indexing in python is 0-based but sequences are 1-based)
        self.coverage = []
        self.coverage_sd = []

        # set that will hold all positions of Ns in this contig (1-based)
        self.positions_of_Ns = None # used in --> add_base_coverage() function


        self.genes = [] # fill with gene names in order to be able to access them from all_genes

        self.gene_coverage_mean = []
        self.gene_coverage_sd = []

        self.gene_lengths_mean = None
        self.gene_lengths_sd = None

        self.gc_content = None

        self.gene_gc_mean = None
        self.gene_gc_sd = None

        self.tetranuc_freqs = None
        self.trinuc_freqs = None
        self.dinuc_freqs = None

        self.z_score_vector = None



    def get_name(self):
        return self.name

    def get_length(self):
        return self.length

    def get_length_of_covered_bases(self):
        '''Return how many bases are covered.'''
        pbc_without_nans = [[cov for cov in coverages if not np.isnan(cov)] for coverages in self.per_base_coverages]
        return [len(coverages) for coverages in pbc_without_nans]

    def get_centre(self):
        return ((self.length / 2) + 0.5) # returns the number of the base at the centre
                                        # or, at even contig length, the pos. inbetween


    def set_centre_corrector(self):
        ''' This method sets self.centre_corrector.
        It is needed for the calculations of the absolute gene positions.
        These are calculated with respect to the central base in the contig.

        For a contig of odd length, there is an actual central base, e.g.
        1-----7-----13
        The gene positions are calculated using 1&7 and 7&13, for the left and right part, resp.

        For a contig of even length, there is no actual central base, e.g.
        1----6[6.5]7----12             (here the self-centre = 6.5)
        The gene positions will be calculated using 1&6 and 7&12.
        That's why we want self-centre_corrector (which by default is set to 0)
        to be set to 0.5 (in order to add / substract it from self.centre)
        when we have a contig of even length. '''

        if (self.length % 2 == 0): # if contig is of even length
            self.centre_corrector = 0.5


    def add_gene(self, gene):
        self.genes.append(gene)

    def get_genes(self):
        return self.genes

    def get_number_of_genes(self):
        return len(self.genes)


    def set_gene_positions(self):
        '''
        Computes the absolute positions for all genes of this contig.
        Scaling to 1-----(0.5)-----0-----(0.5)-----1

        Also overwrites self.genes array to hold the genes in their order on this contig.
        '''

        # get the middle base of this contig
        contig_centre = self.get_centre()

        # list that holds tuples of (gene, gene_start_pos)
        # start_pos info will be used to sort the corresponding genes for this contig
        tmp_gene_list = []

        # for each gene name in the array that holds all genes of this contig
        for gene in self.genes: # this fills tmp_list with the tuples (gene, gene_start_pos)
            gene_start_pos = all_genes[gene].get_start_pos()
            gene_centre = all_genes[gene].get_centre()
            gene_end_pos = all_genes[gene].get_end_pos()


            # if the gene is located in the left half of the contig
            if (gene_centre <= contig_centre):
                # (for even length contigs:
                #  correct the contig reference to the base left to the centre)
                # (nothing happens for odd contig lengths)
                corr_contig_centre = contig_centre - self.centre_corrector
                # set gene position to a value between 1 (starts contig start)
                # and > 0 (starts near contig centre)
                gene_position = (gene_start_pos - corr_contig_centre) / (1 - corr_contig_centre)


#             elif (gene_centre == contig_centre):
#                 gene_position = 0
        # this has been omitted (and < contig_centre changed to <= contig_centre)
        # because a central gene that makes up most of the contig
        # would get the position 0 (even though it most likely ends near both termini)
        # --> distorted depiction of actual distance to contig ends

            # if the gene is located in the right half of the contig
            elif (gene_centre > contig_centre):
                # (for even length contigs:
                # correct the contig reference to the base right to the centre)
                # (nothing happens for odd contig lengths)
                corr_contig_centre = contig_centre + self.centre_corrector
                # set gene position to a value between > 0 (ends near the centre of the contig)
                # and 1 (ends at contig end)
                gene_position = (gene_end_pos - corr_contig_centre) / (self.length - corr_contig_centre)
                # note: self.length is equal to the number of the base at the contig end

            # if none of the above statements evaluates to true
            # either gene_centre or contig_centre is erroneous
            else:
                print("""ERROR IN GENE POS CALCULATION. \n
                gene = """, gene, ", gene_centre = ", gene_centre, ", contig_centre= ", contig_centre)


            all_genes[gene].set_absolute_pos(gene_position)

            # append gene and its start pos info to this list
            tmp_gene_list.append((gene, gene_start_pos))


        # sort list holding (gene, start_pos) by second item in the tuple, i.e. gene_start_pos
        sorted_gene_list = sorted(tmp_gene_list, key=itemgetter(1))

        # substitute contig's gene list with list of genes sorted by start_pos
        self.genes = [i[0] for i in sorted_gene_list]


    def set_positions_of_Ns(self, N_positions):
        self.positions_of_Ns = N_positions


    def add_base_coverage(self, counter, base, coverage):
        '''Adds coverage info for the base i to per_base_coverages[i-1].
        Note: indices of per_base_coverage are 0-based,
        while coverage file is 1-based.'''

        # if bases are skipped in the PBC file (because no coverage is given)
        # add dummy values to the coverage array (as to avoid a shift in positions)
        while len(self.per_base_coverages[counter]) < (base -1):
            # if the base without coverage is an N
            if base in self.positions_of_Ns:
                self.per_base_coverages[counter].append(np.nan) # add a NaN
                # --> the base (and its missing cov) will be ignored in calculating cov metrics
            else:
                # if the base is not an N
                # add a 0 --> the base will be included in the cov metrics calculation
                self.per_base_coverages[counter].append(0)

        # check whether coverage info will be inserted at the right place
        # e.g. for coverage info for base 1, we want per_base_coverages to have a length of 0
        if len(self.per_base_coverages[counter]) == (base-1):
            self.per_base_coverages[counter].append(coverage)

        else:
            print("ERROR! Unexpected base position in per base coverage info") # THROW ERROR
            print("Contig: ", self.get_name(),
                  ", Expected Base: ", (len(self.per_base_coverages[counter])+1), ", Got Base: ", base)


    def compute_own_coverage_info(self):
        # since, while filling the per_base_coverage array,
        # missing base coverages are substituted with NaNs
        # compute mean and SD while igonring NaNs
        self.coverage = [np.nanmean(coverages) if (len(coverages) > 0) else np.nan for coverages in self.per_base_coverages]
        self.coverage_sd = [np.nanstd(coverages) if (len(coverages) > 0) else np.nan for coverages in self.per_base_coverages]

    def no_coverage_info(self):
        if len(self.per_base_coverages) == 0:
            return True
        else:
            return False

    def get_coverage(self):
        return self.coverage

    def get_coverage_sd(self):
        return self.coverage_sd


    def compute_gene_coverage_info(self):
        '''Compute the coverage mean & SD and set them
        for each Gene object associated with this contig'''

        for gene_name in self.genes:
            gene = all_genes[gene_name]

            # because bp are 1-based but python lists are 0-based:
            cov_start = gene.get_start_pos() - 1

            # because list slicing ([cov_start:cov_end]) stops BEFORE the given index:
            cov_end = gene.get_end_pos()
            # (one would want to have gene_end+1, but since the array is 0-based
            # it's already "shifted" +1 to the right

            # since, while filling the per_base_coverage array,
            # missing base coverages are substituted with NaNs
            # get array WITHOUT NaNs
            gene_pbc_without_nans = []
            for pbc in range(num_pbc):
                gene_pbc_without_nans.append([cov for cov in self.per_base_coverages[pbc][cov_start:cov_end]
                                                if not np.isnan(cov)])

            gene_coverage = [np.mean(coverages) for coverages in gene_pbc_without_nans]
            gene_coverage_sd = [np.std(coverages, ddof=1) for coverages in gene_pbc_without_nans]
            # set the coverage info for the gene
            gene.set_coverage_info(gene_coverage, gene_coverage_sd)
            # also set info on how many bases are covered in this gene
            gene.set_length_of_covered_bases([len(coverages) for coverages in gene_pbc_without_nans])




    def set_gene_coverage_mean_and_sd(self, gene_cov_mean, gene_cov_sd):
        self.gene_coverage_mean = gene_cov_mean
        self.gene_coverage_sd = gene_cov_sd

    def get_gene_coverage_mean(self):
        return self.gene_coverage_mean

    def get_gene_coverage_sd(self):
        return self.gene_coverage_sd

    def get_coverage_array(self, start_pos, end_pos):
        '''Return the coverage for all bases starting with
        (including) start_pos and ending with end_pos'''
        return self.coverages[(start_pos-1), end_pos]
    # A) why "start_pos-1"? Because sequences are 1-based
    # but array indexing is 0-based
    # B) why not "end_pos-1"? Because slicing stops BEFORE the given index
    # (one would want to take gene_end+1, but since the array is 0-based
    # it's already "shifted" +1 to the right


    def set_gene_gc_mean_and_sd(self, gene_gc_mean, gene_gc_sd):
        self.gene_gc_mean = gene_gc_mean
        self.gene_gc_sd = gene_gc_sd

    def get_gene_gc_mean(self):
        return self.gene_gc_mean

    def get_gene_gc_sd(self):
        return self.gene_gc_sd

    def compute_gene_positional_info(self):
        '''Compute and set all positional information for this gene
        i.e. absolute position, no right/left neighbour info, single-gene contig info'''
        # set absolute positions of genes
        if len(self.genes) > 0 :
            self.set_gene_positions()
            # self.genes get sorted by start_pos in this step

            # flag first gene as without left neighbour
            first_gene = all_genes[self.genes[0]]
            first_gene.set_no_left_neighbour(1)

            # flag last gene as without right neighbour
            last_gene = all_genes[self.genes[(len(self.genes) - 1)]]
            last_gene.set_no_right_neighbour(1)

            # if there is only one gene in self.genes
            if len(self.genes) == 1:
                # flag as only gene
                all_genes[self.genes[0]].set_single_gene_info(1)

    def single_gene_info(self):
        '''Returns pseudobool 1 if this contig has only one gene
        (0 otherwise).'''

        if len(self.genes) == 1:
            return 1
        else:
            return 0

    def geneless_info(self):
        '''Returns pseudobool 1 if this contig has no genes
        (0 otherwise).'''

        if len(self.genes) == 0:
            return 1
        else:
            return 0



    def set_gene_lengths_mean_and_sd(self, gene_lengths_mean, gene_lengths_sd):
        self.gene_lengths_mean = gene_lengths_mean
        self.gene_lengths_sd = gene_lengths_sd

    def get_gene_lengths_mean(self):
        return self.gene_lengths_mean

    def get_gene_lengths_sd(self):
        return self.gene_lengths_sd

    def covdev_from_overall(self, mean_ref, sd_ref, counter):
        '''Indicates how much the contig cov deviates from mean contig cov,
        in units of contig cov SD (overall).'''

        # if coverage is uniform (e.g. because of mock coverage)
        # there is no SD in coverage; thus deviation in SD can not be computed
        if sd_ref == 0:
            return np.nan
        # deviation in units of 1 SD (overall)
        dev_in_sd = abs(mean_ref - self.coverage[counter]) / sd_ref

        # if own cov is smaller than overall contig cov mean
        if self.coverage[counter] < mean_ref:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd

    def set_gc_content(self, gc_content):
        self.gc_content = gc_content

    def get_gc_content(self):
        return self.gc_content

    def gcdev_from_overall(self, mean_ref, sd_ref):
        '''Indicates how much the contig GC deviates from mean contig GC,
        in units of contig GC SD (overall).'''

        # deviation in units of 1 SD (overall)
        dev_in_sd = abs(mean_ref - self.gc_content) / sd_ref

        # if own GC is smaller than overall contig GC mean
        if self.gc_content < mean_ref:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd

    def set_oligofreqs(self, tetranuc_freqs, trinuc_freqs, dinuc_freqs):
        self.tetranuc_freqs = tetranuc_freqs
        self.trinuc_freqs = trinuc_freqs
        self.dinuc_freqs = dinuc_freqs

    def get_oligofreqs(self):
        return self.tetranuc_freqs, self.trinuc_freqs, self.dinuc_freqs

    def set_z_score_vector(self, z_score_vector):
        self.z_score_vector = z_score_vector

    def get_z_score_vector(self):
        return self.z_score_vector


#  ██████  ███████ ███    ██ ███████
# ██       ██      ████   ██ ██
# ██   ███ █████   ██ ██  ██ █████
# ██    ██ ██      ██  ██ ██ ██
#  ██████  ███████ ██   ████ ███████

class Gene:

    def __init__(self, name, start_pos, end_pos, contig, source, score, strand, attributes):
        self.name = name
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.contig = contig
        self.source = source
        self.score = score
        self.strand = strand
        self.attributes = attributes

        # set length
        self.length = self.end_pos - self.start_pos + 1
        self.percentage_of_contig_length = self.compute_percentage_of_contig_length()
        self.length_of_covered_bases = []

        self.coverage = [np.nan for pbc in range(num_pbc)]
        self.coverage_sd = []


        self.absolute_pos = None
        self.no_left_neighbour = 0 # "false" (= L neighbour exists) by default
        self.no_right_neighbour = 0 # "false" (= R neighbour exists) by default

        self.single_gene = 0 # "false" by default

        self.gc_content = None

        self.tetranuc_freqs = None
        self.trinuc_freqs = None
        self.dinuc_freqs = None

        self.z_score_vector = None


    def get_name(self):     # the name equals the ID in the attributes field
        return self.name    # (the ID is required to be unique within the scope of the GFF)

    def get_start_pos(self):
        return self.start_pos

    def get_end_pos(self):
        return self.end_pos

    def get_source(self):
        return self.source

    def get_contig(self):
        return self.contig

    def get_score(self):
        return self.score

    def get_strand(self):
        return self.strand

    def get_length(self):
        return self.length

    def get_attributes(self):
        return self.attributes # list of all attributes

    def get_centre(self):
        return (self.start_pos + (self.length / 2) - 0.5)
            # works for even lengths:
            # e.g. if gene covers positions 3, 4, 5, 6 in the contig (centre = 4.5)
            # --> length=4, half length =2 --> centre = 3+ 2 -0.5 = 4.5

            # and for odd lengths:
            # e.g. if gene covers positions 3, 4, 5, 6, 7 in the contig (centre = 5)
            # --> length = 5, half length = 2.5 --> centre = 3 + 2.5 - 0.5 = 5

    def set_coverage_info(self, coverage, coverage_sd):
        self.coverage = coverage
        self.coverage_sd = coverage_sd

    def get_coverage(self):
        return self.coverage

    def get_coverage_sd(self):
        return self.coverage_sd

    def set_absolute_pos(self, position):
        self.absolute_pos = position

    def get_absolute_pos(self):
        return self.absolute_pos

    def set_no_left_neighbour(self, pseudobool):
        self.no_left_neighbour = pseudobool

    def get_no_left_neighbour(self):
        return self.no_left_neighbour

    def set_no_right_neighbour(self, pseudobool):
        self.no_right_neighbour = pseudobool

    def get_no_right_neighbour(self):
        return self.no_right_neighbour

    def get_terminal_info(self):
        if self.no_left_neighbour == 1 or self.no_right_neighbour == 1:
            return 1
        else:
            return 0

    def set_single_gene_info(self, pseudobool):
        self.single_gene = pseudobool

    def get_single_gene_info(self):
        return self.single_gene


    def set_length_of_covered_bases(self, length):
        self.length_of_covered_bases = length

    def get_length_of_covered_bases(self):
        return self.length_of_covered_bases

    def compute_percentage_of_contig_length(self):
        contig_length = all_contigs[self.contig].get_length()
        return ((self.length / contig_length) * 100)

    def get_percentage_of_contig_length(self):
        return self.percentage_of_contig_length

    def set_gc_content(self, gc_content):
        self.gc_content = gc_content

    def get_gc_content(self):
        return self.gc_content


    def lendev_from_contig(self):
        '''Indicates how much the gene length deviates from the mean gene length
        on the contig, in units of gene length SD (contig).'''

        # if this is the only gene on the contig
        if self.single_gene == 1:
            # no comparisons possible
            return np.nan


        gene_length_sd_contig = all_contigs[self.contig].get_gene_lengths_sd()
        # if all genes in the contig have the same length
        if gene_length_sd_contig == 0:
            # no deviation from contig
            return np.nan

        gene_length_mean_contig = all_contigs[self.contig].get_gene_lengths_mean()

        # deviation in units of 1 SD (contig)
        dev_in_sd = abs(gene_length_mean_contig - self.length) / gene_length_sd_contig

        # if own length is smaller than gene length mean on contig
        if self.length < gene_length_mean_contig:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd


    def lendev_from_overall(self, mean_ref, sd_ref):
        '''Indicates how much the gene length deviates from overall mean gene length,
        in units of gene length SD (overall).'''
        # deviation in units of 1 SD (overall)
        dev_in_sd = abs(mean_ref - self.length) / sd_ref

        # if own length is smaller than gene length mean
        if self.length < mean_ref:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd


    def covdev_from_contig(self):
        '''Indicates how much the gene cov deviates from the mean gene coverage
        on the contig, in units of gene cov SD (contig).'''

        # if this is the only gene on the contig
        if self.single_gene == 1:
            # no comparisons possible
            return [np.nan for pbc in range(num_pbc)]

        covdevs = []
        for pbc_index in range(num_pbc):
            gene_cov_sd_contig = all_contigs[self.contig].get_gene_coverage_sd()[pbc_index]
            if gene_cov_sd_contig == 0:
                # no deviation from contig
                covdevs.append(np.nan)
                continue
            gene_cov_mean_contig = all_contigs[self.contig].get_gene_coverage_mean()[pbc_index]


            # deviation in units of 1 SD (contig)
            dev_in_sd = abs(gene_cov_mean_contig - self.coverage[pbc_index]) / gene_cov_sd_contig

            # if own cov is smaller than gene cov mean on contig
            if self.coverage[pbc_index] < gene_cov_mean_contig:
                # return negative deviation
                covdevs.append(-dev_in_sd)
            else:
                covdevs.append(-dev_in_sd)
        return covdevs


    def covdev_from_overall(self, mean_ref, sd_ref, pbc_index):
        '''Indicates how much the gene cov deviates from overall mean gene cov,
        in units of gene cov SD (overall).'''

        # if coverage is uniform (e.g. because of mock coverage)
        # there is no SD in coverage; thus deviation in SD can not be computed
        if sd_ref == 0:
            return np.nan
        # deviation in units of 1 SD (overall)
        dev_in_sd = abs(mean_ref - self.coverage[pbc_index]) / sd_ref

        # if own cov is smaller than gene cov mean
        if self.coverage[pbc_index] < mean_ref:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd


    def gcdev_from_contig(self):
        '''Indicates how much the gene GC deviates from the mean gene GC content
        on the contig, in units of gene GC SD (contig).'''
        # if this is the only gene on the contig
        if self.single_gene == 1:
            # no comparisons possible
            return np.nan


        gene_gc_sd_contig = all_contigs[self.contig].get_gene_gc_sd()
        # if all genes in the contig have the same gc content
        if gene_gc_sd_contig == 0:
            # no deviation from contig
            return np.nan

        gene_gc_mean_contig = all_contigs[self.contig].get_gene_gc_mean()

        # deviation in units of 1 SD (contig)
        dev_in_sd = abs(gene_gc_mean_contig - self.gc_content) / gene_gc_sd_contig

        # if own GC is smaller than gene GC mean on contig
        if self.gc_content < gene_gc_mean_contig:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd

    def gcdev_from_overall(self, mean_ref, sd_ref):
        '''Indicates how much the gene GC content deviates from overall mean gene GC,
        in units of gene GC SD (overall).'''

        # deviation in units of 1 SD (overall)
        dev_in_sd = abs(mean_ref - self.gc_content) / sd_ref

        # if own GC is smaller than gene GC mean
        if self.gc_content < mean_ref:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd


    def set_oligofreqs(self, tetranuc_freqs, trinuc_freqs, dinuc_freqs):
        self.tetranuc_freqs = tetranuc_freqs
        self.trinuc_freqs = trinuc_freqs
        self.dinuc_freqs = dinuc_freqs

    def get_oligofreqs(self):
        return self.tetranuc_freqs, self.trinuc_freqs, self.dinuc_freqs

    def set_z_score_vector(self, z_score_vector):
        self.z_score_vector = z_score_vector

    def get_z_score_vector(self):
        return self.z_score_vector




# ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
# ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
# █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
# ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
# ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████

# ==============================================================
# ========================= FUNCTIONS ==========================
# ==============================================================
def compare(given, relation, crit):
    '''Returns True, if the given value relates to the critical value
    in the way that is described by relation.
    i.e. if relation=">" and crit=3, then
    True is returned e.g. if given=2, False is returned if given=4.'''
    operators = {">" : operator.gt, "above" : operator.gt,
                 "<" : operator.lt, "below" : operator.lt,
                ">=" : operator.ge, "<=" : operator.le,
                "==" : operator.eq, "eq" : operator.eq}
    return operators[relation](given, crit)


def round_if_float(x, y):
    '''Returns rounded x (if it is a float).
    x = value to be rounded
    y = digits after the decimal point'''

    if isinstance(x, float) and not np.isnan(x):
        return round(x,y)
    else:
        return x


def invert(x, min_value, max_value):
    '''Invert the given value on its respective scale,
    as defined by the minimum and the maximum value of the scale
    (if x is not NaN).'''

    if np.isnan(x):
        return x
    else:
        return ((max_value + min_value) - x)

def rescale(x, old_min, old_max, new_min, new_max, nansub):
    '''Assigns the given value (x), which is located on a scale
    ranging from old_min to old_max, to a corresponding new value y,
    which is located on a scale ranging from new_min to new_max.
    NaN values are imputed with the value given in the variable nansub.

    Note: The current implementation features a normalisation to range [0,1]
    If you wish to change this, instructions are given below.'''


    if np.isnan(x):
        y = ((nansub - old_min)/(old_max - old_min))
    else:
        y = ((x - old_min)/(old_max - old_min))
    return y


# ██ ███    ███ ██████  ██    ██ ████████ ███████      █████  ██████  ██████   █████  ██    ██
# ██ ████  ████ ██   ██ ██    ██    ██    ██          ██   ██ ██   ██ ██   ██ ██   ██  ██  ██
# ██ ██ ████ ██ ██████  ██    ██    ██    █████       ███████ ██████  ██████  ███████   ████
# ██ ██  ██  ██ ██      ██    ██    ██    ██          ██   ██ ██   ██ ██   ██ ██   ██    ██
# ██ ██      ██ ██       ██████     ██    ███████     ██   ██ ██   ██ ██   ██ ██   ██    ██

def impute_array(data_array):
    '''This function takes a data_array, i.e. an array of arrays describing a
    matrix where each row represents a gene and each column a variable.
    NaN values in a variable are imputed with the mean of this variable
    and all values are rescaled to the range [0,1].
    '''


    # define the indices of the variables to be transformed
    index_of = {
        "c_genelensd" : 6,
        "c_genecovsd" : 11, # e.g. 12th column in the matrix holds c_genecovsd
        "g_lendev_c" : 17,
        "g_covdev_c" : 24,
        "g_gcdev_c" : 31
    }

    # we're going to loop over all genes and store all of the
    # variables to be imputed
    # --> initialise these arrays
    c_genecovsds = [[] for pbc in range(num_pbc)]
    c_genelensds = []

    g_covdev_cs = [[] for pbc in range(num_pbc)]
    g_lendev_cs = []
    g_gcdev_cs = []

    # for each gene array in the given matrix:
    for gene_array in data_array:
        for pbc_index in range(num_pbc):

            # get all gene coverage SDs which are not NaNs
            if not np.isnan(gene_array[index_of["c_genecovsd"]][pbc_index]):          # (expl.: index_of["c_genecovsd"] returns 11
                # and add them to the respective array                              # gene_array[11] will give the c_genecovsd value
                c_genecovsds[pbc_index].append(gene_array[index_of["c_genecovsd"]][pbc_index]) # for the current gene)

            # get all gene deviations from contig mean cov
            if not np.isnan(gene_array[index_of["g_covdev_c"]][pbc_index]):
                g_covdev_cs[pbc_index].append(gene_array[index_of["g_covdev_c"]][pbc_index])

        # get all gene length SDs
        if not np.isnan(gene_array[index_of["c_genelensd"]]):
            c_genelensds.append(gene_array[index_of["c_genelensd"]])

        # get all gene deviations from contig mean GCss
        if not np.isnan(gene_array[index_of["g_gcdev_c"]]):
            g_gcdev_cs.append(gene_array[index_of["g_gcdev_c"]])

        # get all gene deviations from contig mean lengths
        if not np.isnan(gene_array[index_of["g_lendev_c"]]):
            g_lendev_cs.append(gene_array[index_of["g_lendev_c"]])
    # --> at the end of this for loop, we have populated the 5 arrays with the (not-NaN) values of the genes
    # thus, each array represents all (not-NaN) values in one variable


    # for each "variable"
    # get minimum and maximum of all values (needed for inversion and rescaling (see below))
    min_c_genecovsds = [min(coverages, default=np.nan) for coverages in c_genecovsds]
    max_c_genecovsds = [max(coverages, default=np.nan) for coverages in c_genecovsds]

    min_c_genelensds = min(c_genelensds, default=np.nan)
    max_c_genelensds = max(c_genelensds, default=np.nan)

    min_g_covdev_cs = [min(coverages, default=np.nan) for coverages in g_covdev_cs]
    max_g_covdev_cs = [max(coverages, default=np.nan) for coverages in g_covdev_cs]

    min_g_lendev_cs = min(g_lendev_cs, default=np.nan)
    max_g_lendev_cs = max(g_lendev_cs, default=np.nan)

    min_g_gcdev_cs = min(g_gcdev_cs, default=np.nan)
    max_g_gcdev_cs = max(g_gcdev_cs, default=np.nan)

    # get the respective mean of each variable (ignoring NaNs)
    my_mean1 = [np.nanmean(coverages) if (len(coverages) > 0) else np.nan for coverages in c_genecovsds]
    my_mean2 = np.nanmean(c_genelensds)
    my_mean3 = [np.nanmean(coverages) if len(coverages) > 0 else np.nan for coverages in g_covdev_cs]
    my_mean4 = np.nanmean(g_gcdev_cs)
    my_mean5 = np.nanmean(g_lendev_cs)

    # loop over all genes, impute NaNs with mean of variable and rescale to given range (here: [0,1])
    for gene_array in data_array:

        if include_coverage:
            for pbc_index in range(num_pbc):
                gene_array[index_of["c_genecovsd"]][pbc_index] = rescale(gene_array[index_of["c_genecovsd"]][pbc_index], \
                                min_c_genecovsds[pbc_index], max_c_genecovsds[pbc_index], 0, 1, my_mean1[pbc_index]) #1.0, 10.0

                gene_array[index_of["g_covdev_c"]][pbc_index] = rescale(gene_array[index_of["g_covdev_c"]][pbc_index], \
                                min_g_covdev_cs[pbc_index], max_g_covdev_cs[pbc_index], 0, 1, my_mean3[pbc_index]) #1.0, 10.0

        gene_array[index_of["c_genelensd"]] = \
        rescale(gene_array[index_of["c_genelensd"]], min_c_genelensds, max_c_genelensds, 0, 1, my_mean2) #1.0, 10.0

        gene_array[index_of["g_gcdev_c"]] = \
        rescale(gene_array[index_of["g_gcdev_c"]], min_g_gcdev_cs, max_g_gcdev_cs, 0, 1, my_mean4) #1.0, 10.0

        gene_array[index_of["g_lendev_c"]] = \
        rescale(gene_array[index_of["g_lendev_c"]], min_g_lendev_cs, max_g_lendev_cs, 0, 1, my_mean5) #1.0, 10.0


    return data_array


def get_lgth_corr_stats(values, lengths):
    '''Returns mean and standard deviation of the given list of values
    (e.g. coverages, GC contents, ...),
    weighted with the given list of corresponding lengths.
    * values should be a list holding the values (cov, GC, ...) of the sequences
    of interest (e.g. all genes of one contig, or all contigs in the set)
    * lengths should be a list holding the sequence lengths
    The corresponding elements should have the same indices.'''

    # if weights sum up to 0, weighted average can not be computed
    if sum(lengths) == 0:
        return np.nan, np.nan

    # compute weighted mean of values
    weighted_mean = np.average(values, weights = lengths)

    # compute weighted variance
    sum_of_squared_deviations = 0
    for (value, length) in zip(values, lengths):
        sum_of_squared_deviations += length*((value - weighted_mean)**2)
    weighted_variance = sum_of_squared_deviations / np.sum(lengths)

    # compute weighted standard deviation
    weighted_sd = weighted_variance**0.5

    return weighted_mean, weighted_sd



def get_lgth_corr_stats_weighted_by_sq_lgth(values, lengths):
    '''Returns mean and standard deviation of the given list of values
    (e.g. coverages, GC contents, ...),
    weighted by the squared corresponding lengths of the sequences.
    * values should be a list holding the values (cov, GC, ...) of the sequences
    of interest (e.g. all genes of one contig, or all contigs in the set)
    * lengths should be a list holding the sequence lengths
    The corresponding elements should have the same indices.'''


    # compute weighted mean of coverages
    squared_lengths = [length**2 for length in lengths]
    weighted_mean = np.average(values, weights = squared_lengths)

    # compute weighted variance
    sum_of_squared_deviations = 0
    for (value, length) in zip(values, squared_lengths):
        sum_of_squared_deviations += length*((value - weighted_mean)**2)
    weighted_variance = sum_of_squared_deviations / np.sum(squared_lengths)

    # compute weighted standard deviation
    weighted_sd = weighted_variance**0.5

    return weighted_mean, weighted_sd



def get_sum_of_contig_lengths(contig_list):
    '''Returns the sum of the lengths of contigs in contig_list,
    e.g. if these are all contigs originally contained in the assembly,
    then this means that the total assembly length is returned).'''

    contig_lengths = []

    for contig_name in contig_list:
        contig = all_contigs[contig_name]
        contig_lengths.append(contig.get_length())

    return sum(contig_lengths)

def get_Nx(x, contig_lengths):
    '''Returns Nx of given contig set, e.g. if x=50 the N50.'''

    # either contig_lengths is a specific list given
    # or no list is given and all contigs are input,
    # then generate list:
    if contig_lengths == "all":
        contig_lengths = []
        for contig_name, contig in all_contigs.items():
            tmp_contig_length = contig.get_length()
            contig_lengths.append(tmp_contig_length)


    # if contig_lengths is a given list,
    # create a new list from it as to not manipulate the original one
    my_contig_lengths = list(contig_lengths)

    # --------- Calculate Nx ---------
    # sort in descending order
    my_contig_lengths.sort(reverse=True)
    # get total length of contig set
    total_length = sum(my_contig_lengths)
    # determine which percentage
    percentage = x / 100
    # determine minimum length contig is required to have
    min_length = total_length * percentage

    # start summing contig lengths until X% of total assembly length is reached
    cumulative_sum = 0
    for length in my_contig_lengths:
        cumulative_sum += length
        if cumulative_sum >= min_length:
            return length



#  ██████  ██████  ██    ██         █████  ███    ██ ██████          ██████   ██████  ███████         ██ ███    ██ ███████  ██████
# ██      ██    ██ ██    ██        ██   ██ ████   ██ ██   ██         ██   ██ ██    ██ ██              ██ ████   ██ ██      ██    ██
# ██      ██    ██ ██    ██        ███████ ██ ██  ██ ██   ██         ██████  ██    ██ ███████         ██ ██ ██  ██ █████   ██    ██
# ██      ██    ██  ██  ██         ██   ██ ██  ██ ██ ██   ██         ██      ██    ██      ██         ██ ██  ██ ██ ██      ██    ██
#  ██████  ██████    ████  ███████ ██   ██ ██   ████ ██████  ███████ ██       ██████  ███████ ███████ ██ ██   ████ ██       ██████

def compute_coverage_and_positional_info():
    '''Uses per base coverage information to compute & set
    coverage mean and SD for all contigs and genes.
    Computes & sets info on gene length mean & SD.
    Computes & sets absolute positions of genes.
    Computes & sets info on gene GC mean & SD.'''


    all_genes_coverage = []

    # for ALL contigs (with coverage info)
    for contig_name, contig in all_contigs.items():

        # compute coverage info for this contig:
        contig.compute_own_coverage_info()

        #------Length, GC & Coverage Info Genes-------
        # if contig is geneless
        # set its gene coverage mean & SD to np.nan
        if contig.geneless_info() == 1:
            contig.set_gene_coverage_mean_and_sd([np.nan for coverages in range(num_pbc)], [np.nan for coverages in range(num_pbc)])
            contig.set_gene_lengths_mean_and_sd(np.nan, np.nan)

        else:
            # compute coverage info for the genes on this contig
            #  (also set length_of_covered_bases for each gene)
            contig.compute_gene_coverage_info()
            # compute info on absolute positions of genes
            contig.compute_gene_positional_info()



            # to store all ACTUAL gene lengths on this contig
            contig_gene_lengths = []
            # to store how many bases are covered in each gene on this contig
            contig_gene_covered_base_lengths = []
            # to store each gene's mean GC value
            contig_gene_gcs = []
            # to store each gene's mean coverage
            contig_gene_covs = []



            gene_names = contig.get_genes()

            for gene_name in gene_names:
                gene = all_genes[gene_name]

                contig_gene_lengths.append(gene.get_length())
                contig_gene_covered_base_lengths.append(gene.get_length_of_covered_bases())
                contig_gene_gcs.append(gene.get_gc_content())
                contig_gene_covs.append(gene.get_coverage())



            # lengths do not have to be length corrected


            # compute length corrected gene coverage mean and SD for this contig
            contig_weighted_gene_cov_mean, contig_weighted_gene_cov_sd = [], []
            for pbc_index in range(num_pbc): #for each coverage file
                cov_mean, cov_sd = get_lgth_corr_stats( \
                [coverages[pbc_index] for coverages in contig_gene_covs], \
                [coverages[pbc_index] for coverages in contig_gene_covered_base_lengths])
                contig_weighted_gene_cov_mean.append(cov_mean)
                contig_weighted_gene_cov_sd.append(cov_sd)



            # compute length corrected gene GC mean and SD for this contig
            contig_weighted_gene_gc_mean, contig_weighted_gene_gc_sd \
            = get_lgth_corr_stats(contig_gene_gcs, contig_gene_lengths)


            # if there is only one gene on this contig
            if len(gene_names) == 1:
                # set length corrected mean as computed, but SD to NaN
                contig.set_gene_lengths_mean_and_sd(np.mean(contig_gene_lengths),
                                                np.nan) # set standard deviation to NaN
                contig.set_gene_coverage_mean_and_sd(contig_weighted_gene_cov_mean,
                                                 [np.nan for pbc in range(num_pbc)]) # set standard deviation to NaN
                contig.set_gene_gc_mean_and_sd(contig_weighted_gene_gc_mean,
                                              np.nan)


            # if there is more than one gene, also set standard deviation
            else:
                # compute & set info on gene length mean and SD in contig
                contig.set_gene_lengths_mean_and_sd(np.mean(contig_gene_lengths),
                                                np.std(contig_gene_lengths))

                # set info on weighted gene cov mean and SD on contig:
                contig.set_gene_coverage_mean_and_sd(contig_weighted_gene_cov_mean,
                                                     contig_weighted_gene_cov_sd)

                # set info on weighted gene GC mean and SD on contig:
                contig.set_gene_gc_mean_and_sd(contig_weighted_gene_gc_mean,
                                              contig_weighted_gene_gc_sd)



# -------- START: FUNCTIONS FOR TETRANUCLEOTIDE Z-SCORE / CORRELATION COEFF -----


def get_observed_oligonuc_freq(given_seq, oligonuc):
    '''Given a sequence and a oligonucleotide of interest, '''
    return BioPython_Seq(given_seq).count_overlap(oligonuc)


def get_tetranuc_freqs(given_seq):
    '''Returns a dictionary mapping each of the 4^4 = 256 possible trinucleotides
    to its observed frequency in the given sequence.'''
    return {tetranuc : get_observed_oligonuc_freq(given_seq, tetranuc) for tetranuc in tetranucleotides}

def get_trinuc_freqs(given_seq):
    '''Returns a dictionary mapping each of the 4^3 = 64 possible trinucleotides
    to its observed frequency in the given sequence.'''
    return {trinuc : get_observed_oligonuc_freq(given_seq, trinuc) for trinuc in trinucleotides}


def get_dinuc_freqs(given_seq):
    '''Returns a dictionary mapping each of the 4^2 = 16 possible dinucleotides
    to its observed frequency in the given sequence.'''
    return {dinuc : get_observed_oligonuc_freq(given_seq, dinuc) for dinuc in dinucleotides}



def update_observed_oligofreqs(old_observed_tetranuc_freqs, old_observed_trinuc_freqs, old_observed_dinuc_freqs,
                            additional_tetranuc_freqs, additional_trinuc_freqs, additional_dinuc_freqs):
    '''Takes six input dictionaries (string to integer mapping),
    the first 3 mapping all possible tetranucleotides/trinucleotides/dinucleotides
    to their respective frequencies observed in the dataset parsed SO FAR (total count)
    and last 3 mapping all possible tetranucleotides/trinucleotides/dinucleotides
    to their respective frequencies observed in the added sequences (additional count).

    Returns 3 dictionaries, mapping the tetra-/tri-/dinucleotides
    to their updated frequencies (new total counts).'''

    updated_observed_tetranucs = {tetranuc : (old_observed_tetranuc_freqs[tetranuc] + additional_tetranuc_freqs[tetranuc])
                               for tetranuc in tetranucleotides}
    updated_observed_trinucs = {trinuc : (old_observed_trinuc_freqs[trinuc] + additional_trinuc_freqs[trinuc])
                               for trinuc in trinucleotides}
    updated_observed_dinucs = {dinuc : (old_observed_dinuc_freqs[dinuc] + additional_dinuc_freqs[dinuc])
                               for dinuc in dinucleotides}


    return updated_observed_tetranucs, updated_observed_trinucs, updated_observed_dinucs


def get_expected_tetranuc_freq(tetranuc, trinuc_freqs, dinuc_freqs):
    '''Returns the expected frequency for the given tetranucleotide.'''
    # get prefix of len 3 of tetranucleotide (n1n2n3)
    first3bases = tetranuc[:3]
    # and corresponding frequency
    obsfreq_first3bases = trinuc_freqs[first3bases]

    # get suffix of len 3 of tetranucleotide (n2n3n4)
    last3bases = tetranuc[-3:]
    # and corresponding frequency
    obsfreq_last3bases = trinuc_freqs[last3bases]

    # get the two bases in the middle of the tetranuc (n2n3)
    middle2bases = tetranuc[1:-1]
    # and corresponding frequency
    obsfreq_middle2bases = dinuc_freqs[middle2bases]

    # to prevent division by 0:
    if obsfreq_middle2bases == 0:
        return np.nan
    else:
        expected_tetranuc_freq = \
        (obsfreq_first3bases*obsfreq_last3bases)/obsfreq_middle2bases
        return expected_tetranuc_freq

def get_approximate_variance_tetranuc(tetranuc, expfreq_tetranuc, trinuc_freqs, dinuc_freqs):
    '''Returns an approximation of the variance of the
    observed frequency of the given tetranucleotide.'''


    # get prefix of len 3 of tetranucleotide (n1n2n3)
    first3bases = tetranuc[:3]
    # and corresponding frequency
    obsfreq_first3bases = trinuc_freqs[first3bases]

    # get suffix of len 3 of tetranucleotide (n2n3n4)
    last3bases = tetranuc[-3:]
    # and corresponding frequency
    obsfreq_last3bases = trinuc_freqs[last3bases]

    # get the two bases in the middle of the tetranuc (n2n3)
    middle2bases = tetranuc[1:-1]
    # and corresponding frequency
    obsfreq_middle2bases = dinuc_freqs[middle2bases]

    # to prevent division by 0:
    if obsfreq_middle2bases == 0:
        return np.nan
    else:
        approx_variance = (expfreq_tetranuc *
                    (((obsfreq_middle2bases - obsfreq_first3bases)
                      *(obsfreq_middle2bases - obsfreq_last3bases))
                     /(obsfreq_middle2bases**2)))

        return approx_variance


def get_z_score(obsfreq_tetranuc, expfreq_tetranuc, approxvar_tetranuc):
    '''Given the observed frequency of a given tetranucleotide in a given sequence,
    the corresponding expected frequency and its variance,
    this function returns a z-score for this specific tetranucleotide
    in the given sequence.'''

    if np.isnan(expfreq_tetranuc) or np.isnan(approxvar_tetranuc) or approxvar_tetranuc == 0:
        z_score = np.nan
    else:
        z_score = (obsfreq_tetranuc - expfreq_tetranuc)/(approxvar_tetranuc**0.5)

    return z_score



def get_all_obsfreqs(given_seq):
    '''Returns dictionaries mapping
    - all possible 256 tetranucleotides
    - all possible 64 trinucleotides
    - all possible 16 dinucleotides
     to their frequencies observed in the given sequence
     (string to integer mapping).'''
    tetranuc_freqs = get_tetranuc_freqs(given_seq)
    trinuc_freqs = get_trinuc_freqs(given_seq)
    dinuc_freqs = get_dinuc_freqs(given_seq)

    return tetranuc_freqs, trinuc_freqs, dinuc_freqs




def calculate_z_score_vector(tetranuc_freqs, trinuc_freqs, dinuc_freqs):
    '''Given all tetranucleotide-to-frequency, trinucleotide-to-frequency
    and dinucleotide-to-frequency mappings (dictionaries with key:string, value:int)
    observed in a given sequence, this function returns an array that represents
    the z-scores of all possible 256 tetranucleotides for this sequence.'''


    z_scores = []

    for tetranuc in tetranucleotides:

        obsfreq_tetranuc = tetranuc_freqs[tetranuc]

        expfreq_tetranuc = get_expected_tetranuc_freq(tetranuc, trinuc_freqs, dinuc_freqs)

        approxvar_tetranuc = get_approximate_variance_tetranuc(tetranuc, expfreq_tetranuc,
                                                            trinuc_freqs, dinuc_freqs)

        z_score = get_z_score(obsfreq_tetranuc, expfreq_tetranuc, approxvar_tetranuc)

        z_scores.append(z_score)

    return z_scores



def pearsonr_with_nans_omitted(vector1, vector2):
    '''Whenever at least one of the given z-score vectors contains a NaN at an index,
    this element is removed, as well as the corresponding element in the other vector.
    Correlation of these NaN-free vectors is computed. Returns array holding Pearson's R and p-value.'''

    nanfree_vector1 = []
    nanfree_vector2 = []

    for (i, j) in zip(vector1, vector2):
        if not (np.isnan(i) or np.isnan(j)):
            nanfree_vector1.append(i)
            nanfree_vector2.append(j)

    # pearsonr can only be computed if vectors are at least 3 elements long
    # nanfree_vectors are too short when gene has repetive sequence,
    # because the z-score vector is full of nans then
    if len(nanfree_vector1) < 2:
        return [np.nan,np.nan]
    else:
        return stats.pearsonr(nanfree_vector1, nanfree_vector2)

# -------- END: FUNCTIONS FOR TETRANUCLEOTIDE Z-SCORE / CORRELATION COEFF -------


#  ██████  ██████  ███    ███ ██████  ██    ██ ████████ ███████     ███████ ████████  █████  ████████ ███████
# ██      ██    ██ ████  ████ ██   ██ ██    ██    ██    ██          ██         ██    ██   ██    ██    ██
# ██      ██    ██ ██ ████ ██ ██████  ██    ██    ██    █████       ███████    ██    ███████    ██    ███████
# ██      ██    ██ ██  ██  ██ ██      ██    ██    ██    ██               ██    ██    ██   ██    ██         ██
#  ██████  ██████  ██      ██ ██       ██████     ██    ███████     ███████    ██    ██   ██    ██    ███████

def compute_stats(contig_list, lgth_corr_function):
    '''Computes and returns basic positional, length, coverage and GC info.
    This function MUST be run in order to obtain cov, GC & length statistics!
    (i.e. mean and SD over all contigs and genes)'''


    contig_gcs = [] # gather mean GC value for each contig
    contig_coverages = [] # gather mean coverage for each contig
    contig_lengths = [] # gather all ACTUAL contig lengths
    # gather how many bases are COVERED for each contig
    contig_covered_bases_lengths = []



    # in the considered assembly part, initialise all possible
    # tetranucleotides, trinucleotides, dinucleotides with frequency 0
    considered_tetranuc_freqs = {nuc : 0 for nuc in tetranucleotides}
    considered_trinuc_freqs = {nuc : 0 for nuc in trinucleotides}
    considered_dinuc_freqs = {nuc : 0 for nuc in dinucleotides}


    gene_gcs = [] # gather mean GC value for each gene
    gene_coverages = [] # gather mean coverage for each gene
    gene_lengths = [] # gather all ACTUAL gene lengths
    # gather how many bases are COVERED for each gene
    gene_covered_bases_lengths = []



    for contig_name in contig_list:

        contig = all_contigs[contig_name]

        #------ Gather Length, Coverage & GC Info Contig ----------------
        # get contig length + add it to list of all contig lengths
        contig_length = contig.get_length()
        contig_lengths.append(contig_length)
        # get number of covered bases for this contig + add it to the list
        contig_covered_bases_length = contig.get_length_of_covered_bases()
        contig_covered_bases_lengths.append(contig_covered_bases_length)
        # get contig cov + add it to list of all contig covs
        contig_cov = contig.get_coverage()
        contig_coverages.append(contig_cov)
        # get contig GC + add it to list of all contig GC contents
        contig_gc = contig.get_gc_content()
        contig_gcs.append(contig_gc)



        # for this contig,  initialise all possible
        # tetranucleotides, trinucleotides, dinucleotides with frequency 0
        contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs = \
        contig.get_oligofreqs()

        # update the (considered) overall frequencies with the info from this contig
        considered_tetranuc_freqs, considered_trinuc_freqs, considered_dinuc_freqs = \
        update_observed_oligofreqs(considered_tetranuc_freqs, considered_trinuc_freqs, considered_dinuc_freqs,
                                  contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)



        # ------ Gather Length, Coverage & GC Info Genes -----------------------
        for gene_name in contig.get_genes():
            gene = all_genes[gene_name]

            # add lengths of genes on this contig to all gene lengths
            gene_length = gene.get_length()
            gene_lengths.append(gene_length)
            # get number of covered bases for this gene + add it to the list
            gene_covered_bases_length = gene.get_length_of_covered_bases()
            gene_covered_bases_lengths.append(gene_covered_bases_length)
            # add covs of genes on this contig to all gene covs
            gene_cov = gene.get_coverage()
            gene_coverages.append(gene_cov)
            # add GCs of genes on this contig to all gene GCs
            gene_gc = gene.get_gc_content()
            gene_gcs.append(gene_gc)




    # --------- Compute Contig Stats ------------

    # compute min and max contig length,
    # N50, N90 and total length for all considered contigs
    considered_assembly_length = get_sum_of_contig_lengths(contig_list)
    considered_contigs_n50 = get_Nx(50, contig_lengths)
    considered_contigs_n90 = get_Nx(90, contig_lengths)
    considered_contigs_min_length = min(contig_lengths)
    considered_contigs_max_length = max(contig_lengths)


    # calculate z-score vector for complete assembly (only taking specified contigs into account)
    considered_z_score_vector = \
    calculate_z_score_vector(considered_tetranuc_freqs, considered_trinuc_freqs, considered_dinuc_freqs)


    # compute min and max contig length,
    # N50, N90 and total length for the complete assembly
    total_assembly_length = sum(all_contig_lengths)
    all_contigs_n50 = get_Nx(50, all_contig_lengths)
    all_contigs_n90 = get_Nx(90, all_contig_lengths)
    all_contigs_min_length = min(all_contig_lengths)
    all_contigs_max_length = max(all_contig_lengths)

    # compute min and max, weighted mean and SD over all contig covs
    contig_min_cov = [min([cov[i] for cov in contig_coverages]) for i in range(len(contig_coverages[0]))]
    contig_max_cov = [max([cov[i] for cov in contig_coverages]) for i in range(len(contig_coverages[0]))]
    # WEIGHTING
    all_contigs_mean_cov, all_contigs_cov_sd = [], []
    for pbc_index in range(num_pbc):
        mean_cov, cov_sd = lgth_corr_function( \
                            [cov[pbc_index] for cov in contig_coverages], \
                            [cov[pbc_index] for cov in contig_covered_bases_lengths])
        all_contigs_mean_cov.append(mean_cov)
        all_contigs_cov_sd.append(cov_sd)




    # compute min and max, weighted mean and SD over all contig GCs
    contig_min_gc = min(contig_gcs)
    contig_max_gc = max(contig_gcs)
    # WEIGHTING
    all_contigs_mean_gc, all_contigs_gc_sd \
    = lgth_corr_function(contig_gcs, contig_lengths)


    # --------- Compute Gene Stats ------------
    # compute mean, SD, min and max over all gene lengths
    all_genes_mean_length = np.mean(gene_lengths)
    all_genes_length_sd = np.std(gene_lengths)
    gene_min_length = min(gene_lengths)
    gene_max_length = max(gene_lengths)


    # compute min and max, weighted mean and SD over all gene covs
    gene_min_cov = [min([cov[i] for cov in gene_coverages]) for i in range(len(gene_coverages[0]))]
    gene_max_cov = [max([cov[i] for cov in gene_coverages]) for i in range(len(gene_coverages[0]))]
    # WEIGHTING
    all_genes_mean_cov, all_genes_cov_sd = [], []
    for pbc_index in range(num_pbc):
        mean_cov, cov_sd = lgth_corr_function( \
                            [cov[pbc_index] for cov in gene_coverages], \
                            [cov[pbc_index] for cov in gene_covered_bases_lengths])
        all_genes_mean_cov.append(mean_cov)
        all_genes_cov_sd.append(cov_sd)
    # UNWEIGHTED COV MEAN AND SD
    all_genes_mean_cov_unweighted, all_genes_cov_sd_unweighted = [], []
    for pbc_index in range(num_pbc):
        mean_cov = np.mean([cov[pbc_index] for cov in gene_coverages])
        cov_sd = np.std([cov[pbc_index] for cov in gene_coverages])
        all_genes_mean_cov_unweighted.append(mean_cov)
        all_genes_cov_sd_unweighted.append(cov_sd)

    # compute min and max, weighted man and SD over all gene GCs
    gene_min_gc = min(gene_gcs)
    gene_max_gc = max(gene_gcs)
    # WEIGHTING
    all_genes_mean_gc, all_genes_gc_sd \
    = lgth_corr_function(gene_gcs, gene_lengths)


    # ---------  OUTPUT  ------------

    output_dict = {"contig gc mean" : all_contigs_mean_gc,
        "contig gc sd" : all_contigs_gc_sd,
        "contig gc min" : contig_min_gc,
        "contig gc max" : contig_max_gc,

         "z-score vector (considered contigs)" : considered_z_score_vector,
         "z-score vector (all contigs)" : total_z_score_vector,


         "total assembly length" : total_assembly_length,
         "all contigs n50" : all_contigs_n50,
         "all contigs n90" : all_contigs_n90,
         "all contigs len min" : all_contigs_min_length,
         "all contigs len max" : all_contigs_max_length,

        "considered assembly length" : considered_assembly_length,
         "considered contigs n50" : considered_contigs_n50,
         "considered contigs n90" : considered_contigs_n90,
         "considered contigs len min" : considered_contigs_min_length,
         "considered contigs len max" : considered_contigs_max_length,

         "gene gc mean" : all_genes_mean_gc,
         "gene gc sd" : all_genes_gc_sd,
         "gene gc min" : gene_min_gc,
         "gene gc max" : gene_max_gc,

          "gene len mean" : all_genes_mean_length,
          "gene len sd" : all_genes_length_sd,
          "gene len min" : gene_min_length,
          "gene len max" : gene_max_length,

          "contig cov mean" : all_contigs_mean_cov,
          "contig cov sd" : all_contigs_cov_sd,
          "contig cov min" : contig_min_cov,
          "contig cov max" : contig_max_cov,

          "gene cov mean" : all_genes_mean_cov,
          "gene cov sd" : all_genes_cov_sd,
          "gene cov min" : gene_min_cov,
          "gene cov max" : gene_max_cov,
          "gene cov mean unweighted" : all_genes_mean_cov,
          "gene cov sd unweighted" : all_genes_cov_sd,
     }
    # for i in range(len(pbc_paths)):
    #     cov_dict = {"contig cov mean "+str(i) : all_contigs_mean_cov[i],
    #         "contig cov sd "+str(i) : all_contigs_cov_sd[i],
    #          "contig cov min "+str(i) : contig_min_cov[i],
    #          "contig cov max "+str(i) : contig_max_cov[i],
    #
    #          "gene cov mean "+str(i) : all_genes_mean_cov[i],
    #          "gene cov sd "+str(i) : all_genes_cov_sd[i],
    #          "gene cov min "+str(i) : gene_min_cov[i],
    #          "gene cov max "+str(i) : gene_max_cov[i]}
    #     output_dict.update(cov_dict)

    return output_dict




# ███████ ██ ██   ████████ ███████ ██████       ██████  ██████  ███    ██ ████████ ██     ██████  ███████ ███    ██ ███████ ███████
# ██      ██ ██      ██    ██      ██   ██     ██      ██    ██ ████   ██    ██    ██    ██       ██      ████   ██ ██      ██
# █████   ██ ██      ██    █████   ██████      ██      ██    ██ ██ ██  ██    ██ ████████ ██   ███ █████   ██ ██  ██ █████   ███████
# ██      ██ ██      ██    ██      ██   ██     ██      ██    ██ ██  ██ ██    ██ ██  ██   ██    ██ ██      ██  ██ ██ ██           ██
# ██      ██ ███████ ██    ███████ ██   ██      ██████  ██████  ██   ████    ██ ██████    ██████  ███████ ██   ████ ███████ ███████



def filter_contigs_and_genes():
    '''Stores names of all contigs without genes in a list and returns it.
    Moves all contigs without coverage info to a new dict and returns it.
    Moves all genes without coverage info to a new dict and returns it.'''

    # list of names of contigs without genes
    geneless_contigs = []

    # dictionary to which all contigs without cov info will be moved to
    contigs_without_cov = {}

    # dictionary to which all genes without cov info will be moved to
    genes_without_cov = {}


    for contig_name, contig in all_contigs.items():
        # store contig as geneless if true
        if contig.geneless_info() == 1:
            geneless_contigs.append(contig_name)

        # if no per base coverage info is available for this contig
        if contig.no_coverage_info():
            # store it to appropriate dict
            contigs_without_cov[contig_name] = contig

            # and store all its genes to an appropriate dict
            gene_name_list = contig.get_genes()
            for gene_name in gene_name_list:
                genes_without_cov[gene_name] = all_genes[gene_name]


    # delete all contigs & genes that have no cov info
    for gene_name, gene in genes_without_cov.items():
        del all_genes[gene_name]

    for contig_name, contig in contigs_without_cov.items():
        del all_contigs[contig_name]

    return geneless_contigs, contigs_without_cov, genes_without_cov




#  ██████  ███████ ████████     ██████   █████  ██     ██      █████  ██████  ██████   █████  ██    ██
# ██       ██         ██        ██   ██ ██   ██ ██     ██     ██   ██ ██   ██ ██   ██ ██   ██  ██  ██
# ██   ███ █████      ██        ██████  ███████ ██  █  ██     ███████ ██████  ██████  ███████   ████
# ██    ██ ██         ██        ██   ██ ██   ██ ██ ███ ██     ██   ██ ██   ██ ██   ██ ██   ██    ██
#  ██████  ███████    ██        ██   ██ ██   ██  ███ ███      ██   ██ ██   ██ ██   ██ ██   ██    ██

def get_raw_array(contig_name_list, stats_ref):
    '''Loops over all  contigs requested in contig_name_list
    and returns their metrics in an (n x p) "matrix" (data type: array of arrays)
    where each of the n arrays (rows) describes a gene
    on p variables (columns; i.e. array has p elements).
    output_array contains "raw" data, i.e. prior to NaN imputatioin and rescaling.'''

    output_array = []

    for contig_name in contig_name_list:

        contig = all_contigs[contig_name]

        c_ngenes = contig.get_number_of_genes()


        # only if the contig has any genes:
        if c_ngenes > 0:

            c_len = contig.get_length()
            c_pct_assembly_len = percentage_assembly_length(c_len, stats_ref["considered assembly length"])
            # to get the percentage of only the contigs considered, this line should be:
            #c_pct_assembly_len = percentage_assembly_length(c_len, stats_ref["considered assembly length"])

            # -------oligo info--------
            # compute correlation between the z-score vectors of the oligonucleotide frequencies
            # of this contig and of all contigs (only those considered)
            c_z_score_vector = contig.get_z_score_vector()
            c_pearson_test = pearsonr_with_nans_omitted(c_z_score_vector,
                                                  stats_ref["z-score vector (considered contigs)"])
            c_pearson_r = c_pearson_test[0] # get the correlation coefficient
            c_pearson_p = c_pearson_test[1] # get the p-value of the correlation
            # -------oligo info--------


            c_gc_content = contig.get_gc_content()
            c_gcdev = contig.gcdev_from_overall(stats_ref["contig gc mean"],
                                                   stats_ref["contig gc sd"])


            c_cov = contig.get_coverage()
            c_covsd = contig.get_coverage_sd()
            c_covdev = [contig.covdev_from_overall(stats_ref["contig cov mean"][pbc_index], stats_ref["contig cov sd"][pbc_index], pbc_index) for pbc_index in range(num_pbc)]
            c_genecovm = contig.get_gene_coverage_mean()
            c_genecovsd = contig.get_gene_coverage_sd()
            c_genelenm = contig.get_gene_lengths_mean()
            c_genelensd = contig.get_gene_lengths_sd()


            gene_names = contig.get_genes()

            for gene_name in gene_names:

                gene = all_genes[gene_name]


                # compute correlation between the z-score vectors of the oligonucleotide frequencies ...
                #        ... of this gene and of all contigs (only those considered)
                g_z_score_vector = gene.get_z_score_vector()
                g_pearson_test_o = pearsonr_with_nans_omitted(g_z_score_vector,
                                                      stats_ref["z-score vector (considered contigs)"])
                g_pearson_r_o = g_pearson_test_o[0] # get the correlation coefficient
                g_pearson_p_o = g_pearson_test_o[1] # get the p-value of the correlation

                #        ... of this gene and of its contig
                g_pearson_test_c = pearsonr_with_nans_omitted(g_z_score_vector, c_z_score_vector)
                g_pearson_r_c = g_pearson_test_c[0] # get the correlation coefficient
                g_pearson_p_c = g_pearson_test_c[1] # get the p-value of the correlation
                # ------------------



                g_cov_z_score_vector = []
                for pbc_index in range(num_pbc):
                    if (gene.get_coverage()[pbc_index] - stats_ref["gene cov mean unweighted"][pbc_index]) == 0 \
                    or stats_ref["gene cov sd unweighted"][pbc_index] == 0 \
                    or (gene.get_coverage()[pbc_index] == np.nan or stats_ref["gene cov mean unweighted"][pbc_index] == np.nan or stats_ref["gene cov sd unweighted"][pbc_index] == np.nan):
                        g_cov_z_score_vector.append(0)
                    else:
                        g_cov_z_score = (gene.get_coverage()[pbc_index] - stats_ref["gene cov mean unweighted"][pbc_index]) / stats_ref["gene cov sd unweighted"][pbc_index]
                        g_cov_z_score_vector.append(g_cov_z_score)

                above_2 = 0
                below_m2 = 0
                within_2 = 0
                for z_score in g_cov_z_score_vector:
                    if z_score >= 2: # coverage of gene is more than 2SDs higher than mean
                        above_2 += 1
                    elif z_score <= -2: # coverage of gene is more than 2SDs lower than mean
                        below_m2 += 1
                    else: # coverage of gene is within 2SDs of mean
                        within_2 += 1

                if above_2 == num_pbc:
                    g_cov_z_score_bool = 1 # gene coverage is higher than average in all coverage profiles
                elif below_m2 == num_pbc:
                    g_cov_z_score_bool = 1 # gene coverage is lower than average in all coverage profiles
                elif within_2 == num_pbc:
                    g_cov_z_score_bool = 1 # gene coverage is within 2SDs of mean in all coverage profiles
                else:
                    g_cov_z_score_bool = 0 # gene has differential coverage profile


                gene_array = [

                    gene_name, # index 0
                    contig_name, # index 1

                    # spatial contig variables:
                    c_ngenes, # index 2
                    c_len, # index 3
                    c_pct_assembly_len, # c_lendev # index 4
                    c_genelenm, # index 5
                    c_genelensd, # index 6

                    # coverage-related contig variables:
                    c_cov,  # index 7
                    c_covsd, # index 8
                    c_covdev, # index 9
                    c_genecovm, # index 10
                    c_genecovsd, # index 11

                    # compositional contig variables:
                    c_pearson_r, # index 12
                    c_pearson_p, # index 13
                    c_gc_content, # index 14
                    c_gcdev, # index 15


                    # spatial gene variables:
                    gene.get_length(), # index 16
                    gene.lendev_from_contig(), # index 17
                    gene.lendev_from_overall(stats_ref["gene len mean"],
                                            stats_ref["gene len sd"]), # index 18
                    gene.get_absolute_pos(), # index 19
                    gene.get_terminal_info(), # index 20
                    gene.get_single_gene_info(), # index 21

                    # coverage-related gene variables:
                    gene.get_coverage(), # index 22
                    gene.get_coverage_sd(), # index 23
                    gene.covdev_from_contig(), # index 24
                    [gene.covdev_from_overall(stats_ref["gene cov mean"][pbc_index], stats_ref["gene cov sd"][pbc_index], pbc_index) \
                        for pbc_index in range(num_pbc)], # index 25

                    g_cov_z_score_vector, # index 33
                    g_cov_z_score_bool, # index 34

                    # compositional gene variables:
                    g_pearson_r_o, # index 26
                    g_pearson_p_o, # index 27
                    g_pearson_r_c, # index 28
                    g_pearson_p_c, # index 29

                    gene.get_gc_content(), # index 30
                    gene.gcdev_from_contig(), # index 31
                    gene.gcdev_from_overall(stats_ref["gene gc mean"],
                                            stats_ref["gene gc sd"]) # index 32





                ]

                output_array.append(gene_array)

    return output_array




#  ██████  ██    ██ ████████ ██████  ██    ██ ████████     ████████  █████  ██████  ██      ███████
# ██    ██ ██    ██    ██    ██   ██ ██    ██    ██           ██    ██   ██ ██   ██ ██      ██
# ██    ██ ██    ██    ██    ██████  ██    ██    ██           ██    ███████ ██████  ██      █████
# ██    ██ ██    ██    ██    ██      ██    ██    ██           ██    ██   ██ ██   ██ ██      ██
#  ██████   ██████     ██    ██       ██████     ██           ██    ██   ██ ██████  ███████ ███████

def output_table(raw_array, filename):
    '''Prints raw_array including header to a file.'''

    path1_raw_array = output_dir + "/" + filename + ".txt"
    path2_raw_array = output_dir + "/" + filename + ".csv"
    out1_raw_array = open(path1_raw_array, "w")
    out2_raw_array = open(path2_raw_array, "w")


    c_cov_header = []
    g_cov_header = []
    for elem in ["c_cov", "c_covsd", "c_covdev", "c_genecovm", "c_genecovsd"]:
        for pbc_index in range(num_pbc):
            c_cov_header.append(elem+"_"+str(pbc_index))
    for elem in (["g_cov", "g_covsd", "g_covdev_c", "g_covdev_o", "g_cov_zscore"]):
        for pbc_index in range(num_pbc):
                g_cov_header.append(elem+"_"+str(pbc_index))


    header = [
        "g_name", "c_name",
        "c_num_of_genes", "c_len", "c_pct_assemby_len", "c_genelenm", "c_genelensd"] + \
        c_cov_header + \
        ["c_pearson_r", "c_pearson_p", "c_gc_cont", "c_gcdev",
        "g_len", "g_lendev_c", "g_lendev_o", "g_abspos", "g_terminal", "g_single"] + \
        g_cov_header + \
        ["g_cov_z_bool", "g_pearson_r_o", "g_pearson_p_o", "g_pearson_r_c", "g_pearson_p_c",
        "g_gc_cont", "g_gcdev_c", "g_gcdev_o"
    ]


    out1_raw_array.write(("\t".join(x for x in header)) + "\n")
    out2_raw_array.write((",".join(x for x in header)) + "\n")

    for gene_array in raw_array:
        flat_gene_array = []
        for item in gene_array:
            if type(item) == list:
                for subitem in item:
                    flat_gene_array.append(subitem)
            else:
                flat_gene_array.append(item)
        if len(header) != len(flat_gene_array):
            print(header)
            print(len(header), len(flat_gene_array))
            print(flat_gene_array)
            print(gene_array)
        out1_raw_array.write(("\t".join(str(round_if_float(x, 4)) for x in flat_gene_array)) + "\n")
        out2_raw_array.write((",".join(str(round_if_float(x, 4)) for x in flat_gene_array)) + "\n")

    out1_raw_array.close()
    out2_raw_array.close()




def complementary_assembly_length(c_len, total_assembly_length):
    '''Returns how much percent of the total assembly length is
    NOT made up by this contig.
    E.g.   total length 100, contig length 10 --> will return 0.9
           total length 100, contig length 95 --> will return 0.05
    By this, larger contigs will be penalised less than smaller ones.'''

    pct_of_total = c_len / total_assembly_length
    inv_pct_of_total = 1 - pct_of_total

    return inv_pct_of_total

def percentage_assembly_length(c_len, total_assembly_length):

    '''Returns the percentage of total assembly length
    covered by this contig.'''
    pct_of_total = c_len / total_assembly_length

    return pct_of_total

def deviation_from_n50(c_len, n50):
    '''Returns the contig's deviation from N50,
    normalised by N50.
    If the contig is smaller than N50 this will be a positive value between 0 and 1.
    If the contig is larger than N50 this will be a negative value.
    By adding the value (NOT the absolute value) to the sum of metrics,
    smaller contigs are penalised while larger contigs are rewarded.'''

    deviation = n50 - c_len
    normalised_deviation = deviation / n50

    return normalised_deviation






# ███████ ██    ██ ███    ███ ███    ███  █████  ██████  ██    ██     ███████ ██ ██      ███████
# ██      ██    ██ ████  ████ ████  ████ ██   ██ ██   ██  ██  ██      ██      ██ ██      ██
# ███████ ██    ██ ██ ████ ██ ██ ████ ██ ███████ ██████    ████       █████   ██ ██      █████
#      ██ ██    ██ ██  ██  ██ ██  ██  ██ ██   ██ ██   ██    ██        ██      ██ ██      ██
# ███████  ██████  ██      ██ ██      ██ ██   ██ ██   ██    ██        ██      ██ ███████ ███████

def output_summary_file(stats_ref, filename):
    path_summaryfile = output_dir + "/" + filename
    summaryfile = open(path_summaryfile, "w")

    contig_cov_summary = ""
    for pbc_index in range(num_pbc):
        contig_cov_summary += "contig cov mean\t" + str(round(stats_ref["contig cov mean"][pbc_index], 2)) + "\n" + \
         "contig cov sd\t" + str(round(stats_ref["contig cov sd"][pbc_index], 2)) + "\n" + \
         "contig cov min\t" + str(round(stats_ref["contig cov min"][pbc_index], 2)) + "\n" + \
         "contig cov max\t" + str(round(stats_ref["contig cov max"][pbc_index], 2)) + "\n"
    contig_cov_summary += "\n"

    gene_cov_summary = ""
    for pbc_index in range(num_pbc):
        gene_cov_summary += "gene cov mean\t" + str(round(stats_ref["gene cov mean"][pbc_index], 2))  + "\n" + \
          "gene cov sd\t" + str(round(stats_ref["gene cov sd"][pbc_index], 2)) + "\n" + \
          "gene cov min\t" + str(round(stats_ref["gene cov min"][pbc_index], 2))  + "\n" + \
          "gene cov max\t" + str(round(stats_ref["gene cov max"][pbc_index], 2)) + "\n"
    gene_cov_summary += "\n"


    summaryfile.write("total # of contigs:\t" + str(len(all_contigs) + len(contigs_without_cov)) + "\n" +
                      "\t# geneless contigs:\t" + str(len(geneless_contigs))  + "\n" +
                      "\t# contigs w/ genes:\t"
                          + str(len(all_contigs) + len(contigs_without_cov) - len(geneless_contigs)) + "\n" +

                      "\t# contigs w/o cov info:\t" + str(len(contigs_without_cov)) + "\n" +
                      "\t# contigs w/ cov info\t" + str(len(all_contigs)) + "\n" +

                    "\t# geneless contigs w/o cov info:\t"
                                  + str(len(set(contigs_without_cov).intersection(geneless_contigs))) + "\n" +
                    "\t# contigs with genes and w/ cov info:\t"
                              + str(len(all_contigs) - len(geneless_contigs)) + "\n\n" +


                "total assembly length\t" + str(stats_ref["total assembly length"]) + "\n" +
                 "all contigs N50\t" + str(round(stats_ref["all contigs n50"], 2)) + "\n" +
                  "all contigs N90\t" + str(round(stats_ref["all contigs n90"], 2)) + "\n" +
                  "all contigs len min\t" + str(stats_ref["all contigs len min"]) + "\n" +
                  "all contigs len max\t" + str(stats_ref["all contigs len max"]) + "\n\n" +



                "considered assembly length\t" + str(stats_ref["considered assembly length"]) + "\n" +
                "considered contigs N50\t" + str(round(stats_ref["considered contigs n50"], 2)) + "\n" +
                "considered contigs N90\t" + str(round(stats_ref["considered contigs n90"], 2)) + "\n" +
                "considered contigs len min\t" + str(stats_ref["considered contigs len min"]) + "\n" +
                "considered contigs len max\t" + str(stats_ref["considered contigs len max"]) + "\n\n" +


                contig_cov_summary +


                "contig gc mean\t" + str(round(stats_ref["contig gc mean"], 2)) + "\n" +
                 "contig gc sd\t" + str(round(stats_ref["contig gc sd"], 2)) + "\n" +
                 "contig gc min\t" + str(round(stats_ref["contig gc min"], 2)) + "\n" +
                 "contig gc max\t" + str(round(stats_ref["contig gc max"], 2)) + "\n\n" +


                "total # genes:\t" + str(len(all_genes) + len(genes_without_cov)) + "\n" +
                "\t# genes w/o cov info:\t" + str(len(genes_without_cov)) + "\n" +
                "\t# genes w/ cov info\t" + str(len(all_genes))  + "\n\n" +


                 "gene len mean\t" + str(round(stats_ref["gene len mean"], 2))  + "\n" +
                  "gene len sd\t" + str(round(stats_ref["gene len sd"], 2)) + "\n" +
                  "gene len min\t" + str(stats_ref["gene len min"]) + "\n" +
                  "gene len max\t" + str(stats_ref["gene len max"]) + "\n\n" +

                  gene_cov_summary +

                  "gene gc mean\t" + str(round(stats_ref["gene gc mean"], 2))  + "\n" +
                  "gene gc sd\t" + str(round(stats_ref["gene gc sd"], 2)) + "\n" +
                  "gene gc min\t" + str(round(stats_ref["gene gc min"], 2))  + "\n" +
                  "gene gc max\t" + str(round(stats_ref["gene gc max"], 2)) + "\n"


                 )

    summaryfile.close()


# !!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT PART !!!!!!!!!!!!!!!!!!!!!!!!!!!

# ==============================================================
# ====================== GLOBAL VARIABLES ======================
# ==============================================================

all_contigs = {} # dict to hold all contigs parsed; key=contig name, value=object contig
all_genes = {} # dict to hold all genes parsed; key=gene name, value=object gene

all_contig_lengths = []


#  ██████  ███████ ███████
# ██       ██      ██
# ██   ███ █████   █████
# ██    ██ ██      ██
#  ██████  ██      ██

# ==============================================================
# ======================= READ GFF FILE ========================
# ==============================================================
fasta_tag = False   # boolean flag used to control parsing
                    # will be set to true when ##FASTA block is reached

with open(gff_path, 'r') as gff: # open GFF file
    for line in gff:
        if not fasta_tag: # if FASTA block has not been reached yet

            if line.startswith('#'):
                if "#FASTA" in line: # if FASTA block has been reached
                    fasta_tag = True           # set boolean flag to true
                    break                      # and stop parsing


            else: # if the line holds annotations

                tmp_array = line.split()

                if tmp_array[2] == contig_tag:
                    # if tmp_array[2] == "contig":
                    contig_name = tmp_array[0]
                    contig_length = int(tmp_array[4])
                    # add this length to all_contig_lengths
                    # to store the size of the initial (unfiltered) assembly
                    all_contig_lengths.append(contig_length)

                    # initialise contig
                    contig = Contig(contig_name, contig_length)
                    # add to dictionary of contigs
                    all_contigs[contig_name] = contig

                elif (tmp_array[2] == "gene" or (tmp_array[2] == "pseudogene" and include_pseudogenes == True)):
                    associated_contig = tmp_array[0]
                    source = tmp_array[1]
                    start_pos = int(tmp_array[3])
                    end_pos = int(tmp_array[4])

                    if tmp_array[5] == '.':          # if there is no score
                        score = tmp_array[5]         # simply add the '.'
                    else:                                 # if a score is given
                        score = float(tmp_array[5])  # add it as a float

                    strand = tmp_array[6]
                    attributes = tmp_array[8].split(';')
                    id_attribute = attributes[0].split('=') # something like ID=gene1
                    gene_name = id_attribute[1] # get only the name after the equals sign

                    # initialise gene
                    gene = Gene(gene_name, start_pos, end_pos, associated_contig, source, score, strand, attributes)

                    # add to dictionary of genes
                    all_genes[gene_name] = gene

                    # add to list of genes in associated contig
                    all_contigs[associated_contig].add_gene(gene_name)

                else:
                    pass


# ███████  █████  ███████ ████████  █████
# ██      ██   ██ ██         ██    ██   ██
# █████   ███████ ███████    ██    ███████
# ██      ██   ██      ██    ██    ██   ██
# ██      ██   ██ ███████    ██    ██   ██

# ==============================================================
# ===================== READ FASTA FILE =====================
# ==============================================================


# global variable needed to keep track of currently parsed contig
current_contig = None




# ---- start: global variables needed for oligonucleotide frequencies ------
# define nucleotide alphabet
nuc_alphabet = ['A', 'C', 'G', 'T']

# generate an array holding all possible (4^4 = 256) tetranucleotides as strings
tetranucleotides = [''.join(i) for i in itertools_product(nuc_alphabet, repeat = 4)]
# generate an array holding all possible (4^3 = 64) trinucleotides as strings
trinucleotides = [''.join(i) for i in itertools_product(nuc_alphabet, repeat = 3)]
# generate an array holding all possible (4^2 = 16) dinucleotides as strings
dinucleotides = [''.join(i) for i in itertools_product(nuc_alphabet, repeat = 2)]

# intitialise dictionaries to keep count of overall oligo frequencies
total_tetranuc_freqs = {nuc : 0 for nuc in tetranucleotides}
total_trinuc_freqs = {nuc : 0 for nuc in trinucleotides}
total_dinuc_freqs = {nuc : 0 for nuc in dinucleotides}
# ---- end: global variables needed for oligonucleotide frequencies ------


def return_positions_of_Ns(sequence):
    '''For a given sequence (e.g. scaffold / contig / gene) this function
    will return a set holding all indices (1-)'''
    return {(i+1) for i, base in enumerate(sequence) if base == "N"}

# needs to be newline-free FASTA
with open(fasta_path, 'r') as fasta: # open FASTA file
    for line in fasta:
        line = line.strip()

        # if line contains fasta header
        if line.startswith(">"):
            tmp_array = line.split()
            tmp_contig_name = tmp_array[0][1:]
            if tmp_contig_name in all_contigs.keys():
                current_contig = all_contigs[tmp_contig_name]
            else: #contig in FASTA but not in GFF -> geneless, skip line with sequence
                next(fasta)
        else:
            # just in case (for case-insensitive counting)
            orig_seq = line.upper() # convert all bases to upper-case

            # get a set holding the (1-based) positions of all Ns in this contig
            contig_N_pos = return_positions_of_Ns(orig_seq)
            # pass the set to the contig
            current_contig.set_positions_of_Ns(contig_N_pos)

            contig_seq = orig_seq.replace("N", "") # remove all ambiguous characters
            # CAUTION!! LENGTH OF CONTIGS / GENES IS *WITH* AMBIGUOUS CHARACTERS

            # ----- OLIGOTEST
            # count all tetra-/tri-/dinucleotides in the contig sequence and assign to contig
            contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs = \
            get_all_obsfreqs(contig_seq)
            current_contig.set_oligofreqs(contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)

            # calculate and set the vector of the 256 oligonucleotide freq z-scores of this contig
            contig_z_score_vector =\
            calculate_z_score_vector(contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)
            current_contig.set_z_score_vector(contig_z_score_vector)

            # add oligofrequencies in contig to total oligofrequencies observed so far
            total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs = \
            update_observed_oligofreqs(total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs,
                                      contig_tetranuc_freqs, contig_trinuc_freqs, contig_dinuc_freqs)
            # ----- OLIGOTEST



            # compute GC content of this contig
            contig_gc = (contig_seq.count('G')  / len(contig_seq)) \
            + (contig_seq.count('C')  / len(contig_seq))

            current_contig.set_gc_content(contig_gc)


            # for each gene on this contig
            for gene_name in current_contig.get_genes():

                current_gene = all_genes[gene_name]

                # ---- get start and end position ------------
                # because bp are 1-based but python lists are 0-based:
                gene_seq_start = current_gene.get_start_pos() - 1

                # because list slicing ([cov_start:cov_end]) stops BEFORE the given index:
                gene_seq_end = current_gene.get_end_pos()
                # (one would want to have gene_end+1, but since the array is 0-based
                # it's already "shifted" +1 to the right
                # --------------------------------------------

                # regard only gene region on contig
                gene_seq = orig_seq[gene_seq_start:gene_seq_end].replace("N", "")

                # count all tetra-/tri-/dinucleotides in the gene sequence and assign to gene
                gene_tetranuc_freqs, gene_trinuc_freqs, gene_dinuc_freqs = \
                get_all_obsfreqs(gene_seq)
                current_gene.set_oligofreqs(gene_tetranuc_freqs, gene_trinuc_freqs, gene_dinuc_freqs)

                # calculate and set the vector of the 256 oligonucleotide freq z-scores of this gene
                gene_z_score_vector =\
                calculate_z_score_vector(gene_tetranuc_freqs, gene_trinuc_freqs, gene_dinuc_freqs)
                current_gene.set_z_score_vector(gene_z_score_vector)

                # compute GC content of current gene
                gene_gc = (gene_seq.count('G') / len(gene_seq)) \
                + (gene_seq.count('C') / len(gene_seq))

                current_gene.set_gc_content(gene_gc)




total_z_score_vector = calculate_z_score_vector(total_tetranuc_freqs, total_trinuc_freqs, total_dinuc_freqs)

# ██████  ██████   ██████
# ██   ██ ██   ██ ██
# ██████  ██████  ██
# ██      ██   ██ ██
# ██      ██████   ██████

# ==============================================================
# ===================== READ COVERAGE FILE =====================
# ==============================================================
current_contig = None # contig (object) to which coverages are currently added
current_contig_name = "dummy" # name of contig to which coverages are currently added


if include_coverage:
    for pbc_index, pbc_path in enumerate(pbc_paths):
        with open(pbc_path, 'r') as pbc: # open PBC file
            for line in pbc:

                pbc_info_array = line.split()

                if not pbc_info_array[0] == current_contig_name: # if new contig has been reached
                    current_contig_name = pbc_info_array[0]      # update name
                    if current_contig_name in all_contigs:       # if this is a contig with genes
                        current_contig = all_contigs[current_contig_name] # update current contig
                    else: # if current contig is not in GFF contig list (i.e. has no genes)
                        current_contig = None # set this to None so that all other info for this contig gets ignored
                        '''print("No coverage information for contig ", current_contig_name,
                              " stored since it has no genes.", file=open("error.txt", "a"))'''
                if current_contig is not None: # only add to contigs which exist in GFF file
                    base = int(pbc_info_array[1])
                    coverage = int(float(pbc_info_array[2]))
                    current_contig.add_base_coverage(pbc_index, base, coverage)
else: # if user does not want to include coverage information or it is not available
    mock_coverage = 1 # coverage that is set for all bases
    for contig_name, contig in all_contigs.items():
        if contig.geneless_info() == 0: # only store "info" if contig has genes
            for base in range(1, contig.get_length()+1): # set coverage for every base to mock coverage
                for pbc_index in range(num_pbc):
                    contig.add_base_coverage(pbc_index, base, mock_coverage)




# store names of geneless contigs to a list
# exclude all contigs and genes without coverage from further processing
# (i.e. remove them from all_contigs and store them to respective own dicts)
geneless_contigs, contigs_without_cov, genes_without_cov = filter_contigs_and_genes()
# AT THIS POINT, all_contigs WILL ONLY CONTAIN CONTIGS WITH COVERAGE INFO
# (BUT MAY STILL CONTAIN GENELESS CONTIGS)


# use per base coverage information to compute coverage mean & SD for all genes and contigs
# also compute absolute positions of genes on contigs
compute_coverage_and_positional_info()






# ----------------------- GET STATS -> values for reference ----------------------


# get stats taking ALL contigs and genes into account (except those without coverage)
stats_for_all_contigs = compute_stats(all_contigs.keys(), get_lgth_corr_stats)


# # get stats taking all contigs & genes into account but WEIGHTING BY SQUARED LENGTH
# stats_sq_weighting = compute_stats(all_contigs.keys(), get_lgth_corr_stats_weighted_by_sq_lgth)

# # long contigs (> 1000 or >600, resp.) only --------------------------
# contigs1000 = []
# contig1000_genes = []
# for contig_name, contig in all_contigs.items():
#     if contig.get_length() > 1000:
#         contigs1000.append(contig_name)
#         contig1000_genes.extend(contig.get_genes())

# stats_for_contigs1000 = compute_stats(contigs1000, get_lgth_corr_stats)

# ----------------------------------------------------------------------------------


output_dir = output_path + "gene_info"
pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)


raw_array = get_raw_array(all_contigs.keys(), stats_for_all_contigs)
output_table(raw_array, "raw_gene_table")
imputed_array = impute_array(raw_array)
output_table(imputed_array, "imputed_gene_table")


output_summary_file(stats_for_all_contigs, "summary.txt")
