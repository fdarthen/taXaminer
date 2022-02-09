
import numpy as np
from operator import itemgetter


class Config:
    """
    Object to hold all configuration options

    Attributes
    ----------
    fasta_path : str
      path to the FASTA file of the assembly
    gff_path: str

    output_path: str

    taxon_id: str

    read_paths: str

    bam_paths: str

    pbc_paths: str

    cov_set_exists: str

    include_coverage: str

    compute_coverage: str

    insert_size: str

    proteins_path: str

    extract_proteins: str

    assignment_mode: str

    quick_mode_search_rank: str

    quick_mode_match_rank: str

    tax_assignment_path: str

    compute_tax_assignment: str

    database_path: str

    taxon_exclude: str

    exclusion_rank: str

    update_plots: str

    num_groups_plot: str

    merging_labels: str

    output_pdf: str

    output_png: str

    include_pseudogenes: str

    gff_source: str

    input_variables: str

    perform_parallel_analysis: str

    num_pcs: str

    coverage_cutoff_mode: str

    perform_kmeans: str

    kmeans_k: str

    perform_hclust: str

    hclust_k: str

    perform_mclust: str

    mclust_k: str

    perform_dbscan: str

    dbscan_groups: str

    custom_eps: str

    custom_minPts: str

    script_dir: str

    usr_cfg_path: str

    cfg_path: str

    gff_source: str

    gff_gene_tag: str

    gff_fasta_header_type: str

    gff_fasta_header_attr: str

    gff_gene_connection: str

    gff_parent_child_types: str

    gff_parent_child_attr: str

    gff_gene_attr: str

    gff_transcript_tag: str

    gff_cds_tag: str

    """

    def __init__(self, cfg_dict):
        """
        Constructs configuration object.

        Parameters
        ----------
            cfg_dict : dict
                dictionary holding the preprocessed configuration parameters
        """
        self.fasta_path = cfg_dict.get('fasta_path')
        self.gff_path = cfg_dict.get('gff_path')
        self.output_path = cfg_dict.get('output_path')
        self.taxon_id = cfg_dict.get('taxon_id')

        self.gff_ta_path = cfg_dict.get('gff_ta_path')
        self.threads = cfg_dict.get('threads')

        self.read_paths = cfg_dict.get('read_paths')
        self.bam_paths = cfg_dict.get('bam_paths')
        self.pbc_paths = cfg_dict.get('pbc_paths')
        self.cov_set_exists = cfg_dict.get('cov_set_exists')
        self.include_coverage = cfg_dict.get('include_coverage')
        self.compute_coverage = cfg_dict.get('compute_coverage')
        self.min_insert = cfg_dict.get('min_insert')
        self.max_insert = cfg_dict.get('max_insert')
        self.read_orientation = cfg_dict.get('read_orientation')

        self.proteins_path = cfg_dict.get('proteins_path')
        self.extract_proteins = cfg_dict.get('extract_proteins')
        self.assignment_mode = cfg_dict.get('assignment_mode')
        self.quick_mode_search_rank = cfg_dict.get('quick_mode_search_rank')
        self.quick_mode_match_rank = cfg_dict.get('quick_mode_match_rank')
        self.tax_assignment_path = cfg_dict.get('tax_assignment_path')
        self.compute_tax_assignment = cfg_dict.get('compute_tax_assignment')
        self.database_path = cfg_dict.get('database_path')
        self.taxon_exclude = cfg_dict.get('taxon_exclude')
        self.exclusion_rank = cfg_dict.get('exclusion_rank')

        self.update_plots = cfg_dict.get('update_plots')
        self.num_groups_plot = cfg_dict.get('num_groups_plot')
        self.merging_labels = cfg_dict.get('merging_labels')
        self.output_pdf = cfg_dict.get('output_pdf')
        self.output_png = cfg_dict.get('output_png')

        self.include_pseudogenes = cfg_dict.get('include_pseudogenes')
        self.gff_source = cfg_dict.get('gff_source')
        self.input_variables = cfg_dict.get('input_variables')
        self.perform_parallel_analysis = cfg_dict.get('perform_parallel_analysis')
        self.num_pcs = cfg_dict.get('num_pcs')
        self.coverage_cutoff_mode = cfg_dict.get('coverage_cutoff_mode')

        self.perform_kmeans = cfg_dict.get('perform_kmeans')
        self.kmeans_k = cfg_dict.get('kmeans_k')
        self.perform_hclust = cfg_dict.get('perform_hclust')
        self.hclust_k = cfg_dict.get('hclust_k')
        self.perform_mclust = cfg_dict.get('perform_mclust')
        self.mclust_k = cfg_dict.get('mclust_k')
        self.perform_dbscan = cfg_dict.get('perform_dbscan')
        self.dbscan_groups = cfg_dict.get('dbscan_groups')
        self.custom_eps = cfg_dict.get('custom_eps')
        self.custom_minPts = cfg_dict.get('custom_minPts')

        self.script_dir = cfg_dict.get('script_dir')
        self.usr_cfg_path = cfg_dict.get('usr_cfg_path')
        self.cfg_path = cfg_dict.get('cfg_path')

        # GFF parsing rule
        self.gff_source = cfg_dict.get('gff_dict').get('source')
        self.gff_gene_tag = cfg_dict.get('gff_dict').get('gene_tag')
        self.gff_fasta_header_type = cfg_dict.get('gff_dict').get('fasta_header_type') \
                            if 'fasta_header_type' in cfg_dict.get('gff_dict').keys() \
                            else cfg_dict.get('gff_dict').get('transcript_tag')
        self.gff_fasta_header_attr = cfg_dict.get('gff_dict').get('fasta_header_attr')
        self.gff_gene_connection = cfg_dict.get('gff_dict').get('gene_connection')
        self.gff_parent_child_types = cfg_dict.get('gff_dict').get('parent_child_types')
        self.gff_parent_child_attr = cfg_dict.get('gff_dict').get('parent_child_attr')
        self.gff_gene_attr = cfg_dict.get('gff_dict').get('gene_attr')
        self.gff_transcript_tag = cfg_dict.get('gff_dict').get('transcript_tag') \
                            if 'transcript_tag' in cfg_dict.get('gff_dict').keys() \
                            else cfg_dict.get('gff_dict').get('fasta_header_type')
        self.gff_cds_tag = cfg_dict.get('gff_dict').get('cds_tag')


class Assembly:
    """
    Object holding information regarding the assembly

    Attributes
    ----------
    gff_path : str

    fasta_path : str

    pbc_paths : dict

    num_pbc : str

    output_path : str

    output_dir : str

    total_z_score_vector : str

    contigs : dict

    all_contig_lengths : list

    genes : dict

    transcripts : dict

    geneless_contigs : dict

    contigs_without_cov : dict

    genes_without_cov : dict

    """

    def __init__(self, gff_path, fasta_path, pbc_paths, output_path):
        self.gff_path = gff_path
        self.fasta_path = fasta_path
        self.pbc_paths = pbc_paths
        self.num_pbc = max(1, len(pbc_paths))
        self.output_path = output_path
        self.output_dir = self.output_path + "gene_info"

        self.total_z_score_vector = None

        self.contigs = {}
        self.all_contig_lengths = []
        self.genes = {}
        self.transcripts = {}

        self.geneless_contigs = {} # dictionary of contigs without genes
        self.contigs_without_cov = {} # dictionary to which all contigs without cov info will be moved to
        self.genes_without_cov = {} # dictionary to which all genes without cov info will be moved to
        self.partial_genes = {}

    def get_gff_path(self):
        return self.gff_path

    def get_fasta_path(self):
        return self.fasta_path

    def get_pbc_paths(self):
        return self.pbc_paths

    def get_num_pbc(self):
        return self.num_pbc

    def get_output_path(self):
        return self.output_path

    def get_output_dir(self):
        return self.output_dir

    def get_total_z_score_vector(self):
        return self.total_z_score_vector

    def set_total_z_score_vector(self, total_z_score_vector):
        self.total_z_score_vector = total_z_score_vector


    def add_contig(self, contig_name, contig):
        self.contigs[contig_name] = contig

    def remove_contig(self, contig):
        del self.contigs[contig]

    def get_contigs(self):
        return self.contigs

    def get_contig(self, contig):
        return self.contigs.get(contig)


    def add_contig_length(self, contig_length):
        self.all_contig_lengths.append(contig_length)

    def get_all_contig_lengths(self):
        return self.all_contig_lengths


    def add_gene(self, gene_name, gene):
        self.genes[gene_name] = gene

    def remove_gene(self, gene):
        del self.genes[gene]

    def get_genes(self):
        return self.genes

    def get_gene(self, gene):
        return self.genes.get(gene)


    def add_geneless_contig(self, contig_name, contig):
        self.geneless_contigs[contig_name] = contig

    def remove_geneless_contig(self, contig):
        del self.geneless_contigs[contig]

    def get_geneless_contigs(self):
        return self.geneless_contigs

    def get_geneless_contig(self, contig):
        return self.geneless_contigs[contig]


    def add_contig_without_cov(self, contig_name, contig):
        self.contigs_without_cov[contig_name] = contig

    def remove_contig_without_cov(self, contig):
        del self.contigs_without_cov[contig]

    def get_contigs_without_cov(self):
        return self.contigs_without_cov

    def get_contig_without_cov(self, contig):
        return self.contigs_without_cov[contig]


    def add_gene_without_cov(self, gene_name, gene):
        self.genes_without_cov[gene_name] = gene

    def remove_gene_without_cov(self, gene):
        del self.genes_without_cov[gene]

    def get_genes_without_cov(self):
        return self.genes_without_cov

    def get_gene_without_cov(self, gene):
        return self.genes_without_cov[gene]


    def add_partial_gene(self, gene_name, gene):
        if gene_name in self.partial_genes.keys():
            self.partial_genes[gene_name].append(gene)
        else:
            self.partial_genes[gene_name] = [gene]

    def remove_partial_gene(self, gene):
        del self.partial_genes[gene]

    def get_partial_genes(self):
        return self.partial_genes

    def get_partial_gene(self, gene):
        return self.partial_genes[gene]

########################## CONTIG #################################

class Contig:
    """
    Object to store information for each contig

    Attributes
    ----------
    name : str

    length : str

    centre_corrector : str

    per_base_coverages : list

    coverage : list

    coverage_sd : list

    positions_of_Ns : str

    genes : list

    gene_coverage_mean : list

    gene_coverage_sd : list

    gene_lengths_mean : str

    gene_lengths_sd : str

    gc_content : str

    gene_gc_mean : str

    gene_gc_sd : str

    tetranuc_freqs : str

    trinuc_freqs : str

    dinuc_freqs : str

    z_score_vector : str

    """

    def __init__(self, name, length, a):
        self.name = name
        self.length = length

        self.centre_corrector = 0 # see method set_centre_corrector() for detailed explanation
        self.set_centre_corrector()

        self.per_base_coverages = {pbc_index: [] for pbc_index in a.get_pbc_paths().keys()} # array that holds the coverage per base in this contig
                                # note that the element at index i is the coverage at base i+1
                                # (indexing in python is 0-based but sequences are 1-based)
        self.coverage = {}
        self.coverage_sd = {}

        # set that will hold all positions of Ns in this contig (1-based)
        self.positions_of_Ns = None # used in --> add_base_coverage() function


        self.genes = [] # fill with gene names in order to be able to access them from all_genes

        self.gene_coverage_mean = {}
        self.gene_coverage_sd = {}

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

    def get_length_of_covered_bases(self, pbc_index):
        '''Return how many bases are covered.'''
        pbc_without_nans = [cov for cov in self.per_base_coverages[pbc_index] if not np.isnan(cov)]
        return len(pbc_without_nans)

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


    def set_gene_positions(self, a):
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
            gene_start_pos = a.get_gene(gene).get_start_pos()
            gene_centre = a.get_gene(gene).get_centre()
            gene_end_pos = a.get_gene(gene).get_end_pos()


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


            a.get_gene(gene).set_absolute_pos(gene_position)

            # append gene and its start pos info to this list
            tmp_gene_list.append((gene, gene_start_pos))


        # sort list holding (gene, start_pos) by second item in the tuple, i.e. gene_start_pos
        sorted_gene_list = sorted(tmp_gene_list, key=itemgetter(1))

        # substitute contig's gene list with list of genes sorted by start_pos
        self.genes = [i[0] for i in sorted_gene_list]


    def set_positions_of_Ns(self, N_positions):
        self.positions_of_Ns = N_positions


    def add_base_coverage(self, pbc_index, base, coverage):
        '''Adds coverage info for the base i to per_base_coverages[i-1].
        Note: indices of per_base_coverage are 0-based,
        while coverage file is 1-based.'''

        # if bases are skipped in the PBC file (because no coverage is given)
        # add dummy values to the coverage array (as to avoid a shift in positions)
        while len(self.per_base_coverages.get(pbc_index)) < (base -1):
            # if the base without coverage is an N
            if base in self.positions_of_Ns:
                self.per_base_coverages[pbc_index].append(np.nan) # add a NaN
                # --> the base (and its missing cov) will be ignored in calculating cov metrics
            else:
                # if the base is not an N
                # add a 0 --> the base will be included in the cov metrics calculation
                self.per_base_coverages[pbc_index].append(0)

        # check whether coverage info will be inserted at the right place
        # e.g. for coverage info for base 1, we want per_base_coverages to have a length of 0
        if len(self.per_base_coverages.get(pbc_index)) == (base-1):
            self.per_base_coverages[pbc_index].append(coverage)

        else:
            print("ERROR! Unexpected base position in per base coverage info") # THROW ERROR
            print("Contig: ", self.get_name(),
                  ", Expected Base: ", (len(self.per_base_coverages.get(pbc_index))+1), ", Got Base: ", base)


    def compute_own_coverage_info(self):
        # since, while filling the per_base_coverage array,
        # missing base coverages are substituted with NaNs
        # compute mean and SD while igonring NaNs
        self.coverage = {pbc_index: (np.nanmean(coverages) if coverages else np.nan) for pbc_index, coverages in self.per_base_coverages.items()}
        self.coverage_sd = {pbc_index: (np.nanstd(coverages) if coverages else np.nan) for pbc_index, coverages in self.per_base_coverages.items()}

    def no_coverage_info(self):
        if sum([len(cov) if cov else 0 for cov in self.per_base_coverages.values()]) == 0:
            return True
        else:
            return False

    def get_coverage(self, pbc_index):
        return self.coverage.get(pbc_index)

    def get_coverage_sd(self, pbc_index):
        return self.coverage_sd.get(pbc_index)


    def compute_gene_coverage_info(self, a):
        '''Compute the coverage mean & SD and set them
        for each Gene object associated with this contig'''

        for gene_name in self.genes:
            gene = a.get_gene(gene_name)

            # because bp are 1-based but python lists are 0-based:
            cov_start = gene.get_start_pos() - 1

            # because list slicing ([cov_start:cov_end]) stops BEFORE the given index:
            cov_end = gene.get_end_pos()
            # (one would want to have gene_end+1, but since the array is 0-based
            # it's already "shifted" +1 to the right

            # since, while filling the per_base_coverage array,
            # missing base coverages are substituted with NaNs
            # get array WITHOUT NaNs
            gene_pbc_without_nans = {}
            for pbc_index in a.get_pbc_paths().keys():
                gene_pbc_without_nans[pbc_index] = [cov for cov in self.per_base_coverages[pbc_index][cov_start:cov_end]
                                                if not np.isnan(cov)]

            gene_coverage = {pbc_index: (np.mean(coverages) if coverages else np.nan) for pbc_index, coverages in gene_pbc_without_nans.items()}
            gene_coverage_sd = {pbc_index: (np.std(coverages, ddof=1) if coverages else np.nan) for pbc_index, coverages in gene_pbc_without_nans.items()}
            # set the coverage info for the gene
            gene.set_coverage_info(gene_coverage, gene_coverage_sd)
            # also set info on how many bases are covered in this gene
            gene.set_length_of_covered_bases({pbc_index: len(coverages) for pbc_index, coverages in gene_pbc_without_nans.items()})




    def set_gene_coverage_mean_and_sd(self, gene_cov_mean, gene_cov_sd):
        self.gene_coverage_mean = gene_cov_mean
        self.gene_coverage_sd = gene_cov_sd

    def get_gene_coverage_mean(self, pbc_index):
        return self.gene_coverage_mean[pbc_index]

    def get_gene_coverage_sd(self, pbc_index):
        return self.gene_coverage_sd[pbc_index]

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

    def compute_gene_positional_info(self, a):
        '''Compute and set all positional information for this gene
        i.e. absolute position, no right/left neighbour info, single-gene contig info'''
        # set absolute positions of genes
        if self.genes: # check for emtpy list
            self.set_gene_positions(a)
            # self.genes get sorted by start_pos in this step

            # flag first gene as without left neighbour
            first_gene = a.get_gene(self.genes[0])
            first_gene.set_no_left_neighbour(1)

            # flag last gene as without right neighbour
            last_gene = a.get_gene(self.genes[(len(self.genes) - 1)])
            last_gene.set_no_right_neighbour(1)

            # if there is only one gene in self.genes
            if len(self.genes) == 1:
                # flag as only gene
                a.get_gene(self.genes[0]).set_single_gene_info(1)

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

    def covdev_from_overall(self, mean_ref, sd_ref, pbc_index):
        '''Indicates how much the contig cov deviates from mean contig cov,
        in units of contig cov SD (overall).'''

        # if coverage is uniform (e.g. because of mock coverage)
        # there is no SD in coverage; thus deviation in SD can not be computed
        if sd_ref == 0:
            return np.nan
        # deviation in units of 1 SD (overall)
        dev_in_sd = abs(mean_ref - self.coverage.get(pbc_index)) / sd_ref

        # if own cov is smaller than overall contig cov mean
        if self.coverage.get(pbc_index) < mean_ref:
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

############################ GENE ###################################

class Gene:
    """
    Object to store information for each gene

    Attributes
    ----------
    name : str
    start_pos : str
    end_pos : str
    contig : str
    source : str
    score : str
    strand : str
    attributes : str

    # set length
    length : int
    percentage_of_contig_length : str
    length_of_covered_bases : str

    coverage : str
    coverage_sd : str


    absolute_pos : str
    no_left_neighbour : str
    no_right_neighbour : str

    single_gene : str

    gc_content : str

    tetranuc_freqs : str
    trinuc_freqs : str
    dinuc_freqs : str

    z_score_vector : str

    """

    def __init__(self, name, start_pos, end_pos, contig, source, score, strand, attributes, a):
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
        self.percentage_of_contig_length = self.compute_percentage_of_contig_length(a)
        self.length_of_covered_bases = {}

        self.coverage = {pbc_index: np.nan for pbc_index in a.get_pbc_paths().keys()} # list of mean coverage for each cov profile
        self.coverage_sd = {} # dict of coverage SD for each coverage profile


        self.absolute_pos = None
        self.no_left_neighbour = 0      # "false" (= L neighbour exists) by default
        self.no_right_neighbour = 0     # "false" (= R neighbour exists) by default

        self.single_gene = 0            # "false" by default

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

    def get_coverage(self, pbc_index):
        return self.coverage.get(pbc_index)

    def get_coverage_sd(self, pbc_index):
        return self.coverage_sd.get(pbc_index)

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

    def get_length_of_covered_bases(self, pbc_index):
        return self.length_of_covered_bases.get(pbc_index)

    def compute_percentage_of_contig_length(self, a):
        contig_length = a.get_contig(self.contig).get_length()
        return ((self.length / contig_length) * 100)

    def get_percentage_of_contig_length(self):
        return self.percentage_of_contig_length

    def set_gc_content(self, gc_content):
        self.gc_content = gc_content

    def get_gc_content(self):
        return self.gc_content


    def lendev_from_contig(self, a):
        '''Indicates how much the gene length deviates from the mean gene length
        on the contig, in units of gene length SD (contig).'''

        # if this is the only gene on the contig
        if self.single_gene == 1:
            # no comparisons possible
            return np.nan


        gene_length_sd_contig = a.get_contig(self.contig).get_gene_lengths_sd()
        # if all genes in the contig have the same length
        if gene_length_sd_contig == 0:
            # no deviation from contig
            return np.nan

        gene_length_mean_contig = a.get_contig(self.contig).get_gene_lengths_mean()

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


    def covdev_from_contig(self, a, pbc_index):
        '''Indicates how much the gene cov deviates from the mean gene coverage
        on the contig, in units of gene cov SD (contig).'''

        # if this is the only gene on the contig
        if self.single_gene == 1:
            # no comparisons possible
            return np.nan

        gene_cov_sd_contig = a.get_contig(self.contig).get_gene_coverage_sd(pbc_index)
        if gene_cov_sd_contig == 0:
            # no deviation from contig
            return np.nan
        gene_cov_mean_contig = a.get_contig(self.contig).get_gene_coverage_mean(pbc_index)


        # deviation in units of 1 SD (contig)
        dev_in_sd = abs(gene_cov_mean_contig - self.coverage[pbc_index]) / gene_cov_sd_contig

        # if own cov is smaller than gene cov mean on contig
        if self.coverage[pbc_index] < gene_cov_mean_contig:
            # return negative deviation
            return -dev_in_sd
        else:
            return dev_in_sd


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


    def gcdev_from_contig(self, a):
        '''Indicates how much the gene GC deviates from the mean gene GC content
        on the contig, in units of gene GC SD (contig).'''
        # if this is the only gene on the contig
        if self.single_gene == 1:
            # no comparisons possible
            return np.nan


        gene_gc_sd_contig = a.get_contig(self.contig).get_gene_gc_sd()
        # if all genes in the contig have the same gc content
        if gene_gc_sd_contig == 0:
            # no deviation from contig
            return np.nan

        gene_gc_mean_contig = a.get_contig(self.contig).get_gene_gc_mean()

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
