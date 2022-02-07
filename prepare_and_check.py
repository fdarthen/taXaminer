"""Prepare config parameters and check data validity.

Processes the users config file, sets default values for those not
given by the user, checks file existences to assess what needs to be
computed and makes sanity checks for input data

Expects path to config file
"""

import pathlib # to create directories
import yaml # read config file
import sys # parse command line arguments
import logging

from classes import Config


def check_assembly_ids(config_obj):
    """Check if headers in FASTA of assembly match IDs in GFF.

    Args:
      config_obj(obj): Config class object of config parameters
    """

    fasta_ids = set()
    gff_ids = set()

    with open(config_obj.output_path+'tmp/tmp.MILTS.fasta.fai', 'r') as gfai:
        # collect IDs of scaffolds from FASTA
        for line in gfai:
            fasta_ids.add(line.split()[0])

    with open(config_obj.gff_path, 'r') as gff:
        # collect IDs of scaffolds from GFF
        for line in gff:
            if not line.startswith('#'):
                gff_ids.add(line.split()[0])

    # match the two sets of IDs and save which could not be matched
    missing = set()
    for f_id in fasta_ids:
        if not f_id in gff_ids:
            missing.add(f_id)

    # geneless scaffolds may not be annotated in GFF
    # thus check if sets of IDs in FASTA and missing are unequal
    if len(fasta_ids) == len(missing):
        # if the sets are equally long, the IDs could not be matched at all
        sys.exit('Error: FASTA headers and scaffold IDs in GFF do not match')


def set_gff_parsing_rules(gff_source):
    """Parsing which features in GFF to read.

    Set parsing information based on default presets or user input.
    Defines which features are considered as genes, which as
    transcript, where header for FASTA can be found and how gene and
    transcript can be linked

    Args:
      gff_source(str): either name of a preset or path to file with parsing rule

    Returns:
      (dict) dictionary with parameters for parsing the GFF
    """

    # TODO: switch path with parsing to additional parameters in config file
    ## presets for parsing of GFF
    if gff_source == "default":
        rule_dict = {'source': 'default',
            'gene_tag': 'gene',
            'fasta_header_type': 'mRNA,CDS',
            'fasta_header_attr': 'ID',
            'gene_connection': 'parent/child',
            'parent_child_types': 'mRNA,CDS',
            'parent_child_attr': {'parent': 'Parent', 'child': 'ID'},
            'transcript_tag': 'mRNA',
            'cds_tag': 'CDS'}
    elif gff_source == "maker":
        rule_dict = {'source': 'maker',
            'gene_tag': 'gene',
            'fasta_header_type': 'mRNA',
            'fasta_header_attr': 'ID',
            'gene_connection': 'parent/child',
            'parent_child_types': 'mRNA',
            'parent_child_attr': {'parent': 'Parent', 'child': 'ID'},
            'transcript_tag': 'mRNA',
            'cds_tag': 'CDS'}
    elif gff_source == "augustus_masked":
        rule_dict = {'source': 'augustus_masked',
            'gene_tag': 'match',
            'fasta_header_type': 'match',
            'fasta_header_attr': 'Name',
            'gene_connection': 'inline',
            'gene_attr': 'ID',
            'cds_tag': 'match_part'}

    # read file with information on how to parse GFF
    else:
        rule_file_path = pathlib.Path(gff_source)
        if rule_file_path.is_file():
            rule_dict = {}
            with open(rule_file_path, 'r') as rule_file:
                # create dict out of lines
                for line in rule_file:
                    rule_dict[line.split(":")[0].strip()] = ":".join(line.split(":")[1:]).strip()
            # parse the dict with child and parent attributes
            if "parent_child_attr" in rule_dict.keys():
                pc_dict = {}
                for kv_pair in rule_dict.get("parent_child_attr").strip("{}").split(","):
                    pc_dict[kv_pair.split(":")[0].strip()] = kv_pair.split(":")[1].strip()
                rule_dict["parent_child_attr"] = pc_dict
        else:
            logging.error("source type for GFF could not be interpreted. Please check your input.")
            sys.exit()

    return rule_dict


def enumerated_key(config_obj, key_name, pre_keys, *default):
    """Get matches for key_name in user config file.

    Retrieve all matches for the enumerated key_name parameters (they
    have a numbered suffix) and set the paths, either as given by user
    or default values. pre_keys is a list of indexes that were
    identified for previous key_name and need to get default values if
    not defined by the user [Is used for retrieving the coverage
    information, thus if there is a set of read data with suffix 1 and
    the user does not define paths for BAM and PBC, default locations
    for these need to be set]

    Args:
      config_obj(dict): dictionary of config parameters from user file
      key_name(str): name of the parameter to be matched
      pre_keys(list): keys for which values need to be defined
      *default(str): default values for the key

    Returns:
      (dict) {enumerated suffix (index) : value}
    """

    # find matches for key
    matches = []
    for key in config_obj.keys():
        if key_name in key:
            matches.append(key)

    # check if there is a digit in the matches
    # to identify the index of the coverage set
    dict = {}
    for match in matches:
        if config_obj.get(match):
            if match.split('_')[2].isdigit():
                match_num = int(match.split('_')[2])
                dict[match_num] = set_variable_default(config_obj, match, default)

    for pre_key in pre_keys:
        # if there are more values required than given,
        # fill them up with default paths
        if not dict.get(pre_key):
            # if the default value is a path ('/' in there) then alter the
            # path to append the index of the key_name to it
            if '/' in default[0]:
                dict[pre_key] = '.'.join(default[0].split('.')[:-1])+'_'+str(pre_key)+'.'+default[0].split('.')[-1]
            else:
                dict[pre_key] = default[0]

    return dict


def check_dir_slash(path):
    """Appends trailing slash to path if necessary.

    Args:
      path(str): path to be checked for slash

    Returns:
      (str) path with slash at the end

    """
    if path.endswith('/'):
        return path
    else:
        return path+'/'


def check_file_inexistence(paths, and_or):
    """Check if file(s) at (list of) path(s) do not exists.

    Checks for a list of paths or a single one if these files exist. Returns
    'FALSE' if they do exist, else 'TRUE'. The parameter and_or decides if
    all ('AND') or at least one ('OR') of the need to exist

    Args:
      paths(str,list): path or list of paths to check file existence
      and_or(str): 'AND' if all files must exist, 'OR' if only one has to

    Returns:
        (str) 'FALSE' if the required files do exist, else 'TRUE'
    """

    if type(paths) != list:
        paths = [paths]

    counter_e = 0 # counter of how many files exist
    for path in paths:
        file = pathlib.Path(path)
        if file.is_file():
            counter_e += 1

    # check if all or at least one of the files exist
    if (and_or == 'AND' and counter_e == len(paths)) or \
     (and_or == 'OR' and counter_e > 0):
        # returns 'FALSE' if they do as this means no computations
        # need to be performed
        return 'FALSE'
    else:
        return 'TRUE'


def set_variable_default(config_obj, key_name, *default):
    """Check if parameter is given in config file, else set to default.

    Args:
      config_obj(dict): dict of user input for config parameters
      key_name(str): key to look for in config_obj
      *default(str): default parameter for key_name

    Returns:
        (str) if key_name in config_obj then value, else default

    """

    if config_obj.get(key_name):
        return config_obj[key_name]
    elif default: # if default value is given this is evaluated to True
        return default[0]
    else: # key not given in config file and no given default value
        print('Error: no default value available for "{}".\
               Please state specifically.'.format(key_name))


def pca_cov_variables(vars, include, cov_set_num):
    """Define coverage variables to be used in PCA.

    Remove coverage variables from the variables used in the PCA if
    coverage_include == FALSE; else append indices for the different
    coverage sets to the cov_vars

    Args:
      vars(str): comma-separated list of variables to be used in PCA
      include(bool): bool if coverage data is included or not
      cov_set_num: number of coverage sets

    Returns:
        (str) modified list of variables to be used in PCA
    """
    # TODO: same indicies used for cov sets by user and in output

    cov_vars = "c_cov,c_covsd,g_cov,g_covsd,g_covdev_c"
    out_vars = ""

    #TODO: test string or bool?
    # print(include)
    # print(bool(include))

    if include:
        for var in vars.split(','):
            if var in cov_vars:
                # add coverage variable with numbered suffix for each
                # coverage data set
                for i in range(cov_set_num):
                    out_vars = out_vars + ',' + var + '_' + str(i)
            else:
                out_vars = out_vars + ',' + var
    else:
        for var in vars.split(','):
            if var not in cov_vars:
                out_vars = out_vars + ',' + var

    return out_vars[1:]


def set_config_defaults(config_obj):
    """Set config parameters based on user input and default values.

    Check user input config file and set default values for missing ones;
    write results to temporary config file

    Args:
      config_obj(dict): dict of user input for config parameters

    Returns:
      (dict) dict of processed config parameters with default values
    """

    config_vars = {}

    default_pca_vars = "c_name,c_num_of_genes,c_len,c_genelenm,c_genelensd,\
                        g_len,g_lendev_c,g_abspos,g_terminal,c_pearson_r,\
                        g_pearson_r_o,g_pearson_r_c,c_cov,c_covsd,g_cov,\
                        g_covsd,g_covdev_c"

    ### check user input config file and set default values for missing ones
    ## General options
    config_vars['fasta_path'] = set_variable_default(config_obj, 'fasta_path')
    config_vars['gff_path'] = set_variable_default(config_obj, 'gff_path')
    config_vars['output_path'] = check_dir_slash(set_variable_default(config_obj, 'output_path'))
    config_vars['taxon_id'] = set_variable_default(config_obj, 'taxon_id')

    config_vars['threads'] = set_variable_default(config_obj, 'threads', 'auto')

    ## Coverage
    # check user input for coverage data and fill missing with default paths
    config_vars['read_paths'] = enumerated_key(config_obj, 'read_paths', [])
    config_vars['bam_paths'] = enumerated_key(config_obj, 'bam_path', list(config_vars.get('read_paths').keys()), config_vars.get('output_path')+'mapping_sorted.bam')
    config_vars['pbc_paths'] = enumerated_key(config_obj, 'pbc_path', list(config_vars.get('bam_paths').keys()), config_vars.get('output_path')+'pbc.txt')
    # check file existence for each coverage set to determine at which
    # step of the coverage computation the set is
    # dict stores which files to use for next step of computation of PBC
    cov_set_exists = {} # {index: file type to use for futher computation}
    for key in config_vars.get('pbc_paths').keys():
        # reads are required for computation of BAM, BAM for PBC
        # the last file in the list that exists is stored in the dict
        for cov_info in ['read_paths', 'bam_paths', 'pbc_paths']:
            path = config_vars.get(cov_info).get(key)
            if path:
                if check_file_inexistence(path, 'AND') == 'FALSE':
                    # file exists
                    cov_set_exists[key] = cov_info
    config_vars['cov_set_exists'] = cov_set_exists
    # no coverage info available -> no computations
    if len(config_vars.get('cov_set_exists')) == 0:
        config_vars['include_coverage'] = set_variable_default(config_obj, 'include_coverage', 'FALSE')
    else:
        config_vars['include_coverage'] = set_variable_default(config_obj, 'include_coverage', 'TRUE')
    # if for every coverage set the PBC files exist -> no computations
    if list(config_vars.get('cov_set_exists').values()).count('pbc_paths') == len(config_vars.get('cov_set_exists')) or \
        config_vars.get('include_coverage') == 'FALSE':
        config_vars['compute_coverage'] = set_variable_default(config_obj, 'compute_coverage', 'FALSE')
    else:
        config_vars['compute_coverage'] = set_variable_default(config_obj, 'compute_coverage', 'TRUE')
    config_vars['min_insert'] = enumerated_key(config_obj, 'min_insert', list(config_vars.get('read_paths').keys()), '0')
    config_vars['max_insert'] = enumerated_key(config_obj, 'max_insert', list(config_vars.get('read_paths').keys()), '500')
    config_vars['read_orientation'] = enumerated_key(config_obj, 'read_orientation', list(config_vars.get('read_paths').keys()), 'fr')


#Also, this protocol yields pairs where the expected genomic distance from end to end is about 200-500 base pairs.


    ## Taxonomic assignment
    config_vars['proteins_path'] = set_variable_default(config_obj, 'proteins_path', config_vars.get('output_path')+'proteins.faa')
    config_vars['extract_proteins'] = set_variable_default(config_obj, 'extract_proteins', check_file_inexistence(config_vars.get('proteins_path'), 'AND'))
    config_vars['assignment_mode'] = set_variable_default(config_obj, 'assignment_mode', 'exhaustive')
    # make settings for quick taxonomic assignment mode
    if config_vars.get('assignment_mode') == 'quick':
        config_vars['quick_mode_search_rank'] = set_variable_default(config_obj, 'quick_mode_search_rank', 'kingdom')
        config_vars['quick_mode_match_rank'] = set_variable_default(config_obj, 'quick_mode_match_rank', 'order')
    else:
        config_vars['quick_mode_search_rank'] = set_variable_default(config_obj, 'quick_mode_search_rank', None)
        config_vars['quick_mode_match_rank'] = set_variable_default(config_obj, 'quick_mode_match_rank', None)
    config_vars['tax_assignment_path'] = set_variable_default(config_obj, 'tax_assignment_path', config_vars.get('output_path')+'taxonomic_hits.txt')
    if config_vars.get('assignment_mode') == 'quick':
        if type(config_vars.get('tax_assignment_path')) != list:
            if len([config_vars.get('tax_assignment_path')]) != 2:
                path_1 = '.'.join(config_vars.get('tax_assignment_path').split('.')[:-1]) + "_1.txt"
                path_2 = '.'.join(config_vars.get('tax_assignment_path').split('.')[:-1]) + "_2.txt"
                config_vars['tax_assignment_path'] = [path_1, path_2]
    config_vars['compute_tax_assignment'] = set_variable_default(config_obj, 'compute_tax_assignment', check_file_inexistence(config_vars.get('tax_assignment_path'), 'AND'))
    config_vars['database_path'] = set_variable_default(config_obj, 'database_path')
    config_vars['taxon_exclude'] = set_variable_default(config_obj, 'taxon_exclude', 'TRUE')
    config_vars['exclusion_rank'] = set_variable_default(config_obj, 'exclusion_rank', 'species')

    ## Plotting
    config_vars['update_plots'] = set_variable_default(config_obj, 'update_plots', 'FALSE')
    config_vars['num_groups_plot'] = set_variable_default(config_obj, 'num_groups_plot', '25')
    config_vars['merging_labels'] = set_variable_default(config_obj, 'merging_labels', 'None')
    config_vars['output_pdf'] = set_variable_default(config_obj, 'output_pdf', 'TRUE')
    config_vars['output_png'] = set_variable_default(config_obj, 'output_png', 'FALSE')

    ## Gene info
    config_vars['include_pseudogenes'] = set_variable_default(config_obj, 'include_pseudogenes', 'FALSE')
    config_vars['gff_source'] = set_variable_default(config_obj, 'gff_source', 'default')
    config_vars['gff_dict'] = set_gff_parsing_rules(config_vars.get('gff_source'))

    ## PCA
    config_vars['input_variables'] = pca_cov_variables(set_variable_default(config_obj, 'input_variables', default_pca_vars), config_vars.get('include_coverage'), len(config_vars.get('cov_set_exists')))
    config_vars['perform_parallel_analysis'] = set_variable_default(config_obj, 'perform_parallel_analysis', 'FALSE')
    config_vars['num_pcs'] = set_variable_default(config_obj, 'num_pcs', '3')
    config_vars['coverage_cutoff_mode'] = set_variable_default(config_obj, 'coverage_cutoff_mode', 'default')

    ## Clustering
    config_vars['perform_kmeans'] = set_variable_default(config_obj, 'perform_kmeans', 'FALSE')
    config_vars['kmeans_k'] = set_variable_default(config_obj, 'kmeans_k', 'default')
    config_vars['perform_hclust'] = set_variable_default(config_obj, 'perform_hclust', 'FALSE')
    config_vars['hclust_k'] = set_variable_default(config_obj, 'hclust_k', 'default')
    config_vars['perform_mclust'] = set_variable_default(config_obj, 'perform_mclust', 'FALSE')
    config_vars['mclust_k'] = set_variable_default(config_obj, 'mclust_k', 'default')
    config_vars['perform_dbscan'] = set_variable_default(config_obj, 'perform_dbscan', 'FALSE')
    config_vars['dbscan_groups'] = set_variable_default(config_obj, 'dbscan_groups', 'default')
    config_vars['custom_eps'] = set_variable_default(config_obj, 'custom_eps', '0.3')
    config_vars['custom_minPts'] = set_variable_default(config_obj, 'custom_minPts', '10')

    return config_vars


def write_cfg2file(config_vars):
    """Write processed config information to yaml file in tmp dir.

    Args:
      config_vars(dict): processed config {parameter name : value}
    """

    with open(config_vars.get('output_path')+'tmp/tmp.cfg.yml', 'w') as out_cfg:
        for key, value in config_vars.items():
                out_cfg.write('{}: {}\n'.format(key,value))


def write_run_overview(config_path, config_vars):
    """Print relevant config information to console.

    Args:
      config_path(str): path to user config file
      config_vars(dict): processed config {parameter name : value}
    """

    logging.info('')
    logging.info("Config:\t{}".format(config_path))
    logging.info("FASTA:\t{}".format(config_vars.get('fasta_path')))
    logging.info("GFF:\t{}".format(config_vars.get('gff_path')))
    logging.info("Taxon ID:\t{}".format(config_vars.get('taxon_id')))
    logging.info("Output:\t{}\n".format(config_vars.get('output_path')))

    if config_vars.get('include_coverage') == 'TRUE':
        for cov_set, pbc_path in config_vars.get('pbc_paths').items():
            if config_vars.get('cov_set_exists').get(cov_set) == 'pbc_paths':
                logging.info("PBC {}:\t{} [exists]".format(cov_set, pbc_path))
            else:
                if config_vars.get('cov_set_exists').get(cov_set) == 'bam_paths':
                    basis = 'BAM'
                else:
                    basis = 'read mapping'
                logging.info("PBC {}:\t{}".format(cov_set, pbc_path))
                logging.info("  computation based on:\t{} [{}]".format( config_vars.get(config_vars.get('cov_set_exists').get(cov_set)).get(cov_set), basis))
        logging.info('')
    else:
        logging.info("no valid coverage information provided\n")

    if config_vars.get('extract_proteins') == 'FALSE':
        logging.info("Proteins:\t{} [exists]".format(config_vars.get('proteins_path')))
    else:
        logging.info("Proteins:\t{}".format(config_vars.get('proteins_path')))
    if config_vars.get('assignment_mode') == 'quick':
        logging.info("Quick assignment mode selected")
        logging.info("Filtering search performed on level {}".format(config_vars.get('quick_mode_search_rank')))
        logging.info("Hits accepted on level {}".format(config_vars.get('quick_mode_match_rank')))
        if config_vars.get('compute_tax_assignment') == 'FALSE':
            logging.info("Taxonomic hits files [exist]:\n{}{}".format(config_vars.get('tax_assignment_path')[0],config_vars.get('tax_assignment_path')[1]))
        else:
            logging.info("Taxonomic hits files:\n{}{}".format(config_vars.get('tax_assignment_path')[0],config_vars.get('tax_assignment_path')[1]))


    else:
        logging.info("Exhaustive assignment mode selected")
        if config_vars.get('compute_tax_assignment') == 'FALSE':
            logging.info("Taxonomic hits file:{} [exists]".format(config_vars.get('tax_assignment_path')))
        else:
            logging.info("Taxonomic hits file:{}".format(config_vars.get('tax_assignment_path')))

    logging.info("Query taxon excluded:\t{}\n".format(config_vars.get('taxon_exclude')))

    logging.info("Pseudogenes included:\t{}".format(config_vars.get('include_pseudogenes')))
    if config_vars.get('gff_source') != 'default':
        logging.info("Rule for GFF parsing:\t{}".format(config_vars.get('gff_source')))
    logging.info('\n')

    logging.info("PCA variables:\t{}".format(config_vars.get('input_variables')))
    logging.info("Parallel analysis performed:\t{}".format(config_vars.get('perform_parallel_analysis')))
    if config_vars.get('coverage_cutoff_mode') != 'default':
        logging.info("Coverage cutoff:\t{}".format(config_vars.get('coverage_cutoff_mode')))

    for clustering in ['kmeans', 'hclust', 'mclust']:
        if config_vars.get('perform_'+clustering) == 'TRUE':
            logging.info("{} clustering performed with number of groups:\t{}".format(clustering, config_vars.get(clustering+'k')))
    if config_vars.get('perform_dbscan') == 'TRUE':
        logging.info("DBSCAN clustering performed with settings:\t{}".format(config_vars.get('dbscan_groups')))
    logging.info('')


def process_config(config_path, script_dir):
    """Process the users config file and set defaults.

    Processing of user config file, setting defaults, printing overview
    and write processed config parameters to temporary config file

    Args:
      config_path(str): path to user config file
      script_dir(str): path of scripts directory
    """
    # read parameters from config file
    config_obj = yaml.safe_load(open(config_path,'r'))

    config_vars = set_config_defaults(config_obj)
    config_vars["script_dir"] = script_dir
    config_vars["usr_cfg_path"] = config_path
    config_vars["cfg_path"] = config_obj.get('output_path')+'tmp/tmp.cfg.yml'


    # write config to file and to console
    write_cfg2file(config_vars)
    write_run_overview(config_path, config_vars)


def make_checks(config_obj):
    """Make data validity checks.

    Cross checks scaffold IDs in FASTA and GFF

    Args:
      config_obj(obj): Config class object of config parameters
    """
    check_assembly_ids(config_obj)


def cfg2obj(config_path):
    """Make Config class object of config parameters from config file.

    Args:
      config_path(str): path to temporary processed config file

    Returns:
      (obj) Config class object with parameters
    """
    config_obj = yaml.safe_load(open(config_path,'r'))
    cfg = Config(config_obj)

    return cfg


def main():
    """Call module directly with preprocessed config file"""

    config_path = sys.argv[1]
    script_dir = sys.argv[2]

    # read config file and set defaults
    # print to file and summary to console
    process_config(config_path, script_dir)

    # check input data
    make_checks(yaml.safe_load(open(config_path,'r')))


if __name__ == '__main__':
    main()
