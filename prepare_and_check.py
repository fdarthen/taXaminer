# -*- coding: utf-8 -*-

import pathlib # to create directories
import yaml # read config file
import sys # parse command line arguments
import os

def enumerated_key(config_obj, key_name, pre_keys, *default):
    matches = []
    for key in config_obj.keys():
        if key_name in key:
            matches.append(key)

    dict = {}
    for match in matches:
        if match.split('_')[2].isdigit():
            match_num = int(match.split('_')[2])
            dict[match_num] = set_variable_default(config_obj, match, default)

    for pre_key in pre_keys:
        # if there are more values required than given, fill them up with default paths
        if not dict.get(pre_key):
            if '/' in default:
                dict[pre_key] = '.'.join(default[0].split('.')[:-1])+'_'+str(pre_key)+'.'+default[0].split('.')[-1]
            else:
                dict[pre_key] = default[0]

    return dict


def check_dir_slash(path):
    """ check if path ends on slash """
    if path.endswith('/'):
        return path
    else:
        return path+'/'

def check_file_inexistence(paths, and_or):
    """ check for a run option bool (extract_proteins, compute_coverage, etc.) if the required file exists """

    if type(paths) != list:
        paths = [paths]

    counter_e = 0
    for path in paths:
        file = pathlib.Path(path)
        if file.is_file():
            counter_e += 1


    if (and_or == 'AND' and counter_e == len(paths)) or \
     (and_or == 'OR' and counter_e > 0):
        return 'FALSE'
    else:
        return 'TRUE'

def set_variable_default(config_obj, key_name, *default):
    """ check if variable is given in users config file
    else set to default value """

    if config_obj.get(key_name):
        return config_obj[key_name]
    elif default: # if default value is given this is evaluated to True
        return default[0]
    else: # key not given in config file and no given default value
        print('Error: no default value available for "{}". Please state specifically.'.format(key_name))

def pca_cov_variables(vars, include, cov_set_num):
    """ exclude cov_vars if coverage_include == FALSE;
    else add indices to the cov_vars """

    cov_vars = "c_cov,c_covsd,g_cov,g_covsd,g_covdev_c"

    out_vars = ""

    if include:
        for var in vars.split(','):
            if var in cov_vars:
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
    """ check config file to set default values and write to temporary config file """

    config_vars = {}

    default_pca_vars = "c_name,c_num_of_genes,c_len,c_genelenm,c_genelensd,g_len,g_lendev_c,g_abspos,g_terminal,c_pearson_r,g_pearson_r_o,g_pearson_r_c,c_cov,c_covsd,g_cov,g_covsd,g_covdev_c"

    # General options
    config_vars['fasta_path'] = set_variable_default(config_obj, 'fasta_path')
    config_vars['gff_path'] = set_variable_default(config_obj, 'gff_path')
    config_vars['output_path'] = check_dir_slash(set_variable_default(config_obj, 'output_path'))
    config_vars['taxon_id'] = set_variable_default(config_obj, 'taxon_id')

    # Coverage
    config_vars['read_paths'] = enumerated_key(config_obj, 'read_paths', [])
    config_vars['bam_paths'] = enumerated_key(config_obj, 'bam_path', list(config_vars.get('read_paths').keys()), config_vars.get('output_path')+'mapping_sorted.bam')
    config_vars['pbc_paths'] = enumerated_key(config_obj, 'pbc_path', list(config_vars.get('bam_paths').keys()), config_vars.get('output_path')+'pbc.txt')

    cov_set_exists = {}
    for key in config_vars.get('pbc_paths').keys():
        for cov_info in ['read_paths', 'bam_paths', 'pbc_paths']:
            path = config_vars.get(cov_info).get(key)
            if path:
                if check_file_inexistence(path, 'AND') == 'FALSE':
                    # file exists
                    cov_set_exists[key] = cov_info

    config_vars['cov_set_exists'] = cov_set_exists

    if len(config_vars.get('cov_set_exists')) == 0:
        config_vars['include_coverage'] = set_variable_default(config_obj, 'include_coverage', 'FALSE')
    else:
        config_vars['include_coverage'] = set_variable_default(config_obj, 'include_coverage', 'TRUE')

    if list(config_vars.get('cov_set_exists').values()).count('pbc_paths') == len(config_vars.get('cov_set_exists')) or \
        config_vars.get('include_coverage') == 'FALSE':
        config_vars['compute_coverage'] = set_variable_default(config_obj, 'compute_coverage', 'FALSE')
    else:
        config_vars['compute_coverage'] = set_variable_default(config_obj, 'compute_coverage', 'TRUE')

    config_vars['insert_size'] = enumerated_key(config_obj, 'insert_size', list(config_vars.get('read_paths').keys()), '200')
    # Taxonomic assignment
    config_vars['proteins_path'] = set_variable_default(config_obj, 'proteins_path', config_vars.get('output_path')+'proteins.faa')
    config_vars['extract_proteins'] = set_variable_default(config_obj, 'extract_proteins', check_file_inexistence(config_vars.get('proteins_path'), 'AND'))
    config_vars['assignment_mode'] = set_variable_default(config_obj, 'assignment_mode', 'exhaustive')
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
    # Plotting
    config_vars['update_plots'] = set_variable_default(config_obj, 'update_plots', 'FALSE')
    config_vars['num_groups_plot'] = set_variable_default(config_obj, 'num_groups_plot', '25')
    config_vars['merging_labels'] = set_variable_default(config_obj, 'merging_labels', 'None')
    config_vars['output_pdf'] = set_variable_default(config_obj, 'output_pdf', 'TRUE')
    config_vars['output_png'] = set_variable_default(config_obj, 'output_png', 'FALSE')
    # Gene info
    config_vars['include_pseudogenes'] = set_variable_default(config_obj, 'include_pseudogenes', 'FALSE')
    config_vars['gff_source'] = set_variable_default(config_obj, 'gff_source', 'default')

    # PCA
    config_vars['input_variables'] = pca_cov_variables(set_variable_default(config_obj, 'input_variables', default_pca_vars), config_vars.get('include_coverage'), len(config_vars.get('cov_set_exists')))

    config_vars['perform_parallel_analysis'] = set_variable_default(config_obj, 'perform_parallel_analysis', 'FALSE')
    config_vars['num_pcs'] = set_variable_default(config_obj, 'num_pcs', '3')
    config_vars['coverage_cutoff_mode'] = set_variable_default(config_obj, 'coverage_cutoff_mode', 'default')
    # Clustering
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


def write_cfg2file(config_obj, config_vars):

    with open(config_obj.get('output_path')+'tmp/tmp.cfg.yml', 'w') as out_cfg:


        for key, value in config_vars.items():
                out_cfg.write('{}: {}\n'.format(key,value))

def write_run_overview(config_path, config_vars):

    print('')
    print("Config:\t{}".format(config_path))
    print("FASTA:\t{}".format(config_vars.get('fasta_path')))
    print("GFF:\t{}".format(config_vars.get('gff_path')))
    print("Taxon ID:\t{}".format(config_vars.get('taxon_id')))
    print("Output:\t{}\n".format(config_vars.get('output_path')))

    if config_vars.get('include_coverage') == 'TRUE':
        for cov_set, pbc_path in config_vars.get('pbc_paths').items():
            if config_vars.get('cov_set_exists').get(cov_set) == 'pbc_paths':
                print("PBC {}:\t{} [exists]".format(cov_set, pbc_path))
            else:
                if config_vars.get('cov_set_exists').get(cov_set) == 'bam_paths':
                    basis = 'BAM'
                else:
                    basis = 'read mapping'
                print("PBC {}:\t{}".format(cov_set, pbc_path))
                print("  computation based on:\t{} [{}]".format( config_vars.get(config_vars.get('cov_set_exists').get(cov_set)).get(cov_set), basis))
        print('')
    else:
        print("no valid coverage information provided\n")

    if config_vars.get('extract_proteins') == 'FALSE':
        print("Proteins:\t{} [exists]".format(config_vars.get('proteins_path')))
    else:
        print("Proteins:\t{}".format(config_vars.get('proteins_path')))
    if config_vars.get('assignment_mode') == 'quick':
        print("Quick assignment mode selected")
        print("Filtering search performed on level {}".format(config_vars.get('quick_mode_search_rank')))
        print("Hits accepted on level {}".format(config_vars.get('quick_mode_match_rank')))
        if config_vars.get('compute_tax_assignment') == 'FALSE':
            print("Taxonomic hits files [exist]:\n{}{}".format(config_vars.get('tax_assignment_path')[0],config_vars.get('tax_assignment_path')[1]))
        else:
            print("Taxonomic hits files:\n{}{}".format(config_vars.get('tax_assignment_path')[0],config_vars.get('tax_assignment_path')[1]))


    else:
        print("Exhaustive assignment mode selected")
        if config_vars.get('compute_tax_assignment') == 'FALSE':
            print("Taxonomic hits file:{} [exists]".format(config_vars.get('tax_assignment_path')))
        else:
            print("Taxonomic hits file:{}".format(config_vars.get('tax_assignment_path')))

    print("Query taxon excluded:\t{}\n".format(config_vars.get('taxon_exclude')))

    print("Pseudogenes included:\t{}".format(config_vars.get('include_pseudogenes')))
    if config_vars.get('gff_source') != 'default':
        print("Rule for GFF parsing:\t{}".format(config_vars.get('gff_source')))
    print('\n')

    print("PCA variables:\t{}".format(config_vars.get('input_variables')))
    print("Parallel analysis performed:\t{}".format(config_vars.get('perform_parallel_analysis')))
    if config_vars.get('coverage_cutoff_mode') != 'default':
        print("Coverage cutoff:\t{}".format(config_vars.get('coverage_cutoff_mode')))

    for clustering in ['kmeans', 'hclust', 'mclust']:
        if config_vars.get('perform_'+clustering) == 'TRUE':
            print("{} clustering performed with number of groups:\t{}".format(clustering, config_vars.get(clustering+'k')))
    if config_vars.get('perform_dbscan') == 'TRUE':
        print("DBSCAN clustering performed with settings:\t{}".format(config_vars.get('dbscan_groups')))
    print('')


def main():

    config_path = sys.argv[1]

    # read parameters from config file
    config_obj = yaml.safe_load(open(config_path,'r'))

    config_vars = set_config_defaults(config_obj)

    write_cfg2file(config_obj, config_vars)

    write_run_overview(config_path, config_vars)



if __name__ == '__main__':
    main()
