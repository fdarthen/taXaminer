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

    if key_name in config_obj.keys():
        return config_obj[key_name]
    elif default: # if default value is given this is evaluated to True
        return default[0]
    else: # key not given in config file and no given default value
        print('Error: no default value available for "{}". Please state specifically.'.format(key_name))



def set_config_defaults(config_obj):
    """ check config file to set default values and write to temporary config file """

    config_vars = {}

    default_pca_vars = "c_name,c_num_of_genes,c_len,c_genelenm,c_genelensd,g_len,g_lendev_c,g_abspos,g_terminal,c_cov,c_covsd,g_cov,g_covsd,g_covdev_c,c_pearson_r,g_pearson_r_o,g_pearson_r_c"

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

    if len(cov_set_exists) == 0:
        config_vars['include_coverage'] = set_variable_default(config_obj, 'include_coverage', 'FALSE')
    else:
        config_vars['include_coverage'] = set_variable_default(config_obj, 'include_coverage', 'TRUE')

    if list(cov_set_exists.values()).count('pbc_paths') == len(cov_set_exists) or \
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
    config_vars['input_variables'] = set_variable_default(config_obj, 'input_variables', default_pca_vars)
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


def main():

    config_path = sys.argv[1]

    # read parameters from config file
    config_obj = yaml.safe_load(open(config_path,'r'))

    config_vars = set_config_defaults(config_obj)

    write_cfg2file(config_obj, config_vars)


if __name__ == '__main__':
    main()
