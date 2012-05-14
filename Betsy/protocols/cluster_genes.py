#cluster_genes.py


INPUTS = [
    'gse_id',
    'class_label_file',
    'gene_list_file',
    'gse_id_and_platform',
    'cel_files',
    'gpr_files',
    'idat_files',
    'agilent_files',
    'input_signal_file']

OUTPUTS = ['cluster_heatmap','cluster_file']

PARAMETERS = {
    'cluster_alg':['kmeans','pca','hierarchical','som'],
    'distance':['correlation','euclidean'],
    'k':'integer',
    'preprocess':['rma','mas5','loess','illumina_controls',
                  'ilumina','agilent','unknown_preprocess'],
    'gene_center':['mean','median','no_gene_center'],
    'gene_normalize':['variance','sum_of_squares','no_gene_normalize'],
    'quantile':['yes_quantile','no_quantile'],
    'dwd':['yes_dwd','no_dwd'],
    'shiftscale':['yes_shiftscale','no_shiftscale'],
    'combat':['yes_combat','no_combat'],
    'gene_order':['by_sample_ttest','by_gene_list','no_order',
                  'by_class_neighbors'],
    'color':['red_green','blue_yellow'],
    'predataset':['yes_predataset','no_predataset'],
    'filter':'integer',
    'has_missing_value':['no_missing','median_fill','zero_fill',
                         'unknown_missing'],
    'is_logged':['no_logged','logged'],
    'traincontents':'list',
    'testcontents':'list',
    'contents':'list',
    'hm_width':'integer',
    'hm_height':'integer',
    'ill_manifest':'string',
    'ill_chip':'string',
    'ill_bg_mode':['ill_yes','ill_no'],
    'ill_coll_mode':['ill_none','ill_max','ill_median'],
    'ill_clm':'string',
    'ill_custom_chip':'string',
    'ill_custom_manifest':'string',
    'cn_num_neighbors':'integer',
    'cn_num_perm':'integer',
    'cn_user_pval':'float',
    'cn_mean_or_median':['cn_mean','cn_median'],
    'cn_ttest_or_snr':['cn_test','cn_snr'],
    'cn_filter_data':['cn_yes','cn_no'],
    'cn_min_threshold':'float',
    'cn_max_threshold':'float',
    'cn_min_folddiff':'float',
    'cn_abs_diff':'float'}

DEFAULT = {
    'cluster_alg':'kmeans','distance':'correlation','k':5,
    'gene_center':'no_gene_center',
    'gene_normalize':'no_gene_normalize','quantile':'no_quantile',
    'dwd':'no_dwd','shiftscale':'no_shiftscale',
    'combat':'no_combat','filter':0,
    'predataset':'no_predataset',
    'is_logged':'logged',
    'gene_order':'no_order','color':'red_green'}

predicate2arguments={
    'gse_id':([],'[]'),
    'gse_id_and_platform':([],'[]'),
    'idat_files':(['version','unknown_version'],'[]'),
    'gpr_files':(['version','unknown_version'],'[]'),
    'cel_files':(['version','unknown_version'],'[]'),
    'agilent_files':(['version','unknown_version'],'[]'),
    'class_label_file':(['status','given'],'[]'),
    'gene_list_file':([],'[]'),
    'input_signal_file':(['status','given'],'[]')}







