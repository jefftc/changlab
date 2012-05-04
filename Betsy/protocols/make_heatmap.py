#make_heatmap.py


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
    
OUTPUTS = ['cluster_heatmap']

PARAMETERS = {
    'cluster_alg':['no_cluster_alg'],
    'preprocess':['rma','mas5','loess','illumina_controls',
                  'ilumina','agilent','unknown_preprocess'],
    'gene_center':['mean','median','no_gene_center'],
    'gene_normalize':['variance','sum_of_squares','no_gene_normalize'],
    'quantile':['yes_quantile','no_quantile'],
    'dwd':['yes_dwd','no_dwd'],
    'shiftscale':['yes_shiftscale','no_shiftscale'],
    'combat':['yes_combat','no_combat'],
    'gene_order':['by_sample_ttest','by_gene_list','no_order','by_class_neighbors'],
    'color':['red_green','blue_yellow'],
    'predataset':['yes_predataset','no_predataset'],
    'filter':['no_filter','integer between 0 and 100'],
    'has_missing_value':['no_missing','median_fill','zero_fill','unknown_missing'],
    'is_logged':['no_logged','logged'],
    'traincontents':['arbitrary content'],
    'testcontents':['arbitrary content'],
    'contents':['arbitrary content'],
    'hm_width':['number'],
    'hm_height':['number'],
    'ill_manifest':['arbitrary string'],
    'ill_chip':['arbitrary string'],
    'ill_bg_mode':['ill_yes','ill_no'],
    'ill_coll_mode':['ill_none','ill_max','ill_median'],
    'ill_clm':['arbitrary string'],
    'ill_custom_chip':['arbitrary string'],
    'cn_num_neighbors':['number'],
    'cn_num_perm':['number'],
    'cn_user_pval':['number'],
    'cn_mean_or_median':['cn_mean','cn_median'],
    'cn_ttest_or_snr':['cn_test','snr'],
    'cn_filter_data':['cn_yes','cn_no'],
    'cn_min_threshold':['number'],
    'cn_max_threshold':['number'],
    'cn_min_folddiff':['number'],
    'cn_abs_diff':['number']}


DEFAULT = {
    'cluster_alg':'no_cluster_alg',
    'gene_center':'no_gene_center',
    'gene_normalize':'no_gene_normalize','quantile':'no_quantile',
    'dwd':'no_dwd','shiftscale':'no_shiftscale',
    'combat':'no_combat','filter':'no_filter',
    'predataset':'no_predataset','is_logged':'logged',
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

