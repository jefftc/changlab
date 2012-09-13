#differential_expressed_gene_analysis.py

INPUTS = [
    'gse_id',
    'class_label_file',
    'gene_list_file',
    'gse_id_and_platform',
    'cel_files',
    'gpr_files',
    'idat_files',
    'agilent_files',
    'input_signal_file',
    'renanme_list_file']

OUTPUTS = 'make_diffgenes_report'

PARAMETERS = {
    'preprocess':['rma','mas5','loess','illumina_controls',
                  'illumina','agilent','unknown_preprocess'],
    'has_missing_value':['unknown_missing','no_missing','median_fill','zero_fill'],
    'gene_center':['mean','median','no_gene_center'],
    'gene_normalize':['variance','sum_of_squares','no_gene_normalize'],
    'quantile':['yes_quantile','no_quantile'],
    'dwd':['yes_dwd','no_dwd'],
    'bfrm':['yes_bfrm','no_bfrm'],
    'shiftscale':['yes_shiftscale','no_shiftscale'],
    'combat':['yes_combat','no_combat'],
    'gene_order':['t_test_p','t_test_fdr','by_gene_list','no_order','by_class_neighbors'],
    'predataset':['yes_predataset','no_predataset'],
    'filter':'integer',
    'format':['tdf','gct'],
    'is_logged':['no_logged','logged'],
    'traincontents':'list',
    'testcontents':'list',
    'contents':'list',
    'diff_expr':['t_test','sam'],
    'sam_foldchange':'float',
    'sam_delta':'float',
    'ill_manifest':'string',
    'ill_chip':'integer',
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
    'cn_abs_diff':'float',
    'gene_select_threshold':'float',
    'rename_sample':['yes_rename','no_rename'],
    'num_factors':'integer',
    'pca_gene_num':'integer',
    'unique_genes':['average_genes','high_var','first_gene','no_unique_genes'],
    'platform':["'HG_U133_Plus_2'","'HG_U133B'","'HG_U133A'","'HG_U133A_2'",
                "'HG_U95A'","'HumanHT_12'","'HG_U95Av2'","'entrez_ID_human'",
                 "'entrez_ID_symbol_human'","'Hu6800'","'Mouse430A_2'", "'MG_U74Cv2'",
                 "'Mu11KsubB','Mu11KsubA'","'MG_U74Av2'","'Mouse430_2'","'MG_U74Bv2'",
                "'entrez_ID_mouse'","'MouseRef_8'",
                 "'entrez_ID_symbol_mouse'","'RG_U34A'","'RAE230A'",'unknown_platform'],
    'duplicate_probe':['yes_duplicate_probe','high_var_probe','closest_probe'],
    'duplicate_data':['yes_duplicate_data','no_duplicate_data'],
    'missing_probe':['yes_missing_probe','no_missing_probe']
    }


DEFAULT = {
    'gene_center':'no_gene_center',
    'gene_normalize':'no_gene_normalize','quantile':'no_quantile',
    'dwd':'no_dwd','shiftscale':'no_shiftscale',
    'combat':'no_combat','filter':0,'is_logged':'logged',
    'predataset':'no_predataset','bfrm':'no_bfrm',
    'gene_order':'no_order','rename_sample':'no_rename','format':'tdf',
    'duplicate_probe':'yes_duplicate_probe',
    'duplicate_data':'yes_duplicate_data',
    'unique_genes':'no_unique_genes'}

predicate2arguments={
    'gse_id':([],'[]'),
    'gse_id_and_platform':([],'[]'),
    'idat_files':(['version','unknown_version'],'[]'),
    'gpr_files':(['version','unknown_version'],'[]'),
    'cel_files':(['version','unknown_version'],'[]'),
    'agilent_files':(['version','unknown_version'],'[]'),
    'class_label_file':(['status','given'],'[]'),
    'gene_list_file':([],'[]'),
    'input_signal_file':(['status','given'],'[]'),
    'rename_list_file':([],'[]')}

     
