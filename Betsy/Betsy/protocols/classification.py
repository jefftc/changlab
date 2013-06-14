#classification.py
from Betsy.protocol_utils import Parameter

PRETTY_NAME="Classify samples."
COMMON = 'Common Parameters'
NORMALIZE = 'Normalize Parameters'
OPTIONAL = 'Optional Parameters'
ILLUMINA = 'Illumina Normalize Parameters'
CLASS_NEIGHBORS='Class Neighbor Parameters'
CLASSIFY='Classification parameters'
CATEGORIES=[COMMON,NORMALIZE,OPTIONAL,ILLUMINA,CLASS_NEIGHBORS,CLASSIFY]
#input predicates
INPUTS = [
    'gse_id',
    'gse_id_and_platform',
    'cel_files',
    'gpr_files',
    'idat_files',
    'agilent_files',
    'input_signal_file',
    'class_label_file',
    'gene_list_file',
    'renanme_list_file']

#output predicates
OUTPUTS = 'make_classify_report'

#convert predicates to prolog facts
predicate2arguments = {
    'gse_id': ([], '[]'),
    'gse_id_and_platform': ([], '[]'),
    'idat_files': (['version', 'unknown_version'], '[]'),
    'gpr_files': (['version', 'unknown_version'], '[]'),
    'cel_files': (['version', 'unknown_version'], '[]'),
    'agilent_files': (['version', 'unknown_version'], '[]'),
    'class_label_file': (['status', 'given'], '[]'),
    'gene_list_file': ([], '[]'),
    'input_signal_file': (['status', 'given'], '[]'),
    'rename_list_file': ([], '[]')}

#parameter objects
PARAMETERS=[Parameter('preprocess',pretty_name='Preprocess',
                         choices=['rma', 'mas5',
                        'loess', 'illumina_controls',
                        'illumina', 'agilent',
                        'unknown_preprocess'],category=COMMON),
            Parameter('gene_center',pretty_name='Gene Center',
                        default_value='no_gene_center', 
                        choices=['mean', 'median', 'no_gene_center'],
                        category=COMMON),
            Parameter('gene_normalize',pretty_name='Gene Normalize',
                           default_value='no_gene_normalize',
                           choices=['variance', 'sum_of_squares', 'no_gene_normalize'],
                           category=COMMON),
            Parameter('quantile',pretty_name='Quantile',default_value='no_quantile', 
                     choices=['yes_quantile', 'no_quantile'],category=NORMALIZE),
            Parameter('dwd',pretty_name='Dwd',default_value='no_dwd',choices=['yes_dwd', 'no_dwd'],
                category=NORMALIZE),
            Parameter('bfrm',pretty_name='BFRM',default_value='no_bfrm',
                 choices=['yes_bfrm', 'no_bfrm'],category=NORMALIZE),
            Parameter('shiftscale',pretty_name='Shift Scale',default_value='no_shiftscale',
                       choices=['yes_shiftscale', 'no_shiftscale'],category=NORMALIZE),
            Parameter('combat',pretty_name='Combat',default_value='no_combat',
                   choices=['yes_combat', 'no_combat'],category=NORMALIZE),
            Parameter('gene_order',pretty_name='Gene Order',default_value='no_order',
                       choices=['t_test_p', 't_test_fdr', 'by_gene_list', 'no_order',
                   'by_class_neighbors'],category=COMMON),
            Parameter('predataset',pretty_name='Predataset Process',default_value='no_predataset',
                       choices=['yes_predataset', 'no_predataset'],category=COMMON),
            Parameter('filter',pretty_name='Filter',default_value='0',type='integer',
                      category=COMMON),
            Parameter('has_missing_value',pretty_name='How to fill missing value',
                              choices=['no_missing', 'median_fill', 'zero_fill',
                          'unknown_missing'],category=COMMON),
            Parameter('format',pretty_name='File Format',default_value='tdf',
                        choices=['tdf', 'gct'],category=COMMON),
            Parameter('is_logged',pretty_name='Logged or Not',default_value='logged',
                      choices=['no_logged', 'logged'],category=COMMON),
            Parameter('contents',pretty_name='Contents',type='list',category=OPTIONAL),
            Parameter('num_features',pretty_name='Feature Number',default_value='0',
                         type='integer',category=COMMON),
            Parameter('ill_manifest', pretty_name='Illumina Manifest File',
                          type='string',category=ILLUMINA),
            Parameter('ill_chip', pretty_name='Illumina chip File',type='string',
                     category=ILLUMINA),
            Parameter('ill_clm', pretty_name='Illumina clm File',type='string',
                     category=ILLUMINA),
            Parameter('ill_custom_chip', pretty_name='Illumina Custom Chip File',
                             type='string',category=ILLUMINA),
            Parameter('ill_custom_manifest', pretty_name='Illumina Custom Manifest File',
                                 type='string',category=ILLUMINA),
            Parameter('ill_bg_mode', pretty_name='Illumina Background Mode',
                         choices=['ill_yes', 'ill_no'],category=ILLUMINA),
            Parameter('ill_coll_mode', pretty_name='Illumina Coll Mode',
                           choices=['ill_none', 'ill_max', 'ill_median'],category=ILLUMINA),
            Parameter('cn_num_neighbors',pretty_name='Class Neighbors Number',
                             type='integer',category=CLASS_NEIGHBORS),
            Parameter('cn_num_perm',pretty_name='Class Permutation Number',
                        type='integer',category=CLASS_NEIGHBORS),
            Parameter('cn_user_pval',pretty_name='Class User Pval',type='float',
                          category=CLASS_NEIGHBORS),
            Parameter('cn_mean_or_median',pretty_name='Class Neighbors Mean or Median',
                              choices=['cn_mean', 'cn_median'],category=CLASS_NEIGHBORS),
            Parameter('cn_ttest_or_snr',pretty_name='Class Neighbors ttest or snr',
                            choices=['cn_test', 'cn_snr'],category=CLASS_NEIGHBORS),
            Parameter('cn_filter_data',pretty_name='Class Neighbors',
                           choices=['cn_yes', 'cn_no'],category=CLASS_NEIGHBORS),
            Parameter('cn_min_threshold',pretty_name='Class neighbors Min Threshold',
                             type='float',category=CLASS_NEIGHBORS),
            Parameter('cn_max_threshold',pretty_name='Class neighbors Max Threshold',
                             type='float',category=CLASS_NEIGHBORS),
            Parameter('cn_min_folddiff',pretty_name='Class neighbors Min Fold Diff',
                            type='float',category=CLASS_NEIGHBORS),
            Parameter('cn_abs_diff','Class neighbors Abs Diff',None, 'float',
                             None,CLASS_NEIGHBORS,''),
            Parameter('gene_select_threshold',pretty_name='Gene Selection Threshold',
                                  type= 'float',category=OPTIONAL),
            Parameter('rename_sample',pretty_name='Rename Sample',default_value='no_rename',
                          choices=['yes_rename', 'no_rename'], category=OPTIONAL),
            Parameter('num_factors',pretty_name='Number of Factors',type='integer',
                           category=OPTIONAL),
            Parameter('pca_gene_num',pretty_name='PCA Gene Number',type='integer',
                           category=OPTIONAL),
            Parameter('unique_genes',pretty_name='Unique Genes',default_value='no_unique_genes',
                         choices=['average_genes', 'high_var', 'first_gene','no_unique_genes'],
                         category=COMMON),
            Parameter('platform',pretty_name='Platform',choices=["'HG_U133_Plus_2'", "'HG_U133B'", "'HG_U133A'",
                 "'HG_U133A_2'", "'HG_U95A'", "'HumanHT_12'", "'HumanWG_6'","'HG_U95Av2'",
                 "'Entrez_ID_human'", "'Entrez_symbol_human'", "'Hu6800'",
                 "'Mouse430A_2'", "'MG_U74Cv2'", "'Mu11KsubB'", "'Mu11KsubA'",
                 "'MG_U74Av2'", "'Mouse430_2'", "'MG_U74Bv2'",
                 "'Entrez_ID_mouse'", "'MouseRef_8'", "'Entrez_symbol_mouse'",
                 "'RG_U34A'", "'RAE230A'", 'unknown_platform'],category=OPTIONAL),
            Parameter('duplicate_probe',pretty_name='Duplicate Probe',
                            choices=['yes_duplicate_probe', 'high_var_probe','closest_probe'],
                            category=OPTIONAL),
            Parameter('has_annotation_gene_id',pretty_name='Annotation Gene Id',
                                   choices=['yes_gene_id','no_gene_id'],category=OPTIONAL),
            Parameter('group_fc', pretty_name='Fold Change',type='integer',category=COMMON),
            Parameter('classification',pretty_name='Classification Method',
                                   choices=['svm', 'weightedvoting', 'randomforest'],category=CLASSIFY),
            Parameter('svm_kernel',pretty_name='SVM Kernel',
                                   choices=['linear', 'polynomial', 'rbf', 'sigmoid',
                   'precomputed_kernel'],category=CLASSIFY),
            Parameter('wv_feature_stat',pretty_name='Weight Voted Feature Stat',
                                   choices=['wv_snr', 'wv_ttest', 'wv_snr_median',
                        'wv_ttest_median',
                        'wv_snr_minstd', 'wv_ttest_minstd',
                        'wv_snr_median_minstd',
                        'wv_ttest_median_minstd'],category=CLASSIFY),
            Parameter('traincontents',pretty_name='Training contents',type='list',category=CLASSIFY),
            Parameter('testcontents',pretty_name='Test contents',type='list',category=CLASSIFY),
            Parameter('wv_minstd',pretty_name='Weight Voted min std',type='float',category=CLASSIFY)]

















