#batch_effect_remove.py
from Betsy.protocol_utils import Parameter
from Betsy import protocol_utils

PRETTY_NAME="Remove batch effects."
COMMON = protocol_utils.COMMON
NORMALIZE = protocol_utils.NORMALIZE
OPTIONAL = protocol_utils.OPTIONAL
ILLUMINA = protocol_utils.ILLUMINA
CLASS_NEIGHBORS = protocol_utils.CLASS_NEIGHBORS
column_name=[COMMON,NORMALIZE,OPTIONAL,ILLUMINA,CLASS_NEIGHBORS]
#input predicates
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
    'rename_list_file']

#output predicates
OUTPUTS = 'make_batch_report'

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
preprocess = Parameter('preprocess','Preprocess',
                         None, None, ['rma', 'mas5',
                        'loess', 'illumina_controls',
                        'illumina', 'agilent',
                        'unknown_preprocess'],COMMON,
                        '')
gene_center = Parameter('gene_center','Gene Center',
                        'no_gene_center', None,
                        ['mean', 'median', 'no_gene_center'],
                        COMMON,'')
gene_normalize = Parameter('gene_normalize','Gene Normalize',
                           'no_gene_normalize',None,
                           ['variance', 'sum_of_squares', 'no_gene_normalize'],
                           NORMALIZE,'')
quantile = Parameter('quantile','Quantile','no_quantile', None,
                     ['yes_quantile', 'no_quantile'],COMMON,'')
dwd = Parameter('dwd','Dwd','no_dwd',None,['yes_dwd', 'no_dwd'],
                NORMALIZE,'')
bfrm = Parameter('bfrm','BFRM','no_bfrm',None,
                 ['yes_bfrm', 'no_bfrm'],NORMALIZE,'')
shiftscale = Parameter('shiftscale','Shift Scale','no_shiftscale',None,
                       ['yes_shiftscale', 'no_shiftscale'],NORMALIZE,'')
combat = Parameter('combat','Combat','no_combat',None,
                   ['yes_combat', 'no_combat'],NORMALIZE,'')
gene_order = Parameter('gene_order','Gene Order','no_order',None,
                       ['t_test_p', 't_test_fdr', 'by_gene_list', 'no_order',
                   'by_class_neighbors'],COMMON,'')
predataset = Parameter('predataset','Predataset Process','no_predataset',None,
                       ['yes_predataset', 'no_predataset'],COMMON,'')
filtering = Parameter('filter','Filter','0','integer',None,COMMON,'')
has_missing_value = Parameter('has_missing_value','How to fill missing value',
                              None,None,['no_missing', 'median_fill', 'zero_fill',
                          'unknown_missing'],COMMON,'')
format_type = Parameter('format','File Format','tdf',None,['tdf', 'gct'],COMMON,'')
is_logged = Parameter('is_logged','Logged or Not','logged',None,['no_logged', 'logged'],
                      COMMON,'')
contents = Parameter('contents','Contents',None,'list',None,OPTIONAL,'')
num_features = Parameter('num_features','Feature Number','0','integer',None,COMMON,'')
illu_manifest = Parameter('ill_manifest', 'Illumina Manifest File',None,
                          'string',None,ILLUMINA,'')
illu_chip= Parameter('ill_chip', 'Illumina chip File',None,'string',None,ILLUMINA,'')
illu_clm = Parameter('ill_clm', 'Illumina clm File',None,'string',None,ILLUMINA,'')
illu_custom_chip = Parameter('ill_custom_chip', 'Illumina Custom Chip File',None,
                             'string',None,ILLUMINA,'')
illu_custom_manifest = Parameter('ill_custom_manifest', 'Illumina Custom Manifest File',
                                 None,'string',None,ILLUMINA,'')
illu_bg_mode = Parameter('ill_bg_mode', 'Illumina Background Mode',None,None,
                         ['ill_yes', 'ill_no'],ILLUMINA,'')
illu_coll_mode = Parameter('ill_coll_mode', 'Illumina Coll Mode',None,None,
                           ['ill_none', 'ill_max', 'ill_median'],ILLUMINA,'')

cn_num_neighbors = Parameter('cn_num_neighbors','Class Neighbors Number',None,'integer',None,
                          CLASS_NEIGHBORS,'')
cn_num_perm = Parameter('cn_num_perm','Class Permutation Number',None,'integer',None,
                          CLASS_NEIGHBORS,'')
cn_user_pval = Parameter('cn_user_pval','Class User Pval',None,'float',None,
                          CLASS_NEIGHBORS,'')
cn_mean_or_median = Parameter('cn_mean_or_median','Class Neighbors Mean or Median',None,None,
                             ['cn_mean', 'cn_median'],CLASS_NEIGHBORS,'')
cn_ttest_or_snr = Parameter('cn_ttest_or_snr','Class Neighbors ttest or snr',None,None,
                            ['cn_test', 'cn_snr'],CLASS_NEIGHBORS,'')
cn_filter_data = Parameter('cn_filter_data','Class Neighbors',None,None,
                           ['cn_yes', 'cn_no'],CLASS_NEIGHBORS,'')
cn_min_threshold = Parameter('cn_min_threshold','Class neighbors Min Threshold',None, 'float',
                             None,CLASS_NEIGHBORS,'')
cn_max_threshold = Parameter('cn_max_threshold','Class neighbors Max Threshold',None, 'float',
                             None,CLASS_NEIGHBORS,'')
cn_min_folddiff = Parameter('cn_min_folddiff','Class neighbors Min Fold Diff',None, 'float',
                             None,CLASS_NEIGHBORS,'')
cn_abs_folddiff = Parameter('cn_abs_diff','Class neighbors Abs Diff',None, 'float',
                             None,CLASS_NEIGHBORS,'')
gene_select_threshold = Parameter('gene_select_threshold','Gene Selection Threshold',None, 'float',
                             None,OPTIONAL,'')
rename_sample = Parameter('rename_sample','Rename Sample','no_rename',None,
                          ['yes_rename', 'no_rename'], OPTIONAL,'')
num_factors = Parameter('num_factors','Number of Factors',None,'integer',None,
                           OPTIONAL,'')
pca_gene_num = Parameter('pca_gene_num','PCA Gene Number',None,'integer',None,
                           OPTIONAL,'')
unique_genes = Parameter('unique_genes','Unique Genes','no_unique_genes',None,
                         ['average_genes', 'high_var', 'first_gene','no_unique_genes'],
                         COMMON,'')
platform = Parameter('platform','Platform',None,None,["'HG_U133_Plus_2'", "'HG_U133B'", "'HG_U133A'",
                 "'HG_U133A_2'", "'HG_U95A'", "'HumanHT_12'", "'HumanWG_6'","'HG_U95Av2'",
                 "'Entrez_ID_human'", "'Entrez_symbol_human'", "'Hu6800'",
                 "'Mouse430A_2'", "'MG_U74Cv2'", "'Mu11KsubB'", "'Mu11KsubA'",
                 "'MG_U74Av2'", "'Mouse430_2'", "'MG_U74Bv2'",
                 "'Entrez_ID_mouse'", "'MouseRef_8'", "'Entrez_symbol_mouse'",
                 "'RG_U34A'", "'RAE230A'", 'unknown_platform'],OPTIONAL,'')
duplicate_probe = Parameter('duplicate_probe','Duplicate Probe',None,None,
                            ['yes_duplicate_probe', 'high_var_probe','closest_probe'],OPTIONAL,'')
has_annotation_gene_id = Parameter('has_annotation_gene_id','Annotation Gene Id',None,None,
                                   ['yes_gene_id','no_gene_id'],OPTIONAL,'')
group_fc = Parameter('group_fc', 'Fold Change',None,'integer',None,COMMON,'')

PARAMETERS_list = [preprocess,gene_center,gene_normalize,quantile,dwd,bfrm,shiftscale,combat,gene_order,
              predataset,filtering,has_missing_value,format_type,is_logged,contents,num_features,
              illu_manifest,illu_chip,illu_clm,illu_custom_chip,illu_custom_manifest,illu_bg_mode,
              illu_coll_mode,cn_num_neighbors,cn_num_perm,cn_user_pval,cn_mean_or_median,
              cn_ttest_or_snr,cn_filter_data,cn_min_threshold,cn_max_threshold,cn_min_folddiff,
              cn_abs_folddiff,gene_select_threshold,rename_sample,num_factors,pca_gene_num,
              unique_genes,platform,duplicate_probe,has_annotation_gene_id,group_fc]

PARAMETERS = dict()
for parameter in PARAMETERS_list:
    PARAMETERS[parameter.name]=parameter    

##PARAMETERS = {
##    'preprocess': ['rma', 'mas5', 'loess', 'illumina_controls',
##                   'illumina', 'agilent', 'unknown_preprocess'],
##    'gene_center': ['mean', 'median', 'no_gene_center'],
##    'gene_normalize': ['variance', 'sum_of_squares', 'no_gene_normalize'],
##    'quantile': ['yes_quantile', 'no_quantile'],
##    'dwd': ['yes_dwd', 'no_dwd'],
##    'bfrm': ['yes_bfrm', 'no_bfrm'],
##    'shiftscale': ['yes_shiftscale', 'no_shiftscale'],
##    'combat': ['yes_combat', 'no_combat'],
##    'gene_order': ['t_test_p', 't_test_fdr', 'by_gene_list', 'no_order',
##                   'by_class_neighbors'],
##    'predataset': ['yes_predataset', 'no_predataset'],
##    'filter': 'integer',
##    'has_missing_value': ['no_missing', 'median_fill', 'zero_fill',
##                          'unknown_missing'],
##    'format': ['tdf', 'gct'],
##    'is_logged': ['no_logged', 'logged'],
##    'contents': 'list',
##    'ill_manifest': 'string',
##    'ill_chip': 'string',
##    'ill_bg_mode': ['ill_yes', 'ill_no'],
##    'ill_coll_mode': ['ill_none', 'ill_max', 'ill_median'],
##    'ill_clm': 'string',
##    'ill_custom_chip': 'string',
##    'ill_custom_manifest': 'string',
##    'cn_num_neighbors': 'integer',
##    'cn_num_perm': 'integer',
##    'cn_user_pval': 'float',
##    'cn_mean_or_median': ['cn_mean', 'cn_median'],
##    'cn_ttest_or_snr': ['cn_test', 'cn_snr'],
##    'cn_filter_data': ['cn_yes', 'cn_no'],
##    'cn_min_threshold': 'float',
##    'cn_max_threshold': 'float',
##    'cn_min_folddiff': 'float',
##    'cn_abs_diff': 'float',
##    'gene_select_threshold': 'float',
##    'rename_sample': ['yes_rename', 'no_rename'],
##    'num_factors': 'integer',
##    'pca_gene_num': 'integer',
##    'unique_genes': ['average_genes', 'high_var', 'first_gene',
##                     'no_unique_genes'],
##    'platform': ["'HG_U133_Plus_2'", "'HG_U133B'", "'HG_U133A'",
##                 "'HG_U133A_2'", "'HG_U95A'", "'HumanHT_12'","'HumanWG_6'", "'HG_U95Av2'",
##                 "'Entrez_ID_human'", "'Entrez_symbol_human'", "'Hu6800'",
##                 "'Mouse430A_2'", "'MG_U74Cv2'", "'Mu11KsubB'", "'Mu11KsubA'",
##                 "'MG_U74Av2'", "'Mouse430_2'", "'MG_U74Bv2'",
##                 "'Entrez_ID_mouse'", "'MouseRef_8'", "'Entrez_symbol_mouse'",
##                 "'RG_U34A'", "'RAE230A'", 'unknown_platform'],
##    'duplicate_probe': ['yes_duplicate_probe', 'high_var_probe',
##                        'closest_probe'],
##    'duplicate_data': ['yes_duplicate_data', 'no_duplicate_data'],
##    'has_annotation_gene_id':['yes_gene_id','no_gene_id'],
##    'group_fc':'integer'}
##
##
##DEFAULT = {
##    'gene_center': 'no_gene_center',
##    'gene_normalize': 'no_gene_normalize',
##    'filter': '0', 'predataset': 'no_predataset',
##    'is_logged': 'logged', 'rename_sample': 'no_rename',
##    'gene_order': 'no_order', 'format': 'tdf',
##    'duplicate_probe': 'yes_duplicate_probe',
##    'duplicate_data': 'yes_duplicate_data',
##    'unique_genes': 'no_unique_genes'}
##

##common_parameters=['preprocess','gene_center','gene_normalize',
##                   'has_missing_value','format','is_logged', 'unique_genes',
##                   'filter', 'gene_order','predataset','num_features','group_fc']
##
##normalize_parameters = ['quantile','dwd','bfrm','shiftscale','combat']
##
##optional_parameters = ['contents','gene_select_threshold', 'rename_sample',
##                       'num_factors', 'pca_gene_num', 'platform','duplicate_probe',
##                       'duplicate_data', 'has_annotation_gene_id']
