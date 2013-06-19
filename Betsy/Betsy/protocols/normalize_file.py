#normalize_file.py
from Betsy.protocol_utils import Parameter

PRETTY_NAME="Normalize gene expression."
COMMON = 'Common Parameters'
NORMALIZE = 'Normalize Parameters'
OPTIONAL = 'Optional Parameters'
ILLUMINA = 'Illumina Normalize Parameters'


CATEGORIES=[COMMON,NORMALIZE,OPTIONAL,ILLUMINA]

#input predicates
INPUTS = [
    'gse_id',
    'gse_id_and_platform',
    'cel_files',
    'gpr_files',
    'idat_files',
    'agilent_files',
    'input_signal_file',
    ]

#output predicates
OUTPUTS = 'make_normalize_report'

# <input predicate name> : information for converting to fact
predicate2arguments = {
    'gse_id': ([], '[]'),
    'gse_id_and_platform': ([], '[]'),
    'idat_files': (['version', 'unknown_version'], '[]'),
    'gpr_files': (['version', 'unknown_version'], '[]'),
    'cel_files': (['version', 'unknown_version'], '[]'),
    'agilent_files': (['version', 'unknown_version'], '[]'),
    'input_signal_file': (['status', 'given'], '[]'),
    }

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
            Parameter('bfrm',pretty_name='BFRM',default_value='no_bfrm',
                 choices=['yes_bfrm', 'no_bfrm'],category=NORMALIZE),
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
            Parameter('gene_select_threshold',pretty_name='Gene Selection Threshold',
                                  type= 'float',category=OPTIONAL),
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
                                   choices=['yes_gene_id','no_gene_id'],category=OPTIONAL)]

