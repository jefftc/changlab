#batch_effect_remove.py
from Betsy.protocol_utils import Parameter

PRETTY_NAME="Remove batch effects."
COMMON = 'Common Parameters'
OPTIONAL = 'Optional Parameters'
ILLUMINA = 'Illumina Normalize Parameters'
CLASS_NEIGHBORS='Class Neighbor Parameters'
CATEGORIES=[COMMON,OPTIONAL,ILLUMINA,CLASS_NEIGHBORS]
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
    'rename_list_file',]

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
    'input_signal_file': (['status', 'given'], '[]'),
    'rename_list_file': ([], '[]')}

#parameter objects
PARAMETERS=[Parameter('preprocess',pretty_name='Preprocess',
                         choices=['rma', 'mas5',
                        'loess', 'illumina_controls',
                        'illumina', 'agilent',
                        'unknown_preprocess'],category=COMMON),     
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
            Parameter('rename_sample',pretty_name='Rename Sample',default_value='no_rename',
                          choices=['yes_rename', 'no_rename'], category=OPTIONAL),
            Parameter('num_factors',pretty_name='Number of Factors',type='integer',
                           category=OPTIONAL)]
  
