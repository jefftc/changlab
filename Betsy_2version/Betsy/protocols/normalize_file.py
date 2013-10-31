#normalize_file.py
from Betsy.protocol_utils import Parameter

PRETTY_NAME="Normalize gene expression."
COMMON = 'Common Parameters'
NORMALIZE = 'Normalize Parameters'
OPTIONAL = 'Optional Parameters'
ILLUMINA = 'Illumina Normalize Parameters'


CATEGORIES=[COMMON,NORMALIZE,OPTIONAL,ILLUMINA]

#input  datatype
INPUTS = [
    'AgilentFiles',
    'ExpressionFiles',
    'RenameFile',
    'GPRFiles',
    'GEOSeries',
    'IDATFiles',
    'CELFiles',
    'SignalFile',
    'ClassLabelFile',
    'RenameFile',
    'GeneListFile'
    ]

#output datatype
OUTPUTS = 'ReportFile'

#parameter objects
PARAMETERS=[Parameter('preprocess',pretty_name='Preprocess',
                         choices=['rma', 'mas5',
                        'loess', 'illumina_controls',
                        'illumina', 'agilent','rsem',
                        'unknown'],category=COMMON,
                         description='preprocess method'),
            Parameter('gene_center',pretty_name='Gene Center', 
                        choices=['mean', 'median', 'no'],
                        category=COMMON,description='gene center method'),
            Parameter('gene_normalize',pretty_name='Gene Normalize',
                           choices=['variance', 'sum_of_squares', 'no'],
                           category=COMMON,description='gene normalize method'),
            Parameter('quantile_norm',pretty_name='Quantile', 
                     choices=['yes', 'no'],category=NORMALIZE,
                      description='normalize by quanitle'),
            Parameter('combat_norm',pretty_name='Combat', 
                     choices=['yes', 'no'],category=NORMALIZE,
                      description='normalize by combat'),
            Parameter('shiftscale_norm',pretty_name='Shiftscale', 
                     choices=['yes', 'no'],category=NORMALIZE,
                      description='normalize by shiftscale'),
            Parameter('dwd_norm',pretty_name='Dwd', 
                     choices=['yes', 'no'],category=NORMALIZE,
                      description='normalize by dwd'),
            Parameter('bfrm_norm',pretty_name='BFRM',
                 choices=['yes', 'no'],category=NORMALIZE,
                      description='normalize by bfrm'),
            Parameter('gene_order',pretty_name='Gene Order',
                 choices=['no', 't_test_p','t_test_fdr','by_class_neighbors'],category=NORMALIZE,
                      description='gene order'),
            Parameter('predataset',pretty_name='Predataset Process',
                       choices=['yes', 'no'],category=COMMON,
                      description='filter and threshold genes'),
            Parameter('filter',pretty_name='Filter',
                      category=COMMON,description='percentage to filter the missing value(0-100)'),
            Parameter('missing_values',pretty_name='if has missing value',
                              choices=['no', 'yes','unknown'],category=COMMON,
                       description='has the missing value or not'),
            Parameter('missing_algorithm',pretty_name='How to fill missing value',
                      choices=['none','median_fill','zero_fill'],category=COMMON,
                      description='method to the fill the missing value'),
            Parameter('format',pretty_name='File Format',default_value='tdf',
                        choices=['tdf', 'gct'],category=COMMON,description='output file format'),
            Parameter('logged',pretty_name='Logged or Not',default_value='yes',
                      choices=['no', 'yes'],category=COMMON,
                      description='signal is logged or not'),
            Parameter('contents',pretty_name='Contents',category=OPTIONAL,
                      choices=["train0", "train1", "test", 'class0,class1,test',
                  "class0", "class1", "class0,class1",
                  "no"], description='output group information'),
            Parameter('num_features',pretty_name='Gene Number',
                         category=COMMON,description='select num of genes'),
            Parameter('ill_manifest', pretty_name='Illumina Manifest File',
                      choices=[
                    'HumanHT-12_V3_0_R2_11283641_A.txt',
                    'HumanHT-12_V4_0_R2_15002873_B.txt',
                    'HumanHT-12_V3_0_R3_11283641_A.txt',
                    'HumanHT-12_V4_0_R1_15002873_B.txt',
                    'HumanMI_V1_R2_XS0000122-MAP.txt',
                    'HumanMI_V2_R0_XS0000124-MAP.txt',
                    'HumanRef-8_V2_0_R4_11223162_A.txt',
                    'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
                    'HumanRef-8_V3_0_R2_11282963_A.txt',
                    'HumanRef-8_V3_0_R3_11282963_A.txt',
                    'HumanWG-6_V2_0_R4_11223189_A.txt',
                    'HumanWG-6_V3_0_R2_11282955_A.txt',
                    'HumanWG-6_V3_0_R3_11282955_A.txt',
                    'MouseMI_V1_R2_XS0000127-MAP.txt',
                    'MouseMI_V2_R0_XS0000129-MAP.txt',
                    'MouseRef-8_V1_1_R4_11234312_A.txt',
                    'MouseRef-8_V2_0_R2_11278551_A.txt',
                    'MouseRef-8_V2_0_R3_11278551_A.txt',
                    'MouseWG-6_V1_1_R4_11234304_A.txt',
                    'MouseWG-6_V2_0_R2_11278593_A.txt',
                    'MouseWG-6_V2_0_R3_11278593_A.txt',
                    'RatRef-12_V1_0_R5_11222119_A.txt'
                    ],category=ILLUMINA,
                      description='Illumina manifest file in tab-delimited (TXT) format'),
            Parameter('ill_chip', pretty_name='Illumina chip File',choices=[
                    'ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
                    'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
                    'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
                    'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
                    'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
                    'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
                    'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
                    'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
                    'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
                    'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
                    'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
                    'ilmn_RatRef_12_V1_0_R5_11222119_A.chip'
                    ], category=ILLUMINA,description='CHIP file to map probes to genes.'),
            Parameter('ill_clm', pretty_name='Illumina clm File',
                     category=ILLUMINA,description='CLM file to map file names to sample names.'),
            Parameter('ill_custom_chip', pretty_name='Illumina Custom Chip File',
                             category=ILLUMINA,
                      description='Other CHIP file to map probes to genes.'),
            Parameter('ill_custom_manifest', pretty_name='Illumina Custom Manifest File',
                                 category=ILLUMINA,
                      description='Other Illumina manifest file in tab-delimited (TXT) format.'),
            Parameter('ill_bg_mode', pretty_name='Illumina Background Mode',
                         choices=['ill_yes', 'ill_no'],category=ILLUMINA,
                      description='Perform background subtraction.'),
            Parameter('ill_coll_mode', pretty_name='Illumina Coll Mode',
                        choices=['none', 'max', 'median'],category=ILLUMINA,
                      description='Collapse probes to genes based on the manifest or CHIP file (if provided).'),
            Parameter('gene_select_threshold',pretty_name='Gene Selection Threshold',
                        category=OPTIONAL,
                      description='gene select threshold'),
            Parameter('num_factors',pretty_name='Number of Factors',
                           category=OPTIONAL,description='num of factors for bfrm'),
            Parameter('pca_gene_num',pretty_name='PCA Gene Number',
                           category=OPTIONAL,description='gene number in PCA'),
            Parameter('unique_genes',pretty_name='Unique Genes',
                         choices=['average_genes', 'high_var', 'first_gene','no'],
                         category=COMMON,description='how to select unique genes'),
            Parameter('platform',pretty_name='Platform',choices=["'HG_U133_Plus_2'", "'HG_U133B'", "'HG_U133A'",
                 "'HG_U133A_2'", "'HG_U95A'", "'HumanHT_12'", "'HumanWG_6'","'HG_U95Av2'",
                 "'Entrez_ID_human'", "'Entrez_symbol_human'", "'Hu6800'",
                 "'Mouse430A_2'", "'MG_U74Cv2'", "'Mu11KsubB'", "'Mu11KsubA'",
                 "'MG_U74Av2'", "'Mouse430_2'", "'MG_U74Bv2'",
                 "'Entrez_ID_mouse'", "'MouseRef_8'", "'Entrez_symbol_mouse'",
                 "'RG_U34A'", "'RAE230A'", 'unknown_platform'],category=OPTIONAL,
                      description='output platform'),
            Parameter('duplicate_probe',pretty_name='Duplicate Probe',
                            choices=['yes', 'high_var_probe','closest_probe','no'],
                            category=OPTIONAL,description='how to remove duplicate probe'),
            Parameter('annotate',pretty_name='Annotation',
                                   choices=['yes','no'],category=OPTIONAL,
                      description='how to annotate the output file'),
            Parameter('report_type',pretty_name='normalize report',category=COMMON,
                      choices=['normalize'],
                      default_value='normalize',description='make normalize report')]

