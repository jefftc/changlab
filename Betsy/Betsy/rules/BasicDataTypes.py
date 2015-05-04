#BasicDataTypes

from Betsy.bie3 import *

CONTENTS = [
    "train0", "train1","test", "class0,class1,test",
    "class0", "class1", "class0,class1","unspecified",
    "diff_class0", "diff_class1", "diff_class0,diff_class1",
    "diff_unspecified",
    ]

RenameFile = DataType(
    'RenameFile',
    AttributeDef(
        "contents", CONTENTS, 'unspecified', 'unspecified', help="contents"),
    AttributeDef(
        "labels_from", ["title", "description"], 'title', 'title',
        help="labels from title or description"),
    help="A file used to rename the sample name in the gene expression file.")

ExpressionFiles = DataType(
    "ExpressionFiles",
    AttributeDef(
        "contents", CONTENTS, 'unspecified', 'unspecified', help="contents"),
    AttributeDef(
        "filetype", ['unknown', 'cel', 'gpr', 'idat', 'agilent', 'matrix'],
        'unknown', 'unknown', help="filetype"),
    help="Expression file folder, can be CELFiles, IDATFiles, "\
    "AgilentFile, GPRFiles")


GeneListFile=DataType(
    "GeneListFile", 
    AttributeDef(
        'cn_mean_or_median', ['mean', 'median'], 'mean', 'mean',
        help="class neighbors mean or median"),
    AttributeDef(
        'cn_ttest_or_snr',['t_test', 'snr'], 't_test', 't_test',
        help="class neighbors ttest or snr"),
    AttributeDef(
        'cn_filter_data',['yes','no'], 'no', 'no',
        help="class neighbors filter data or not"),
    AttributeDef(
        'gene_order', ['no', "gene_list", "class_neighbors", "t_test_p",
                       "t_test_fdr",'diff_ttest','diff_sam',
                       'diff_ebayes','diff_fold_change'],
        't_test_p', "t_test_p", help="gene order method"),
    AttributeDef(
        "contents", CONTENTS, 'unspecified', 'unspecified', help="contents"),
    help="A file that contains a list of genes.")


list_files = [RenameFile, ExpressionFiles, GeneListFile]

all_modules = [
    ]
