from Betsy.bie3 import *

CONTENTS2NUMCLASSES = {
    "unspecified" : None,   # can be any number of classes
    
    # For 2 class problems.
    "class0" : 1,
    "class1" : 1,
    "class0,class1" : 2,
    "class0,class1,test" : 3,

    # Not sure that this is needed.
    "train0" : 1,
    "train1" : 1,
    "test" : 1,
    
    # Why is diff_unspecified needed?  Is this needed at all?
    "diff_unspecified" : 1,
    "diff_class0" : 1,
    "diff_class1" : 1,
    "diff_class0,diff_class1" : 2,
    }

CONTENTS = CONTENTS2NUMCLASSES.keys()


## PREPROCESS1 = [
##     "unknown", "illumina", "agilent", "mas5", "rma", "loess",
##     "rsem", 'RSEM_genes', 'RSEM_exons', 'humanmethylation450', 'mirnaseq',
##     'rppa', 'clinical', 'affymetrix']
## # What is this for?
## # TODO: get rid of this
## PREPROCESS_WOrma = [
##     "unknown", "illumina", "agilent", "mas5", "loess",
##     "rsem", 'RSEM_genes', 'RSEM_exons', 'humanmethylation450', 'mirnaseq',
##     'rppa', 'clinical', 'affymetrix']
## PREPROCESS = PREPROCESS1 + ['any']

PREPROCESS = [
    "unknown",

    # Gene Expression
    "mas5",
    "rma",
    "agilent",
    "illumina",
    "loess",  # What is this for?  Better in normalization?

    # RNA-Seq
    #"RSEM",
    "tpm",
    "fpkm",
    "counts",
    "cpm",
    
    # Not sure why this is here.
    #'RSEM_genes',
    #'RSEM_exons',

    # Not sure why there are here.
    #'humanmethylation450',
    #'mirnaseq',
    #'rppa',
    #'clinical',
    #'affymetrix',
    ]
ANY_PREPROCESS = PREPROCESS + ["any"]


GENE_ORDER = [
    "none",
    "user_defined",
    "fold_change",
    "ttest_p",
    "ttest_fdr",
    "sam_p",
    "ebayes_fdr",
    "class_neighbors",
    ]
GENE_ORDER_not_none = [x for x in GENE_ORDER if x != "none"]

YESNO = ["yes", "no"]



RenameFile = DataType(
    'RenameFile',
    AttributeDef(
        "contents", CONTENTS, 'unspecified', 'unspecified', help="contents"),
    #AttributeDef(
    #    "labels_from", ["title", "description"], 'title', 'title',
    #    help="labels from title or description"),
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
        "filtered", YESNO, "no", "no",
        help="Whether this gene list has been filtered."),
    AttributeDef(
        "gene_order", GENE_ORDER, "none", "none",
        help="How this gene list is ordered.",
        ),
    AttributeDef(
        "preprocess", PREPROCESS, "unknown", "unknown",
        help="preprocess method"),
    AttributeDef(
        "contents", CONTENTS, 'unspecified', 'unspecified', help="contents"),
    #AttributeDef(
    #    'cn_mean_or_median', ['mean', 'median'], 'mean', 'mean',
    #    help="class neighbors mean or median"),
    #AttributeDef(
    #    'cn_ttest_or_snr', ['t_test', 'snr'], 't_test', 't_test',
    #    help="class neighbors ttest or snr"),
    #AttributeDef(
    #    'cn_filter_data', ['yes','no'], 'no', 'no',
    #    help="class neighbors filter data or not"),
    help="A file that contains a list of genes.",
    )



all_data_types = [
    RenameFile,
    GeneListFile,
    ExpressionFiles,
    ]

all_modules = [
    ]
