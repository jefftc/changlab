#GOAnalysis
from Betsy.bie3 import *
import BasicDataTypes as BDT

DavidFile = DataType(
    'DavidFile',
    #AttributeDef(
    #    "contents", BDT.CONTENTS, "unspecified", "unspecified",
    #    help='contents'),
    help="David result file")

GatherFile = DataType(
    'GatherFile',
    AttributeDef(
        "network", ['yes', 'no'], 'no', 'no',
        help="network option in Gather tool",
        ),
    AttributeDef(
        'homologs', ['yes', 'no'], 'no', 'no',
        help="homologs option in Gather tool",
        ),
    AttributeDef(
        'annot_type', ['kegg', 'gene_ontology', 'transfac'],
        'gene_ontology', 'gene_ontology',
        help="annot type option in Gather tool",
        ),
    #AttributeDef(
    #    "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
    #    help='contents'),
    #AttributeDef(
    #    'gene_order', [
    #        'no', "gene_list", "class_neighbors", "t_test_p", "t_test_fdr",
    #        'diff_ttest', 'diff_sam', 'diff_ebayes', 'diff_fold_change'],
    #    'diff_ttest', "diff_ttest",
    #    help="gene order method for gene list"),
    help="result file using Gather tool")

all_data_types = [
    DavidFile,
    GatherFile,
    ]

all_modules = [
    ModuleNode(
        'annotate_genes_with_david',
        BDT.GeneListFile, DavidFile,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        help="annotate genes with david tool",
        ),
    
    ModuleNode(
        "annotate_genes_with_gather",
        BDT.GeneListFile, GatherFile,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("filtered", MUST_BE, "yes"),
        #Constraint(
        #    "gene_order", CAN_BE_ANY_OF,
        #    ['no', "gene_list", "class_neighbors", "t_test_p", "t_test_fdr",
        #     'diff_ttest', 'diff_sam', 'diff_ebayes', 'diff_fold_change']),
        #Consequence("gene_order", SAME_AS_CONSTRAINT),
        help="annotate genes using Gather tool",
        )
    ]
