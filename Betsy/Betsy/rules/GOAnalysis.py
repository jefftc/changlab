#GOAnalysis
from Betsy.bie3 import *
import Database
import BasicDataTypes

DavidFile = DataType(
    'DavidFile',AttributeDef(
        "contents",
        Database.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="David result file")
GatherFile = DataType(
    'GatherFile',
    AttributeDef("network",['yes','no'],'no','no',help="network option in Gather tool"),
    AttributeDef('homologs',['yes','no'],'no','no',help="homologs option in Gather tool"),
    AttributeDef('annot_type',['kegg','gene_ontology','transfac'],'kegg','kegg',
                 help="annot type option in Gather tool"),
    AttributeDef("contents", Database.CONTENTS,
                               'unspecified','unspecified',help='contents'),
    AttributeDef('gene_order',['no', "gene_list", "class_neighbors",
                               "t_test_p", "t_test_fdr",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change'],
                 't_test_p',"t_test_p",help="gene order method for gene list"),
    help="result file using Gather tool")
    
list_files = [DavidFile,GatherFile]

all_modules = [
    Module(
        'annotate_genes_with_david',
         BasicDataTypes.GeneListFile,
         DavidFile,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help="annotate genes with david tool"),
    Module(
        'annotate_genes_with_gather',
         BasicDataTypes.GeneListFile,GatherFile,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Constraint("gene_order",CAN_BE_ANY_OF,['no', "gene_list", "class_neighbors",
                               "t_test_p", "t_test_fdr",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change']),
         Consequence("gene_order",SAME_AS_CONSTRAINT),
         help="annotate genes using Gather tool"
         )]
