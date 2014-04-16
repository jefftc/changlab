#GatherFile
from Betsy.bie3 import *
import SignalFile_rule
GatherFile = DataType(
    'GatherFile',
    AttributeDef("network",['yes','no'],'no','no'),
    AttributeDef('homologs',['yes','no'],'no','no'),
    AttributeDef('annot_type',['kegg','gene_ontology','transfac'],'kegg','kegg'),
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified'))
list_files = [GatherFile]
all_modules = [
        Module(
        'annotate_genes_with_gather',
         SignalFile_rule.GeneListFile,GatherFile,
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
         )]
