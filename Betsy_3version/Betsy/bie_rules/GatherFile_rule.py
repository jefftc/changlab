#GatherFile
from Betsy.bie3 import *
import SignalFile_rule
GatherFile = DataType(
    'GatherFile',
    AttributeDef("network",['yes','no'],'no','no'),
    AttributeDef('homologs',['yes','no'],'no','no'),
    AttributeDef('annot_type',['kegg','gene_ontology','transfac'],'kegg','kegg'))
list_files = [GatherFile]
all_modules = [
        Module(
        'annotate_genes_with_gather',
         SignalFile_rule.GeneListFile,
         GatherFile)]
