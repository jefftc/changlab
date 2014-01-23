#DavidFile
from Betsy.bie3 import *
import SignalFile_rule
DavidFile = DataType(
    'DavidFile')

list_files = [DavidFile]

all_modules = [
    Module(
        'annotate_genes_with_david',
         SignalFile_rule.GeneListFile,
         DavidFile)]
