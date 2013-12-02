#DavidFile
from Betsy import bie
import SignalFile_rule
DavidFile = bie.DataType(
    'DavidFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True))

list_files = [DavidFile]

all_modules = [
    bie.Module(
        'annotate_genes_with_david',
         SignalFile_rule.GeneListFile,
         DavidFile)]
