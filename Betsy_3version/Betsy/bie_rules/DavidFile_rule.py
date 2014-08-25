#DavidFile
from Betsy.bie3 import *
import SignalFile_rule
DavidFile = DataType(
    'DavidFile',AttributeDef(
        "contents",
        SignalFile_rule.CONTENTS,
        "unspecified", "unspecified",help='contents'),
    help="David result file")

list_files = [DavidFile]

all_modules = [
    Module(
        'annotate_genes_with_david',
         SignalFile_rule.GeneListFile,
         DavidFile,
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help="annotate genes with david tool")]
