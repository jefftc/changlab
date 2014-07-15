#GenesetAnalysis
from Betsy.bie3 import *
import SignalFile_rule
EMTAnalysis=DataType('EMTAnalysis',
                     AttributeDef("contents",SignalFile_rule.CONTENTS,
                                  'unspecified','unspecified'))
    
CellTypeFile=DataType('CellTypeFile',
                     AttributeDef("contents",SignalFile_rule.CONTENTS,
                                  'unspecified','unspecified'))
list_files = [EMTAnalysis,CellTypeFile]

all_modules = [
    Module(
        'analyze_phenotype_for_EMT',
         [SignalFile_rule.SignalFile,CellTypeFile],EMTAnalysis,
         UserInputDef("geneset_value"),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence('contents',SAME_AS_CONSTRAINT,0)),

    ]
