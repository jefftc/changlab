from Betsy.bie3 import *
import BasicDataTypes as BDT
import GeneExpProcessing


ClinicalAnalysis = DataType(
    "ClinicalAnalysis",
    AttributeDef("contents",BDT.CONTENTS,
                               'unspecified','unspecified',help="contents"),
    help="Clincal analysis file")
ClinicalFile = DataType(
    "ClinicalFile",
    AttributeDef("contents",BDT.CONTENTS,
                               'unspecified','unspecified',help="contents"),
    help="clinical file for clinical analysis")

EMTAnalysis=DataType('EMTAnalysis',
                     AttributeDef("contents",BDT.CONTENTS,
                                  'unspecified','unspecified',help="contents"),
                     help="EMT analysis result file")
    
CellTypeFile=DataType('CellTypeFile',
                     AttributeDef("contents",BDT.CONTENTS,
                                  'unspecified','unspecified',help='contents'),
                      help="Cell type file for EMT analysis")

list_files = [ClinicalAnalysis,ClinicalFile,EMTAnalysis,CellTypeFile]

all_modules = [
    ModuleNode(
        'analyze_clinical_outcome',
        [GeneExpProcessing.SignalFile,ClinicalFile],
        ClinicalAnalysis,
        OptionDef("outcome",help="outcome header"),
        OptionDef("dead",help="dead header"),
        OptionDef("genename",help="gene name to analyze"),
        OptionDef("rank_cutoff",0.5,help="number of rank cutoff"),
        OptionDef("zscore_cutoff",'n1,1',help="number of zscore cutoff"),
        Constraint("logged",MUST_BE,'yes',0),
        Constraint("contents",CAN_BE_ANY_OF,BDT.CONTENTS,0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents",SAME_AS_CONSTRAINT, 0),
        help="analyze clincal outcome"
        ),
    ModuleNode(
        'analyze_phenotype_for_EMT',
         [GeneExpProcessing.SignalFile,CellTypeFile],EMTAnalysis,
         OptionDef("geneset_value",help="geneset value for EMT analysis"),
         Constraint("contents",CAN_BE_ANY_OF,BDT.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence('contents',SAME_AS_CONSTRAINT,0),
         help="analyze phenotype for EMT"),
    ]
