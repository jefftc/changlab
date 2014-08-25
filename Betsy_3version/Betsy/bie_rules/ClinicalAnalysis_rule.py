#ClinicalAnalysis_rule.py
from Betsy.bie3 import *
import SignalFile_rule
ClinicalAnalysis = DataType(
    "ClinicalAnalysis",
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified',help="contents"),
    help="Clincal analysis file")
ClinicalFile = DataType(
    "ClinicalFile",
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified',help="contents"),
    help="clinical file for clinical analysis")


list_files = [ClinicalAnalysis,ClinicalFile]

all_modules = [
    Module(
        'analyze_clinical_outcome',
        [SignalFile_rule.SignalFile,ClinicalFile],
        ClinicalAnalysis,
        OptionDef("outcome",help="outcome header"),
        OptionDef("dead",help="dead header"),
        OptionDef("genename",help="gene name to analyze"),
        OptionDef("rank_cutoff",0.5,help="number of rank cutoff"),
        OptionDef("zscore_cutoff",'n1,1',help="number of zscore cutoff"),
        Constraint("logged",MUST_BE,'yes',0),
        Constraint("contents",CAN_BE_ANY_OF,["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents",SAME_AS_CONSTRAINT, 0),
        help="analyze clincal outcome"
        ),
    ]
