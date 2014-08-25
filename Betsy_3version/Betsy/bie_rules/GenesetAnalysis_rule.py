#GenesetAnalysis
from Betsy.bie3 import *
import SignalFile_rule
GenesetAnalysis=DataType(
    'GenesetAnalysis',
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                               'unspecified','unspecified',help="contents"),
    AttributeDef("allgenes",['yes','no'], "no","no",help="analyze all geneset"),
    AttributeDef("automatch",['yes','no'], "no","no",help="automatch in geneset analysis"),
    help="Geneset analysis file"
    
    )
GenesetPlot=DataType(
    'GenesetPlot',
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                               'unspecified','unspecified',help="contents"),
    AttributeDef("allgenes",['yes','no'], "no","no",help="analyze all geneset"),
    AttributeDef("automatch",['yes','no'], "no","no",help="automatch in geneset analysis"),
    help="Geneset plot file"
    )

GenesetFile = DataType('GenesetFile',
                       AttributeDef("contents",SignalFile_rule.CONTENTS,
                               'unspecified','unspecified',help="contents"),)
list_files = [GenesetAnalysis,GenesetPlot,GenesetFile]

all_modules = [
    Module(
        'score_pathway_with_geneset',
        [GenesetFile,
         SignalFile_rule.SignalFile],GenesetAnalysis,
         OptionDef("geneset_value",help="geneset to score pathway"),
         #Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("gene_center",MUST_BE,'mean',1),
         Constraint("gene_normalize",MUST_BE,'variance',1),
         Constraint("annotate",MUST_BE,'yes',1),
         Constraint("unique_genes",MUST_BE,'high_var',1),
         Consequence("allgenes",SET_TO,'no'),
         Consequence("automatch",SET_TO_ONE_OF,['yes','no']),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         help="score pathway with geneset"
        ),
    Module(
        'score_pathway_with_geneset_all',
        [GenesetFile,
         SignalFile_rule.SignalFile],GenesetAnalysis,
         #Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("gene_center",MUST_BE,'mean',1),
         Constraint("gene_normalize",MUST_BE,'variance',1),
         Constraint("annotate",MUST_BE,'yes',1),
         Constraint("unique_genes",MUST_BE,'high_var',1),
         Consequence("allgenes",SET_TO,'yes'),
         Consequence("automatch",SET_TO_ONE_OF,['yes','no']),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         help="score pathway with all the geneset"
         ),
    
    Module(
        'plot_geneset_score_bar',
        GenesetAnalysis,GenesetPlot,
        Constraint("allgenes",CAN_BE_ANY_OF,['yes','no']),
        Constraint("automatch",CAN_BE_ANY_OF,['yes','no']),
        Consequence("allgenes",SAME_AS_CONSTRAINT),
        Consequence("automatch",SAME_AS_CONSTRAINT),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="plot geneset analysis file to bar plot"
        ),
    ]
