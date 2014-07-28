#GenesetAnalysis
from Betsy.bie3 import *
import SignalFile_rule
GenesetAnalysis=DataType(
    'GenesetAnalysis',
    AttributeDef("allgenes",['yes','no'], "no","no"),
    AttributeDef("automatch",['yes','no'], "no","no"),
    
    )
GenesetPlot=DataType(
    'GenesetPlot',
    AttributeDef("allgenes",['yes','no'], "no","no"),
    AttributeDef("automatch",['yes','no'], "no","no"),
    )

GenesetFile = DataType('GenesetFile')
list_files = [GenesetAnalysis,GenesetPlot,GenesetFile]

all_modules = [
    Module(
        'score_pathway_with_geneset',
        [GenesetFile,
         SignalFile_rule.SignalFile],GenesetAnalysis,
         UserInputDef("geneset_value"),
         #Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("gene_center",MUST_BE,'mean',1),
         Constraint("gene_normalize",MUST_BE,'variance',1),
         Constraint("annotate",MUST_BE,'yes',1),
         Constraint("unique_genes",MUST_BE,'high_var',1),
         Consequence("allgenes",SET_TO,'no'),
         Consequence("automatch",SET_TO_ONE_OF,['yes','no']),
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
         ),
    
    Module(
        'plot_geneset_score_bar',
        GenesetAnalysis,GenesetPlot,
        Constraint("allgenes",CAN_BE_ANY_OF,['yes','no']),
        Constraint("automatch",CAN_BE_ANY_OF,['yes','no']),
        Consequence("allgenes",SAME_AS_CONSTRAINT),
        Consequence("automatch",SAME_AS_CONSTRAINT),
        ),
    ]
