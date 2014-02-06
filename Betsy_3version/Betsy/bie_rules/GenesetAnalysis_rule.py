#GenesetAnalysis
from Betsy.bie3 import *
import SignalFile2_rule
GenesetAnalysis=DataType(
    'GenesetAnalysis',
    AttributeDef("geneset",["yes","no"],"no","no"),
    AttributeDef("allgenes",['yes','no'], "no","no"),
    AttributeDef("automatch",['yes','no'], "no","no"),
    )
GenesetPlot=DataType(
    'GenesetPlot',
    AttributeDef("geneset",["yes","no"],"no","no"),
    AttributeDef("allgenes",['yes','no'], "no","no"),
    AttributeDef("automatch",['yes','no'], "no","no"),
    )

GenesetFile = DataType('GenesetFile')
list_files = [GenesetAnalysis,GenesetPlot,GenesetFile]

all_modules = [
    Module(
        'score_pathway_with_geneset',
        [GenesetFile,
         SignalFile2_rule.SignalFile2],GenesetAnalysis,
         UserInputDef("geneset_value"),
         Constraint("logged",MUST_BE,'yes',1),
         Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("gene_center",MUST_BE,'mean',1),
         Constraint("gene_normalize",MUST_BE,'variance',1),
         Constraint("unique_genes",MUST_BE,'high_var',1),
         Consequence("geneset",SET_TO,"yes"),
         Consequence("allgenes",SET_TO_ONE_OF,['yes','no']),
         Consequence("automatch",SET_TO_ONE_OF,['yes','no'])),
    Module(
        'plot_geneset_score_bar',
        GenesetAnalysis,GenesetPlot,
        UserInputDef("geneset_value"),
        Consequence("geneset",SET_TO,"yes"),
        Consequence("allgenes",SET_TO_ONE_OF,['yes','no']),
        Consequence("automatch",SET_TO_ONE_OF,['yes','no'])),
    ]
