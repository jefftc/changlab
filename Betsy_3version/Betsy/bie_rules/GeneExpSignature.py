#GeneExpSignature
from Betsy.bie3 import *
import Database
import GeneExpProcessing
import PcaAnalysis
import BasicDataTypes

GenesetAnalysis=DataType(
    'GenesetAnalysis',
    AttributeDef("contents",Database.CONTENTS,
                               'unspecified','unspecified',help="contents"),
    AttributeDef("allgenes",['yes','no'], "no","no",help="analyze all geneset"),
    AttributeDef("automatch",['yes','no'], "no","no",help="automatch in geneset analysis"),
    help="Geneset analysis file"
    
    )
GenesetPlot=DataType(
    'GenesetPlot',
    AttributeDef("contents",Database.CONTENTS,
                               'unspecified','unspecified',help="contents"),
    AttributeDef("allgenes",['yes','no'], "no","no",help="analyze all geneset"),
    AttributeDef("automatch",['yes','no'], "no","no",help="automatch in geneset analysis"),
    help="Geneset plot file"
    )

GenesetFile = DataType('GenesetFile',
                       AttributeDef("contents",Database.CONTENTS,
                               'unspecified','unspecified',help="contents"),)

SignatureScore = DataType(
    'SignatureScore',AttributeDef("contents",Database.CONTENTS,
                                  'unspecified','unspecified',help="contents"),
    help="Signature score file")

GenesetReportFile = DataType(
    'GenesetReportFile',
    help="Report file for gene set report"
    )
list_files = [GenesetAnalysis,GenesetPlot,GenesetFile,SignatureScore,GenesetReportFile]

all_modules = [
    Module(
        'score_pathway_with_geneset',
        [GenesetFile,
         GeneExpProcessing.SignalFile],GenesetAnalysis,
         OptionDef("geneset_value",help="geneset to score pathway"),
         #Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("gene_center",MUST_BE,'mean',1),
         Constraint("gene_normalize",MUST_BE,'variance',1),
         Constraint("annotate",MUST_BE,'yes',1),
         Constraint("unique_genes",MUST_BE,'high_var',1),
         Consequence("allgenes",SET_TO,'no'),
         Consequence("automatch",SET_TO_ONE_OF,['yes','no']),
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         help="score pathway with geneset"
        ),
    Module(
        'score_pathway_with_geneset_all',
        [GenesetFile,
         GeneExpProcessing.SignalFile],GenesetAnalysis,
         #Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("gene_center",MUST_BE,'mean',1),
         Constraint("gene_normalize",MUST_BE,'variance',1),
         Constraint("annotate",MUST_BE,'yes',1),
         Constraint("unique_genes",MUST_BE,'high_var',1),
         Consequence("allgenes",SET_TO,'yes'),
         Consequence("automatch",SET_TO_ONE_OF,['yes','no']),
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS,0),
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
        Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="plot geneset analysis file to bar plot"
        ),
    Module(
        'score_pathway_with_scoresig',
        [GeneExpProcessing.SignalFile,GeneExpProcessing.SignalFile],SignatureScore,
        OptionDef('platform_value','HG_U133A',help="platform to add"),
        Constraint("format",MUST_BE,'tdf',0),
        Constraint("preprocess",MUST_BE,'rma',0),
        Constraint("quantile_norm",MUST_BE,'yes',0),
        Constraint("logged",MUST_BE,'yes',0),
        Constraint("platform",MUST_BE,"yes",0),
        Constraint("duplicate_probe",MUST_BE,'high_var_probe',0),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("preprocess",MUST_BE,'mas5',1),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("platform",MUST_BE,'yes',1),
        Constraint("duplicate_probe",MUST_BE,'high_var_probe',1),
        Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Consequence('contents',SAME_AS_CONSTRAINT,0),
        help="score pathway iwth scoresig method"
        ),
    
    Module(
        'make_geneset_report',
        [GenesetAnalysis,
         GenesetPlot],GenesetReportFile,
        help="make geneset report"),
    ]
