from Betsy.bie3 import *
import BasicDataTypes as BDT
import SignalFile
#import PcaAnalysis

GenesetAnalysis=DataType(
    'GenesetAnalysis',
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    AttributeDef(
        "allgenes", ['yes', 'no'], "no", "no", help="analyze all geneset"),
    AttributeDef(
        "automatch", ['yes', 'no'], "no", "no",
        help="automatch in geneset analysis"),
    help="Geneset analysis file"
    
    )
GenesetPlot=DataType(
    'GenesetPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    AttributeDef(
        "allgenes", ['yes', 'no'], "no", "no", help="analyze all geneset"),
    AttributeDef(
        "automatch", ['yes', 'no'], "no", "no",
        help="automatch in geneset analysis"),
    help="Geneset plot file"
    )

GenesetFile = DataType(
    'GenesetFile',
    AttributeDef(
        "contents", BDT.CONTENTS,
        'unspecified', 'unspecified', help="contents"),
    )

SignatureScore = DataType(
    'SignatureScore',
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, 'unknown', 'any',
        #"preprocess", SignalFile.PREPROCESS, 'unknown', 'unknown',
        help="preprocess"),
    help="Output path containing results from scoresig analysis.",
    )

GenesetReportFile = DataType(
    'GenesetReportFile',
    help="Report file for gene set report"
    )

ScoreCompareReportFile = DataType(
    'ScoreCompareReportFile',
    help="Report file for score signature comparison report"
    )

ScoreComparePlot = DataType(
    'ScoreComparePlot',
    help="Plot file for score signature comparison"
    )

all_data_types = [GenesetAnalysis, GenesetPlot, GenesetFile, SignatureScore,
              GenesetReportFile, ScoreCompareReportFile, ScoreComparePlot]

all_modules = [
    ModuleNode(
        'score_pathway_with_geneset',
        [GenesetFile,
         SignalFile.SignalFile], GenesetAnalysis,
         OptionDef("geneset_value", help="geneset to score pathway"),
         #Constraint("quantile_norm", MUST_BE, 'yes', 1),
         Constraint("gene_center", MUST_BE, 'mean', 1),
         Constraint("gene_normalize", MUST_BE, 'variance', 1),
         Constraint("annotate", MUST_BE, 'yes', 1),
         Constraint("unique_genes", MUST_BE, 'high_var', 1),
         Consequence("allgenes", SET_TO, 'no'),
         Consequence("automatch", SET_TO_ONE_OF, ['yes', 'no']),
         Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
         Constraint("contents", SAME_AS, 0, 1),
         Consequence("contents", SAME_AS_CONSTRAINT, 0),
         help="score pathway with geneset"
        ),
    ModuleNode(
        'score_pathway_with_geneset_all',
        [GenesetFile,
         SignalFile.SignalFile], GenesetAnalysis,
         #Constraint("quantile_norm", MUST_BE, 'yes', 1),
         Constraint("gene_center", MUST_BE, 'mean', 1),
         Constraint("gene_normalize", MUST_BE, 'variance', 1),
         Constraint("annotate", MUST_BE, 'yes', 1),
         Constraint("unique_genes", MUST_BE, 'high_var', 1),
         Consequence("allgenes", SET_TO, 'yes'),
         Consequence("automatch", SET_TO_ONE_OF, ['yes', 'no']),
         Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
         Constraint("contents", SAME_AS, 0, 1),
         Consequence("contents", SAME_AS_CONSTRAINT, 0),
         help="score pathway with all the geneset"
         ),
    
    ModuleNode(
        'plot_geneset_score_bar',
        GenesetAnalysis, GenesetPlot,
        Constraint("allgenes", CAN_BE_ANY_OF, ['yes', 'no']),
        Constraint("automatch", CAN_BE_ANY_OF, ['yes', 'no']),
        Consequence("allgenes", SAME_AS_CONSTRAINT),
        Consequence("automatch", SAME_AS_CONSTRAINT),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot geneset analysis file to bar plot"
        ),
    ModuleNode(
        'score_pathway_with_scoresig_affymetrix',
        [SignalFile.SignalFile, SignalFile.SignalFile], SignatureScore,
        Constraint("format", MUST_BE, 'tdf', 0),
        Constraint("preprocess", MUST_BE, 'rma', 0),
        Constraint("quantile_norm", MUST_BE, 'yes', 0),
        Constraint("logged", MUST_BE, 'yes', 0),
        Constraint("platform", MUST_BE, "u133A", 0),
        Constraint("duplicate_probe", MUST_BE, 'high_var_probe', 0),
        Constraint("format", MUST_BE, 'tdf', 1),
        Constraint("preprocess", MUST_BE, 'mas5', 1),
        Constraint("logged", MUST_BE, 'yes', 1),
        Constraint("platform", SAME_AS, 0, 1),
        Constraint("duplicate_probe", MUST_BE, 'high_var_probe', 1),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence('contents', SAME_AS_CONSTRAINT, 0),
        Consequence('preprocess', SET_TO, 'rma'),
        help="score pathway iwth scoresig method for affymetrix"
        ),
     ModuleNode(
        'score_pathway_with_scoresig',
        SignalFile.SignalFile, SignatureScore,
        Constraint("format", MUST_BE, 'tdf'),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Constraint("quantile_norm", MUST_BE, 'yes'),
        Constraint("logged", MUST_BE, 'yes'),
        Constraint("duplicate_probe", MUST_BE, 'high_var_probe'),
        Constraint("platform", MUST_BE, "u133A", 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence('contents', SAME_AS_CONSTRAINT),
        Consequence('preprocess', SAME_AS_CONSTRAINT),
        help="score pathway iwth scoresig method"
        ),
    ModuleNode(
        'convert_SignatureScore_preprocess',
        SignatureScore, SignatureScore,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SET_TO, 'any'),
        help='convert SignatureScore preprocess from others to any',
        ),
    
##    ModuleNode(
##        'make_geneset_report',
##        [GenesetAnalysis,
##         GenesetPlot], GenesetReportFile,
##        help="make geneset report"),
    ModuleNode(
        'make_score_signature_report',
        [GenesetAnalysis, GenesetPlot, SignatureScore], GenesetReportFile,
        help="make signature report"),
    
    #ModuleNode(
    #    'compare_signature_predictions',
    #    [SignatureScore, SignatureScore, SignatureScore],
    #    ScoreCompareReportFile,
    #    #Constraint("preprocess", MUST_BE, 'RSEM_genes', 0),
    #    Constraint("preprocess", MUST_BE, 'tpm', 0),
    #    Constraint("preprocess", MUST_BE, 'agilent', 1),
    #    Constraint("preprocess", MUST_BE, 'affymetrix', 2),
    #    help='compare three SignatureScore'),

    ModuleNode(
        'plot_signature_predictions_comparison',
        ScoreCompareReportFile, ScoreComparePlot,
        help='plot the SignatureScore comparison',
        ),
    ]
