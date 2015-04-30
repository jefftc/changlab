#DiffExp
from Betsy.bie3 import *
import Database
import BasicDataTypes
import GeneExpProcessing
import GOAnalysis
import Heatmap
import GSEAAnalysis

DiffExprFile = DataType('DiffExprFile', AttributeDef(
    "gene_order", ['diff_ttest', 'diff_sam', 'diff_ebayes',
                   'diff_fold_change'], "diff_ttest", 'diff_ttest',
    help="differential method"),
                        AttributeDef("contents", Database.CONTENTS,
                                     'diff_unspecified', 'diff_unspecified',
                                     help='contents'),
                        help="Differential expression result file")
DiffReportFile = DataType('DiffReportFile',
                          help="Report file for diff exp report")
list_files = [DiffExprFile, DiffReportFile]

all_modules = [
    Module(
        'calc_diffexp_with_ttest', [
            GeneExpProcessing.ClassLabelFile, GeneExpProcessing.SignalFile
        ], DiffExprFile, OptionDef(
            "diffexp_foldchange_value", 0,
            help="fold change value for differential expression analysis"),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS, 0),
        Constraint("logged", MUST_BE, 'yes', 1),
        Constraint("format", MUST_BE, 'tdf', 1),
        Constraint("gene_order", MUST_BE, 'no', 1),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("gene_order", SET_TO, 'diff_ttest'),
        Consequence('contents', SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess",CAN_BE_ANY_OF, GeneExpProcessing.PREPROCESS),
        #Constraint("preprocess", SAME_AS,0,1),
        help="calculate the differential expression with ttest method"),
    Module(
        'calc_diffexp_with_sam', [GeneExpProcessing.ClassLabelFile,
                                  GeneExpProcessing.SignalFile], DiffExprFile,
        OptionDef("sam_delta_value", 1.0,
                  help="delta value for sam differential expression method"),
        OptionDef(
            "diffexp_foldchange_value", 0,
            help="fold change value for differential expression analysis"),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS, 0),
        Constraint("logged", MUST_BE, 'yes', 1),
        Constraint("format", MUST_BE, 'tdf', 1),
        Constraint("gene_order", MUST_BE, 'no', 1),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("gene_order", SET_TO, 'diff_sam'),
        Consequence('contents', SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess",CAN_BE_ANY_OF, GeneExpProcessing.PREPROCESS),
        #Constraint("preprocess", SAME_AS,0,1),
        help="calculate the differential expression with sam method"),
    Module(
        'calc_diffexp_with_ebayes',
        [GeneExpProcessing.ClassLabelFile,
         GeneExpProcessing.SignalFile], DiffExprFile, OptionDef(
             "diffexp_foldchange_value", 0,
             help="fold change value for differential expression analysis"),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS, 0),
        Constraint("logged", MUST_BE, 'yes', 1),
        Constraint("format", MUST_BE, 'tdf', 1),
        Constraint("gene_order", MUST_BE, 'no', 1),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("gene_order", SET_TO, 'diff_ebayes'),
        Consequence('contents', SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess",CAN_BE_ANY_OF, GeneExpProcessing.PREPROCESS),
        #Constraint("preprocess", SAME_AS,0,1),
        help="calculate the differential expression with ebayes method"),
    Module(
        'calc_diffexp_with_fold_change',
        [GeneExpProcessing.ClassLabelFile,
         GeneExpProcessing.SignalFile], DiffExprFile, OptionDef(
             "diffexp_foldchange_value", 0,
             help="fold change value for differential expression analysis"),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS, 0),
        Constraint("logged", MUST_BE, 'yes', 1),
        Constraint("format", MUST_BE, 'tdf', 1),
        Constraint("gene_order", MUST_BE, 'no', 1),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("gene_order", SET_TO, 'diff_fold_change'),
        Consequence('contents', SAME_AS_CONSTRAINT, 0),
        #Constraint("preprocess",CAN_BE_ANY_OF, GeneExpProcessing.PREPROCESS),
        #Constraint("preprocess", SAME_AS,0,1),
        help="calculate the differential expression with fold change method"),
    Module(
        'generate_genelist_from_diffexprfile', DiffExprFile,
        BasicDataTypes.GeneListFile,
        OptionDef("select_gene_by_foldchange", "",
                  help="select gene by foldchange value,etc.5"),
        OptionDef("select_gene_by_p_value", "",
                  help="select gene by p-value value, etc. 0.05"), OptionDef(
                      "select_gene_by_fdr", "",
                      help="select gene by fdr value, etc.0.05"), Constraint(
                          "gene_order", CAN_BE_ANY_OF,
                          ['diff_ttest', 'diff_sam', 'diff_ebayes',
                           'diff_fold_change']), Constraint("contents",
                                                            CAN_BE_ANY_OF,
                                                            Database.CONTENTS),
        Consequence("gene_order", SAME_AS_CONSTRAINT),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help=
        "generate genelist from the result file of differential expression analysis"),
    Module('make_diffgenes_report',
           [DiffExprFile, DiffExprFile, Heatmap.Heatmap, GOAnalysis.GatherFile,
            GSEAAnalysis.GseaFile], DiffReportFile, OptionDef("hm_width", 20),
           OptionDef("hm_height", 1),
           Constraint("gene_order", MUST_BE, 'diff_ttest', 0),
           Constraint("gene_order", SAME_AS, 0, 3),
           Constraint("contents", MUST_BE, 'unspecified', 0),
           Constraint("contents", MUST_BE, 'unspecified', 1),
           Constraint("gene_order", MUST_BE, 'diff_sam', 1),
           Constraint("cluster_alg", MUST_BE, 'no_cluster_alg', 2)),
]
