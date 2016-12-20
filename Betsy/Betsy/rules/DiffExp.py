from Betsy.bie3 import *
import BasicDataTypes as BDT
import SignalFile
import GOAnalysis
import GSEAAnalysis
#import Clustering
import ExpressionVisualization as EV


DE_ALGORITHM = ["fold_change", "ttest", "sam", "ebayes"]

DiffExprFile = DataType(
    'DiffExprFile',
    AttributeDef(
        "de_algorithm", DE_ALGORITHM, "fold_change", "fold_change", 
        help="Which algorithm used to calculate differential expression."),
    #AttributeDef(
    #    "contents", BDT.CONTENTS, 'diff_unspecified', 'diff_unspecified',
    #    help='contents'),
    help="Differential expression result file")

DiffReportFile = DataType(
    'DiffReportFile',
    help="Report file for diff exp report",
    )

all_data_types = [
    DiffExprFile,
    DiffReportFile,
    ]

all_modules = [
    ModuleNode(
        'calc_diffexp_with_ttest',
        [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExprFile,
        ## OptionDef(
        ##     "fold_change_cutoff", "",
        ##     help="Keep only genes that change by at least this fold change.")
        ## OptionDef(
        ##     "p_cutoff", "",
        ##     help="Keep only genes with p-values at least this significant.")
        ## OptionDef(
        ##     "bonf_cutoff", "",
        ##     help="Keep only genes with adjusted p-values this significant.")
        ## OptionDef(
        ##     "fdr_cutoff", "",
        ##     help="Keep only genes with FDR at least this significant.")

        #Constraint("gene_order", MUST_BE, 'none', 0),
        Consequence("de_algorithm", SET_TO, 'ttest'),
        
        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("logged", MUST_BE, 'yes', 0),
        Constraint("format", MUST_BE, 'tdf', 0),
        Constraint("contents", SAME_AS, 0, 1),
        #Consequence('contents', SAME_AS_CONSTRAINT, 0),
        help="calculate the differential expression with ttest method"),
    
    ModuleNode(
        'calc_diffexp_with_sam',
        [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExprFile,
        ## OptionDef(
        ##     "sam_delta_value", 1.0,
        ##     help="delta value for sam differential expression method"),
        ## OptionDef(
        ##     "diffexp_foldchange_value", "",
        ##     help="fold change value for differential expression analysis"),

        #Constraint("gene_order", MUST_BE, 'no', 0),
        Consequence("de_algorithm", SET_TO, 'sam'),
        
        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("logged", MUST_BE, 'yes', 0),
        Constraint("format", MUST_BE, 'tdf', 0),
        Constraint("contents", SAME_AS, 0, 1),
        #Consequence('contents', SAME_AS_CONSTRAINT, 0),
        help="calculate the differential expression with sam method"),
    
    ModuleNode(
        'calc_diffexp_with_ebayes',
        [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExprFile,
        ## OptionDef(
        ##      "diffexp_foldchange_value", "",
        ##      help="fold change value for differential expression analysis"),

        #Constraint("gene_order", MUST_BE, 'no', 0),
        Consequence("de_algorithm", SET_TO, 'ebayes'),
        
        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("logged", MUST_BE, 'yes', 0),
        Constraint("format", MUST_BE, 'tdf', 0),
        Constraint("contents", SAME_AS, 0, 1),
        #Consequence('contents', SAME_AS_CONSTRAINT, 0),
        help="calculate the differential expression with ebayes method"),
    
    ModuleNode(
        'calc_diffexp_with_fold_change',
        [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExprFile,
        ## OptionDef(
        ##     "fold_change_cutoff", "",
        ##     help="Keep only genes that change by at least this fold change.")

        #Constraint("gene_order", MUST_BE, 'no', 0),
        Consequence("de_algorithm", SET_TO, 'fold_change'),
        
        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("logged", MUST_BE, 'yes', 0),
        Constraint("format", MUST_BE, 'tdf', 0),
        Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        help="calculate the differential expression with fold change method"),
    
    ModuleNode(
        "filter_genes_by_diffexp",
        [SignalFile._SignalFile_Filter, DiffExprFile], BDT.GeneListFile,
        OptionDef(
            "fold_change_cutoff", "2",
            help="Filter genes with less than this fold change."),
        OptionDef("p_value_cutoff", "0.05", help="Filter genes by p-value."),
        OptionDef("select_gene_by_fdr", "", help="Filter genes by FDR."),

        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS, 0),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("filter_missing_values", MUST_BE, "no", 0),
                   
        #Constraint("algorithm", CAN_BE_ANY_OF, DE_ALGORITHM),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("algorithm", SAME_AS_CONSTRAINT),
        Consequence("filtered", SET_TO, "yes"),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        help="Filter the genes from a SignalFile based on the results "
        "of a differential expression analysis.",
        ),
    
    ModuleNode(
        'make_diffgenes_report',
        [DiffExprFile, DiffExprFile, EV.Heatmap, GOAnalysis.GatherFile,
         GSEAAnalysis.GseaFile], DiffReportFile,
        OptionDef("hm_width", 20),
        OptionDef("hm_height", 1),
        
        Constraint("de_algorithm", MUST_BE, 'ttest', 0),
        Constraint("de_algorithm", MUST_BE, 'sam', 1),
        Constraint("cluster_alg", MUST_BE, 'none', 2),
        #Constraint("de_algorithm", SAME_AS, 0, 3),
        #Constraint("contents", MUST_BE, 'unspecified', 0),
        #Constraint("contents", MUST_BE, 'unspecified', 1),
        ),
    ]
