from Betsy.bie3 import *
import BasicDataTypes as BDT
import SignalFile
import GOAnalysis
import GSEAAnalysis
import ExpressionVisualization as EV


DE_ALGORITHMS = [
    "fold_change",
    "ttest",
    "sam",
    "ebayes",
    "deseq2",
    "edger",
    ]

DiffExpAnalysis = DataType(
    "DiffExpAnalysis",
    SignalFile.ATTR_CONTENTS,
    SignalFile.ATTR_PREPROCESS,
    AttributeDef(
        "de_algorithm", DE_ALGORITHMS, "fold_change", "fold_change", 
        help="Which algorithm used to calculate differential expression."),
    help="Differential expression result file",
    )

DiffExpVennComparison = DataType(
    "DiffExpVennComparison",
    SignalFile.ATTR_CONTENTS,
    SignalFile.ATTR_PREPROCESS,
    AttributeDef(
        "de_algorithm", DE_ALGORITHMS, "fold_change", "fold_change", 
        help="Which algorithm used to calculate differential expression."),
    help="Comparison of the results from a differential expression analysis.",
    )


## DiffReportFile = DataType(
##     'DiffReportFile',
##     help="Report file for diff exp report",
##     )

all_data_types = [
    DiffExpAnalysis,
    DiffExpVennComparison,
    #DiffReportFile,
    ]

all_modules = [
##     ModuleNode(
##         'calc_diffexp_with_ttest',
##         [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExpAnalysis,
##         ## OptionDef(
##         ##     "fold_change_cutoff", "",
##         ##     help="Keep only genes that change by at least this fold change.")
##         ## OptionDef(
##         ##     "p_cutoff", "",
##         ##     help="Keep only genes with p-values at least this significant.")
##         ## OptionDef(
##         ##     "bonf_cutoff", "",
##         ##     help="Keep only genes with adjusted p-values this significant.")
##         ## OptionDef(
##         ##     "fdr_cutoff", "",
##         ##     help="Keep only genes with FDR at least this significant.")

##         #Constraint("gene_order", MUST_BE, 'none', 0),
        
##         #Constraint("contents", MUST_BE, "class0,class1", 0),
##         Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
##         Constraint("contents", SAME_AS, 0, 1),
##         Consequence('contents', SAME_AS_CONSTRAINT, 0),
##         Constraint("logged", MUST_BE, "yes", 0),
##         Constraint("format", MUST_BE, "tdf", 0),
##         Consequence("de_algorithm", SET_TO, 'ttest'),
##         help="calculate the differential expression with ttest method"),
    
##     ModuleNode(
##         'calc_diffexp_with_sam',
##         [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExpAnalysis,
##         ## OptionDef(
##         ##     "sam_delta_value", 1.0,
##         ##     help="delta value for sam differential expression method"),
##         ## OptionDef(
##         ##     "diffexp_foldchange_value", "",
##         ##     help="fold change value for differential expression analysis"),

##         #Constraint("gene_order", MUST_BE, 'no', 0),
##         Consequence("de_algorithm", SET_TO, 'sam'),
        
##         Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
##         Constraint("contents", SAME_AS, 0, 1),
##         Consequence('contents', SAME_AS_CONSTRAINT, 0),
##         Constraint("logged", MUST_BE, 'yes', 0),
##         Constraint("format", MUST_BE, 'tdf', 0),
##         #Consequence('contents', SAME_AS_CONSTRAINT, 0),
##         help="calculate the differential expression with sam method"),
    
##     ModuleNode(
##         'calc_diffexp_with_ebayes',
##         [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExpAnalysis,
##         ## OptionDef(
##         ##      "diffexp_foldchange_value", "",
##         ##      help="fold change value for differential expression analysis"),

##         #Constraint("gene_order", MUST_BE, 'no', 0),
##         Consequence("de_algorithm", SET_TO, 'ebayes'),
        
##         Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
##         #Constraint("contents", MUST_BE, "class0,class1", 0),
##         Constraint("logged", MUST_BE, 'yes', 0),
##         Constraint("format", MUST_BE, 'tdf', 0),
##         Constraint("contents", SAME_AS, 0, 1),
##         #Consequence('contents', SAME_AS_CONSTRAINT, 0),
##         help="calculate the differential expression with ebayes method"),
    
    ModuleNode(
        "calc_diffexp_microarray",
        [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExpAnalysis,
        OptionDef(
            "fold_change_cutoff", default="",
            help="Find genes with at least this much fold change.  "
            "Leave blank for no cutoff."),
        OptionDef(
            "fdr_cutoff", default="",
            help="Use this FDR cutoff.  Leave blank for no cutoff."),
        OptionDef(
            "p_cutoff", default="",
            help="Use this p-value cutoff.  Leave blank for no cutoff."),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint(
            "preprocess", CAN_BE_ANY_OF,
            ["unknown", "mas5", "rma", "agilent", "illumina", "tpm", "fpkm"]),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Constraint("logged", MUST_BE, "yes", 0),
        Constraint("format", MUST_BE, "tdf", 0),
        Consequence(
            "de_algorithm", SET_TO_ONE_OF,
            ["fold_change", "ttest", "sam", "ebayes"]),
        help="Calculate differential expression of genes."
        ),

    ModuleNode(
        "calc_diffexp_readcounts",
        [SignalFile.SignalFile, SignalFile.ClassLabelFile], DiffExpAnalysis,
        OptionDef(
            "fold_change_cutoff", default="",
            help="Find genes with at least this much fold change.  "
            "Leave blank for no cutoff."),
        OptionDef(
            "fdr_cutoff", default="",
            help="Use this FDR cutoff.  Leave blank for no cutoff."),
        OptionDef(
            "p_cutoff", default="",
            help="Use this p-value cutoff.  Leave blank for no cutoff."),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("preprocess", MUST_BE, "counts", 0),
        Constraint("preprocess", SAME_AS, 0, 1),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Constraint("logged", MUST_BE, "no", 0),
        Constraint("format", MUST_BE, "tdf", 0),
        Consequence("de_algorithm", SET_TO_ONE_OF, ["deseq2", "edger"]),
        help="Calculate differential expression of genes from "
        "read counts (RNA-Seq data)."
        ),

    ModuleNode(
        "compare_diffexp_venn",
        DiffExpAnalysis, DiffExpVennComparison,
        OptionDef(
            "merge_up_and_down_genes", default="no",
            help="Whether to count genes that go up and down in the "
            "same group."),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("de_algorithm", CAN_BE_ANY_OF, DE_ALGORITHMS),
        Consequence("de_algorithm", SAME_AS_CONSTRAINT),
        help="Compare the genes in a differential expression analysis "
        "in a Venn diagram.",
        ),
        
    

    
    ModuleNode(
        "filter_genes_by_diffexp",
        [SignalFile._SignalFile_Filter, DiffExpAnalysis], BDT.GeneListFile,
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
    
##     ModuleNode(
##         'make_diffgenes_report',
##         [DiffExpAnalysis, DiffExpAnalysis, EV.Heatmap, GOAnalysis.GatherFile,
##          GSEAAnalysis.GseaFile], DiffReportFile,
##         OptionDef("hm_width", 20),
##         OptionDef("hm_height", 1),
        
##         Constraint("de_algorithm", MUST_BE, 'ttest', 0),
##         Constraint("de_algorithm", MUST_BE, 'sam', 1),
##         Constraint("cluster_alg", MUST_BE, 'none', 2),
##         #Constraint("de_algorithm", SAME_AS, 0, 3),
##         #Constraint("contents", MUST_BE, 'unspecified', 0),
##         #Constraint("contents", MUST_BE, 'unspecified', 1),
##         ),
    ]
