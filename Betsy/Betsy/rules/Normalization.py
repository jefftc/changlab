from Betsy.bie3 import *
import PcaAnalysis
import GeneExpProcessing
import BasicDataTypes as BDT
import ArrayPlatforms

YESNO = BDT.YESNO


PreprocessingReport = DataType(
    "PreprocessingReport",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, "unknown", "any",
        help="preprocess for normalize report file"),
    AttributeDef(
        "quantile_norm", YESNO, "no", "no", help="quantile normalization"),
    AttributeDef(
        "combat_norm", YESNO, "no", "no", help="combat normalization"),
    AttributeDef(
        "dwd_norm", YESNO, "no", "no", help="dwd normalization"),
    AttributeDef(
        "shiftscale_norm", YESNO, "no", "no", help="shiftscale normalization"),
    AttributeDef(
        "bfrm_norm", YESNO, "no", "no", help="bfrm normalization"),
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "unknown", "no", help="gene center method"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "unknown", "no", help="gene normalize method"),
    AttributeDef(
        # Why is the u133A?
        "platform", ["yes", "no", 'u133A'], "no", "no",
        help="add platform or not"),
    AttributeDef(
        "annotate", YESNO, "no", "yes", help="annotate file or not"),
    AttributeDef(
        "num_features", ["yes", "no"], "no", "no",
        help="select a num of features or not"),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no", help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "closest_probe", "high_var_probe"],
        "no", "no", help="method to remove duplicated probes"),
    AttributeDef(
        "group_fc", ["yes", "no"], "no", "no",
        help="group fold change or not.  WHAT IS THIS?"),
    help="Preprocess gene expression data with some quality checks.")


BatchEffectReport = DataType(
    'BatchEffectReport',
    help="Report file for batch effect remove report"
    )

all_data_types = [
    PreprocessingReport,
    BatchEffectReport,
    ]


# Illumina has a different preprocessing report.
PREPROCESS_NOT_ILLUMINA = [
    x for x in BDT.PREPROCESS if x != "illumina"]

all_modules = [
    ModuleNode(
        'make_preprocessing_report',
        [
            GeneExpProcessing.SignalFile,
            GeneExpProcessing.IntensityPlot,
            ArrayPlatforms.ControlPlot,
            PcaAnalysis.PcaPlot,
            ArrayPlatforms.ActbPlot,
            PcaAnalysis.PcaPlot,
            ],
        PreprocessingReport,

        #SignalFile
        Constraint(
            'preprocess', CAN_BE_ANY_OF, PREPROCESS_NOT_ILLUMINA, 0),
        Constraint("annotate", MUST_BE, "yes", 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint('quantile_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('combat_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('shiftscale_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('bfrm_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('dwd_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('gene_center', CAN_BE_ANY_OF, ['median', 'mean', 'no'], 0),
        Constraint(
            'gene_normalize', CAN_BE_ANY_OF,
            ['variance', 'sum_of_squares', 'no'], 0),
        Constraint(
            'unique_genes', CAN_BE_ANY_OF,
            ['no', 'average_genes', 'high_var', 'first_gene'], 0),
        Constraint('platform', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('group_fc', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('num_features', CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            'duplicate_probe', CAN_BE_ANY_OF,
            ["no", "closest_probe", "high_var_probe"], 0),
        
        # First PcaPlot.
        Constraint('quantile_norm', SAME_AS, 0, 3),
        Constraint('combat_norm', SAME_AS, 0, 3),
        Constraint('shiftscale_norm', SAME_AS, 0, 3),
        Constraint('bfrm_norm', SAME_AS, 0, 3),
        Constraint('dwd_norm', SAME_AS, 0, 3),
        Constraint('gene_center', SAME_AS, 0, 3),
        Constraint('gene_normalize', SAME_AS, 0, 3),
        Constraint('unique_genes', SAME_AS, 0, 3),
        Constraint('platform', SAME_AS, 0, 3),
        Constraint('group_fc', SAME_AS, 0, 3),
        Constraint('num_features', SAME_AS, 0, 3),
        Constraint('duplicate_probe', SAME_AS, 0, 3),
        
        # Second PcaPlot.
        Constraint('quantile_norm', MUST_BE, 'no', 5),
        Constraint('combat_norm', MUST_BE, 'no', 5),
        Constraint('shiftscale_norm', MUST_BE, 'no', 5),
        Constraint('bfrm_norm', MUST_BE, 'no', 5),
        Constraint('dwd_norm', MUST_BE, 'no', 5),
        Constraint('gene_center', MUST_BE, 'no', 5),
        Constraint('gene_normalize', MUST_BE, 'no', 5),
        Constraint('unique_genes', MUST_BE, 'no', 5),
        Constraint('platform', MUST_BE, 'no', 5),
        Constraint('group_fc', MUST_BE, 'no',5 ),
        Constraint('num_features', MUST_BE, 'no', 5),
        Constraint('duplicate_probe', MUST_BE, 'no', 5),
        
        Constraint('contents', SAME_AS, 0, 1),
        Constraint('contents', SAME_AS, 0, 2),
        Constraint('contents', SAME_AS, 0, 3),
        Constraint('contents', SAME_AS, 0, 4),
        Constraint('contents', SAME_AS, 0, 5),
        Constraint("preprocess", SAME_AS, 0, 1),
        Constraint("preprocess", SAME_AS, 0, 2),
        Constraint("preprocess", SAME_AS, 0, 3),
        Constraint("preprocess", SAME_AS, 0, 4),
        Constraint("preprocess", SAME_AS, 0, 5),
        Consequence('preprocess', SAME_AS_CONSTRAINT),
        help="make preprocessing report"),
    
    ModuleNode(
        'make_preprocessing_report_illumina',
        [
            GeneExpProcessing.SignalFile,     # 0
            GeneExpProcessing.IntensityPlot,  # 1
            ArrayPlatforms.BiotinPlot,        # 2
            PcaAnalysis.PcaPlot,              # 3
            ArrayPlatforms.ActbPlot,          # 4
            PcaAnalysis.PcaPlot,              # 5
            ArrayPlatforms.HousekeepingPlot,  # 6
            ArrayPlatforms.Hyb_barPlot,       # 7
            ArrayPlatforms.ControlFile        # 8
            ],
        PreprocessingReport,
        
        # SignalFile
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("preprocess", MUST_BE, "illumina", 0),
        Constraint('quantile_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('combat_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('dwd_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('shiftscale_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('bfrm_norm', CAN_BE_ANY_OF, YESNO, 0),
        Constraint('gene_center', CAN_BE_ANY_OF, ['median', 'mean', 'no'], 0),
        Constraint(
            'gene_normalize', CAN_BE_ANY_OF,
            ['variance', 'sum_of_squares', 'no'], 0),
        Constraint("platform", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("annotate", MUST_BE, "yes", 0),
        Constraint('num_features', CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            'unique_genes', CAN_BE_ANY_OF,
            ['no', 'average_genes', 'high_var', 'first_gene'], 0),
        Constraint(
            'duplicate_probe', CAN_BE_ANY_OF,
            ["no", "closest_probe", "high_var_probe"], 0),
        Constraint('group_fc', CAN_BE_ANY_OF, YESNO, 0),

        # IntensityPlot
        Constraint("contents", SAME_AS, 0, 1),
        Constraint("preprocess", SAME_AS, 0, 1),

        # BiotinPlot
        Constraint("contents", SAME_AS, 0, 2),
        
        # PcaPlot (preprocessed).
        Constraint("contents", SAME_AS, 0, 3),
        Constraint("preprocess", SAME_AS, 0, 3),
        Constraint('quantile_norm', SAME_AS, 0, 3),
        Constraint('combat_norm', SAME_AS, 0, 3),
        Constraint('dwd_norm', SAME_AS, 0, 3),
        Constraint('shiftscale_norm', SAME_AS, 0, 3),
        Constraint('bfrm_norm', SAME_AS, 0, 3),
        Constraint('gene_center', SAME_AS, 0, 3),
        Constraint('gene_normalize', SAME_AS, 0, 3),
        Constraint('unique_genes', SAME_AS, 0, 3),
        Constraint('platform', SAME_AS, 0, 3),
        Constraint('num_features', SAME_AS, 0, 3),
        Constraint('duplicate_probe', SAME_AS, 0, 3),
        Constraint('group_fc', SAME_AS, 0, 3),

        # ActbPlot
        Constraint("contents", SAME_AS, 0, 4),
        Constraint("preprocess", SAME_AS, 0, 4),
        
        # PcaPlot (No processing)
        Constraint("contents", SAME_AS, 0, 5),
        Constraint("preprocess", SAME_AS, 0, 5),
        Constraint('quantile_norm', MUST_BE, "no", 5),
        Constraint('combat_norm', MUST_BE, "no", 5),
        Constraint('dwd_norm', MUST_BE, "no", 5),
        Constraint('shiftscale_norm', MUST_BE, "no", 5),
        Constraint('bfrm_norm', MUST_BE, "no", 5),
        Constraint('gene_center', MUST_BE, "no", 5),
        Constraint('gene_normalize', MUST_BE, "no", 5),
        Constraint('unique_genes', MUST_BE, "no", 5),
        Constraint('platform', MUST_BE, "no", 5),
        Constraint('num_features', MUST_BE, "no", 5),
        Constraint('duplicate_probe', MUST_BE, "no", 5),
        Constraint('group_fc', MUST_BE, "no", 5),

        # HousekeepingPlot
        Constraint("contents", SAME_AS, 0, 6),

        # Hyb_barPlot
        Constraint("contents", SAME_AS, 0, 7),

        # ControlFile
        Constraint("preprocess", SAME_AS, 0, 8),
        Constraint("contents", SAME_AS, 0, 8),
        Constraint('format', MUST_BE, "gct", 8),
        Constraint("logged", MUST_BE, "no", 8),
        
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("gene_center", SAME_AS_CONSTRAINT, 0),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT, 0),
        Consequence("platform", SAME_AS_CONSTRAINT, 0),
        Consequence("annotate", SAME_AS_CONSTRAINT, 0),
        Consequence("num_features", SAME_AS_CONSTRAINT, 0),
        Consequence("unique_genes", SAME_AS_CONSTRAINT, 0),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT, 0),
        Consequence("group_fc", SAME_AS_CONSTRAINT, 0),
        
        help="make preprocessing report for illumina",
        ),
    
    ModuleNode(
        "make_batch_effect_report",
        [
            GeneExpProcessing.SignalFile,
            PcaAnalysis.PcaPlot,
         
            GeneExpProcessing.SignalFile,
            PcaAnalysis.PcaPlot,
         
            GeneExpProcessing.SignalFile,
            PcaAnalysis.PcaPlot,
         
            GeneExpProcessing.SignalFile,
            PcaAnalysis.PcaPlot,
         
            GeneExpProcessing.SignalFile,
            PcaAnalysis.PcaPlot,
         
            GeneExpProcessing.SignalFile,
            PcaAnalysis.PcaPlot
            ],
        BatchEffectReport,
        
        Constraint("quantile_norm", MUST_BE, "no", 0),
        Constraint("quantile_norm", SAME_AS, 0, 1),
        Constraint("quantile_norm", MUST_BE, "yes", 2),
        Constraint("quantile_norm", SAME_AS, 2, 3),
        
        Constraint("quantile_norm", MUST_BE, "yes", 4),
        Constraint("dwd_norm", MUST_BE, "yes",4),
        Constraint("quantile_norm", SAME_AS, 4, 5),
        Constraint("dwd_norm", SAME_AS, 4, 5),
        
        Constraint("quantile_norm", MUST_BE, "yes",6),
        Constraint("bfrm_norm", MUST_BE, "yes",6),
        Constraint("quantile_norm", SAME_AS,6,7),
        Constraint("bfrm_norm", SAME_AS,6,7),
        
        Constraint("quantile_norm", MUST_BE, "yes",8),
        Constraint("shiftscale_norm", MUST_BE, "yes",8),
        Constraint("quantile_norm", SAME_AS, 8, 9),
        Constraint("shiftscale_norm", SAME_AS, 8, 9),
        
        Constraint("quantile_norm", MUST_BE, "yes", 10),
        Constraint("combat_norm", MUST_BE, "yes", 10),
        Constraint("quantile_norm", SAME_AS, 10, 11),
        Constraint("combat_norm", SAME_AS, 10, 11),
        
        help="make batch effect remove report",
        ),
    
    ModuleNode(
        'convert_preprocessing_report_preprocess',
        PreprocessingReport, PreprocessingReport,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SET_TO, "any"),
        help="convert preprocess from others to any in PreprocessingReport",
        ),
    ]
