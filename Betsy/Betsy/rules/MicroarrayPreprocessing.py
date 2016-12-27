# DataTypes:
# CELFiles
# AgilentFiles
# GPRFiles
# IDATFiles
# ILLUFolder
#
# CompleteExpressionPreprocessing
# ExpressionPreprocessingReport
# BatchEffectReport
# 
# Modules:
# extract_matrix_data
# 
# extract_CEL_files
# detect_CEL_version
# convert_CEL_to_v3
# preprocess_rma
# preprocess_mas5
#
# extract_illumina_idat_files
# preprocess_illumina
# get_illumina_control
# get_illumina_signal
# 
# extract_agilent_files
# preprocess_agilent
# extract_gpr_files
# normalize_with_loess     # move somewhere else?
#
# Variables:
# ILLU_MANIFEST
# ILLU_CHIP
# DEFAULT_MANIFEST
# DEFAULT_CHIP


from Betsy.bie3 import *
import SignalFile
import BasicDataTypes as BDT
import ExpressionVisualization as EV
YESNO = BDT.YESNO

ILLU_MANIFEST = [
    'HumanHT-12_V4_0_R1_15002873_B.txt',
    'HumanHT-12_V4_0_R2_15002873_B.txt',
    'HumanHT-12_V3_0_R2_11283641_A.txt',
    'HumanHT-12_V3_0_R3_11283641_A.txt',
    'HumanMI_V1_R2_XS0000122-MAP.txt',
    'HumanMI_V2_R0_XS0000124-MAP.txt',
    'HumanRef-8_V2_0_R4_11223162_A.txt',
    'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
    'HumanRef-8_V3_0_R2_11282963_A.txt',
    'HumanRef-8_V3_0_R3_11282963_A.txt',
    'HumanWG-6_V2_0_R4_11223189_A.txt',
    'HumanWG-6_V3_0_R2_11282955_A.txt',
    'HumanWG-6_V3_0_R3_11282955_A.txt',
    'MouseMI_V1_R2_XS0000127-MAP.txt',
    'MouseMI_V2_R0_XS0000129-MAP.txt',
    'MouseRef-8_V1_1_R4_11234312_A.txt',
    'MouseRef-8_V2_0_R2_11278551_A.txt',
    'MouseRef-8_V2_0_R3_11278551_A.txt',
    'MouseWG-6_V1_1_R4_11234304_A.txt',
    'MouseWG-6_V2_0_R2_11278593_A.txt',
    'MouseWG-6_V2_0_R3_11278593_A.txt',
    'RatRef-12_V1_0_R5_11222119_A.txt',
    ]

ILLU_CHIP = [
    'ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
    'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
    'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
    'ilmn_HumanRef_8_V3_0_R1_11282963_A_WGDASL.chip',
    'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
    'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
    'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
    'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
    'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
    'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
    'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
    'ilmn_RatRef_12_V1_0_R5_11222119_A.chip',
    ]
DEFAULT_MANIFEST = "HumanHT-12_V4_0_R2_15002873_B.txt"
DEFAULT_CHIP = "ilmn_HumanHT_12_V4_0_R1_15002873_B.chip"

# Illumina has a different preprocessing report.
PREPROCESS_NOT_ILLUMINA = [
    x for x in BDT.PREPROCESS if x != "illumina"]
del x  # Prevent warning: Local variable (x) shadows global defined on ...

CELFiles = DataType(
    "CELFiles",
    SignalFile.ATTR_CONTENTS,
    AttributeDef(
        "version", ["unknown", "cc", "v3_v4"], "unknown", "v3_v4",
        help="cel file version"),
    help="A folder of cel files")

AgilentFiles = DataType(
    "AgilentFiles",
    SignalFile.ATTR_CONTENTS,
    help="A folder of agilent files.")

GPRFiles = DataType(
    "GPRFiles",
    SignalFile.ATTR_CONTENTS,
    help="A folder of GPR files.")

IDATFiles = DataType(
    "IDATFiles",
    SignalFile.ATTR_CONTENTS,
    help="A folder of IDAFiles.")

ILLUFolder = DataType(
    "ILLUFolder",
    #AttributeDef(
    #    "illu_manifest", ILLU_MANIFEST, 'HumanHT-12_V4_0_R2_15002873_B.txt',
    #    'HumanHT-12_V4_0_R2_15002873_B.txt', help="illumina manifest"),
    #AttributeDef(
    #    'illu_chip', ILLU_CHIP, 'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
    #    'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
    #    help="illumina chip type"),
    #AttributeDef(
    #    'illu_bg_mode', ['false', 'true'], "false", "false",
    #    help="illumina background mode"),
    #AttributeDef(
    #    'illu_coll_mode', ['none', 'max', 'median'], "none", "none",
    #    help="illumina coll mode"),
    SignalFile.ATTR_CONTENTS,
    help="A folder generated from preprocess_illumina," \
    "it contains UnprocessedSignalFile and IlluminaControlFile.")


CompleteExpressionPreprocessing = DataType(
    "CompleteExpressionPreprocessing",
    AttributeDef("logged", YESNO, "yes", "yes", help="logged or not"),
    # Handle UNPROC_ATTRIBUTES myself, since "preprocess" can also
    # take value of "any".
    SignalFile.ATTR_CONTENTS,
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, "unknown", "any",
        help="preprocess method"),
    # Handle ANNOTATE_ATTRIBUTES myself, to set default "annotate" to
    # "yes".
    AttributeDef(
        "annotate", YESNO, "no", "yes",
        help="annotate file or not [WHAT DOES THIS MEAN?]"),
    AttributeDef(
        "relabel_sample", YESNO, "no", "no",
        help="rename sample or not"),
    AttributeDef(
        # Why is the u133A?
        "platform", ["yes", "no", 'u133A'], "no", "no",
        help="add platform or not"),
    
    *(SignalFile.IMPUTE_ATTRIBUTES +
      SignalFile.MERGE_ATTRIBUTES +
      SignalFile.NORMALIZE_ATTRIBUTES +
      SignalFile.ORDER_ATTRIBUTES +
      #SignalFile.ANNOTATE_ATTRIBUTES +
      SignalFile.FILTER_ATTRIBUTES),
    help="Preprocess gene expression data with some quality checks."
    )


ExpressionPreprocessingReport = DataType(
    'ExpressionPreprocessingReport',
    AttributeDef(
        "preprocess", BDT.ANY_PREPROCESS, "unknown", "any",
        help="preprocess method"),
    help="Summarize the expression preprocessing in a report.  NOT IMPLEMENTED"
    )


BatchEffectReport = DataType(
    'BatchEffectReport',
    help="Report file for batch effect remove report.  NOT IMPLEMENTED."
    )


all_data_types = [
    CELFiles,
    AgilentFiles,
    GPRFiles,
    IDATFiles,
    ILLUFolder,

    CompleteExpressionPreprocessing,
    ExpressionPreprocessingReport,
    BatchEffectReport,
    ]

all_modules = [
    ModuleNode(
        "extract_matrix_data",
        BDT.ExpressionFiles, SignalFile.UnprocessedSignalFile,
        Constraint("filetype", MUST_BE, 'matrix'),
        Consequence("preprocess", SET_TO, "unknown"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "unknown"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        help="extract SignalFile_Postprocess files from Expression Files"),
    
    #CELFiles
    ModuleNode(
        "extract_CEL_files",
        BDT.ExpressionFiles, CELFiles,
        Constraint("filetype", MUST_BE, 'cel'),
        Consequence("version", SET_TO, "unknown"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract CEL files from Expression Files folder"),
    
    ModuleNode(
        "detect_CEL_version", CELFiles, CELFiles,
        Constraint("version", MUST_BE, "unknown"),
        Consequence("version", BASED_ON_DATA, ["cc", "v3_v4"]),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="detect the version of cel files"),

    ModuleNode(
        "convert_CEL_to_v3",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "cc"),
        Consequence("version", SET_TO, "v3_v4"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="convert cel files to v3 version"),
    
    ModuleNode(
        "preprocess_rma",
        CELFiles, SignalFile.UnprocessedSignalFile,
        Constraint("version", MUST_BE, 'v3_v4'),
        Consequence("logged", SET_TO, "yes"),
        Consequence("preprocess", SET_TO, "rma"),
        Consequence("format", SET_TO, "jeffs"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help=
        "preprocess CELFiles with rma method,generate SignalFile_Postprocess",
        ),
    
    ModuleNode(
        "preprocess_mas5",
        CELFiles, SignalFile.UnprocessedSignalFile,
        Constraint("version", MUST_BE, 'v3_v4'),
        Consequence("logged", SET_TO, "no"),
        Consequence("preprocess", SET_TO, "mas5"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Consequence("format", SET_TO, "jeffs"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="preprocess CELFiles with mas5 method,generate "
        "SignalFile_Postprocess",
        ),
    
    # IDATFiles
    ModuleNode(
        "extract_illumina_idat_files",
        BDT.ExpressionFiles, IDATFiles,
        Constraint("filetype", MUST_BE, 'idat'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract idat files from Expression File folder"),
    
    ModuleNode(
        "preprocess_illumina",
        IDATFiles, ILLUFolder,
        #Consequence('illu_manifest', SET_TO_ONE_OF, ILLU_MANIFEST),
        #Consequence('illu_chip', SET_TO_ONE_OF, ILLU_CHIP),
        #Consequence('illu_bg_mode', SET_TO_ONE_OF, ["false", "true"]),
        #Consequence(
        #    'illu_coll_mode', SET_TO_ONE_OF, ["none", "max", "median"]),

        OptionDef(
            "illu_manifest", DEFAULT_MANIFEST,
            help="manifest: %s" % ", ".join(ILLU_MANIFEST)),
        OptionDef(
            'illu_chip', DEFAULT_CHIP,
            help="CHIP file to map probes to genes: %s" %
            ", ".join(ILLU_CHIP)),
        OptionDef(
            "illu_bg_mode", "false", 
            help="background subtraction (false or true)"),
        OptionDef(
            'illu_coll_mode', "none",
            help="collapse probes to genes (none, max, or median)"),
        
        OptionDef(
            "illu_clm", '', help="CLM file for mapping file to sample names"),
        OptionDef("illu_custom_chip", '', help=""),
        OptionDef("illu_custom_manifest", '', help=""),

        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="Preprocess .idat files into signal values",
        ),
    
    ModuleNode(
        "get_illumina_control",
        ILLUFolder, SignalFile.IlluminaControlFile,
        #Consequence('preprocess', SET_TO, "illumina"),
        #Consequence("logged", SET_TO, "no"),
        #Consequence("format", SET_TO, "gct"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract illumina ControlFile from ILLUFolder"),
    
    ModuleNode(
        "get_illumina_signal",
        ILLUFolder, SignalFile.UnprocessedSignalFile,
        Consequence('preprocess', SET_TO, "illumina"),
        Consequence('format', SET_TO, "gct"),
        Consequence('logged', SET_TO, "no"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract the SignalFile_Postprocess from ILLUFolder"),
    
    # AgilentFiles
    ModuleNode(
        "extract_agilent_files",
        BDT.ExpressionFiles, AgilentFiles,
        Constraint("filetype", MUST_BE, 'agilent'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract agilent files from ExpressionFiles"),
    
    ModuleNode(
        "preprocess_agilent",
        AgilentFiles, SignalFile.UnprocessedSignalFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence('logged', SET_TO, "unknown"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Consequence('preprocess', SET_TO, "agilent"),
        Consequence('format', SET_TO, "tdf"),
        help="preprocess agilent, generate SignalFile_Postprocess"),
    
    # GPRFiles
    ModuleNode(
        "extract_gpr_files",
        BDT.ExpressionFiles, GPRFiles,
        Constraint("filetype", MUST_BE, 'gpr'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract gpr files from ExpressionFiles"),
    
    ModuleNode(
        "normalize_with_loess",
        GPRFiles, SignalFile.UnprocessedSignalFile,
        Consequence("format", SET_TO, "tdf"),
        Consequence("logged", SET_TO, "unknown"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Consequence("preprocess", SET_TO, "loess"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="normalize GPRFiles,generate SignalFile_Postprocess"),

    ModuleNode(
        'do_complete_preprocessing',
        [
            SignalFile.SignalFile,           #  0   no normalization
            SignalFile.SignalFile,           #  1   normalized
            EV.SignalDistributionBoxplot,    #  2
            EV.ActbPlot,                     #  3   no normalization
            EV.ActbPlot,                     #  4   normalized
            EV.PCAPlot,                      #  5   no normalization
            EV.PCAPlot,                      #  6   normalized
            EV.Heatmap,                      #  7   no normalization
            EV.Heatmap,                      #  8   normalized
         ],
        CompleteExpressionPreprocessing,

        # SignalFile (not normalized)
        Constraint("logged", CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            "preprocess", CAN_BE_ANY_OF,
            [x for x in BDT.PREPROCESS if x != "illumina"], 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"], 0),
        Constraint('quantile_norm', MUST_BE, "no", 0),
        Constraint('combat_norm', MUST_BE, "no", 0),
        Constraint('dwd_norm', MUST_BE, "no", 0),
        Constraint('shiftscale_norm', MUST_BE, "no", 0),
        Constraint('bfrm_norm', MUST_BE, "no", 0),
        Constraint('gene_center', MUST_BE, "no", 0),
        Constraint('gene_normalize', MUST_BE, "no", 0),
        Constraint("gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER, 0),
        Constraint('annotate', CAN_BE_ANY_OF, YESNO, 0),
        Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("platform", CAN_BE_ANY_OF, ["no", "yes", 'u133A'], 0),
        Constraint('num_features', CAN_BE_ANY_OF, YESNO, 0),
        Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            'unique_genes', CAN_BE_ANY_OF,
            ['no', 'average_genes', 'high_var', 'first_gene'], 0),
        Constraint(
            'duplicate_probe', CAN_BE_ANY_OF,
            ["no", "closest_probe", "high_var_probe"], 0),
        Constraint('group_fc', CAN_BE_ANY_OF, YESNO, 0),

        # SignalFile (normalized)
        Constraint('logged', SAME_AS, 0, 1),
        Constraint("contents", SAME_AS, 0, 1),
        Constraint("preprocess", SAME_AS, 0, 1),
        Constraint("filter_missing_values", SAME_AS, 0, 1),
        Constraint("missing_algorithm", SAME_AS, 0, 1),
        Constraint('quantile_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('combat_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('dwd_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('shiftscale_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('bfrm_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('gene_center', CAN_BE_ANY_OF, ['median', 'mean', 'no'], 1),
        Constraint(
            'gene_normalize', CAN_BE_ANY_OF,
            ['variance', 'sum_of_squares', 'no'], 1),
        Constraint("gene_order", SAME_AS, 0, 1),
        Constraint("annotate", SAME_AS, 0, 1),
        Constraint("relabel_sample", SAME_AS, 0, 1),
        Constraint("platform", SAME_AS, 0, 1),
        Constraint('num_features', SAME_AS, 0, 1),
        Constraint("filter_and_threshold", SAME_AS, 0, 1),
        Constraint('unique_genes', SAME_AS, 0, 1),
        Constraint('duplicate_probe', SAME_AS, 0, 1),
        Constraint('group_fc', SAME_AS, 0, 1),

        # SignalDistributionBoxplot
        Constraint("logged", MUST_BE, "yes", 2),
        Constraint("contents", SAME_AS, 0, 2),
        Constraint("preprocess", SAME_AS, 0, 2),
        Constraint('relabel_sample', SAME_AS, 0, 2),

        # ActbPlot (No normalization)
        Constraint('logged', MUST_BE, "yes", 3),
        Constraint("preprocess", SAME_AS, 0, 3),
        Constraint("contents", SAME_AS, 0, 3),
        Constraint('quantile_norm', SAME_AS, 0, 3),
        Constraint('combat_norm', SAME_AS, 0, 3),
        Constraint('dwd_norm', SAME_AS, 0, 3),
        Constraint('shiftscale_norm', SAME_AS, 0, 3),
        Constraint('bfrm_norm', SAME_AS, 0, 3),
        Constraint('gene_center', SAME_AS, 0, 3),
        Constraint('gene_normalize', SAME_AS, 0, 3),
        Constraint('annotate', MUST_BE, "yes", 3),
        Constraint('relabel_sample', SAME_AS, 0, 3),
        Constraint('platform', SAME_AS, 0, 3),
        Constraint('num_features', SAME_AS, 0, 3),
        Constraint('filter_and_threshold', MUST_BE, "no", 3),
        Constraint('unique_genes', MUST_BE, "no", 3),
        Constraint('duplicate_probe', MUST_BE, "no", 3),
        Constraint('group_fc', MUST_BE, "no", 3),

        # ActbPlot (Normalized)
        Constraint('logged', MUST_BE, "yes", 4),
        Constraint("contents", SAME_AS, 0, 4),
        Constraint("preprocess", SAME_AS, 0, 4),
        Constraint('quantile_norm', SAME_AS, 1, 4),
        Constraint('combat_norm', SAME_AS, 1, 4),
        Constraint('dwd_norm', SAME_AS, 1, 4),
        Constraint('shiftscale_norm', SAME_AS, 1, 4),
        Constraint('bfrm_norm', SAME_AS, 1, 4),
        Constraint('gene_center', SAME_AS, 1, 4),
        Constraint('gene_normalize', SAME_AS, 1, 4),
        Constraint('annotate', MUST_BE, "yes", 4),
        Constraint('relabel_sample', SAME_AS, 0, 4),
        Constraint('platform', SAME_AS, 0, 4),
        Constraint('num_features', SAME_AS, 0, 4),
        Constraint('filter_and_threshold', MUST_BE, "no", 4),
        Constraint('unique_genes', SAME_AS, 0, 4),
        Constraint('duplicate_probe', SAME_AS, 0, 4),
        Constraint('group_fc', SAME_AS, 0, 4),

        # PcaPlot (No normalization)
        Constraint("logged", MUST_BE, "yes", 5),
        Constraint("contents", SAME_AS, 0, 5),
        Constraint("preprocess", SAME_AS, 0, 5),
        Constraint('quantile_norm', SAME_AS, 0, 5),
        Constraint('combat_norm', SAME_AS, 0, 5),
        Constraint('dwd_norm', SAME_AS, 0, 5),
        Constraint('shiftscale_norm', SAME_AS, 0, 5),
        Constraint('bfrm_norm', SAME_AS, 0, 5),
        Constraint('gene_center', SAME_AS, 0, 5),
        Constraint('gene_normalize', SAME_AS, 0, 5),
        Constraint('annotate', SAME_AS, 0, 5),
        Constraint('relabel_sample', SAME_AS, 0, 5),
        Constraint('platform', SAME_AS, 0, 5),
        Constraint('num_features', SAME_AS, 0, 5),
        Constraint('filter_and_threshold', MUST_BE, "yes", 5),
        Constraint('unique_genes', SAME_AS, 0, 5),
        Constraint('duplicate_probe', SAME_AS, 0, 5),
        Constraint('group_fc', SAME_AS, 0, 5),

        # PCAPlot (normalized).
        Constraint('logged', MUST_BE, "yes", 6),
        Constraint("contents", SAME_AS, 0, 6),
        Constraint("preprocess", SAME_AS, 0, 6),
        Constraint('quantile_norm', SAME_AS, 1, 6),
        Constraint('combat_norm', SAME_AS, 1, 6),
        Constraint('dwd_norm', SAME_AS, 1, 6),
        Constraint('shiftscale_norm', SAME_AS, 1, 6),
        Constraint('bfrm_norm', SAME_AS, 1, 6),
        Constraint('gene_center', SAME_AS, 1, 6),
        Constraint('gene_normalize', SAME_AS, 1, 6),
        Constraint('annotate', SAME_AS, 0, 6),
        Constraint('relabel_sample', SAME_AS, 0, 6),
        Constraint('platform', SAME_AS, 0, 6),
        Constraint('num_features', SAME_AS, 0, 6),
        Constraint('filter_and_threshold', MUST_BE, "yes", 6),
        Constraint('unique_genes', SAME_AS, 0, 6),
        Constraint('duplicate_probe', SAME_AS, 0, 6),
        Constraint('group_fc', SAME_AS, 0, 6),

        # Heatmap (no normalization).
        Constraint('logged', MUST_BE, "yes", 7),
        Constraint("contents", SAME_AS, 0, 7),
        Constraint("preprocess", SAME_AS, 0, 7),
        Constraint('quantile_norm', SAME_AS, 0, 7),
        Constraint('combat_norm', SAME_AS, 0, 7),
        Constraint('dwd_norm', SAME_AS, 0, 7),
        Constraint('shiftscale_norm', SAME_AS, 0, 7),
        Constraint('bfrm_norm', SAME_AS, 0, 7),
        Constraint("gene_center", MUST_BE, "mean", 7),
        Constraint("gene_normalize", MUST_BE, "variance", 7),
        Constraint('annotate', SAME_AS, 0, 7),
        Constraint('relabel_sample', SAME_AS, 0, 7),
        Constraint('platform', SAME_AS, 0, 7),
        Constraint('num_features', SAME_AS, 0, 7),
        Constraint("filter_and_threshold", MUST_BE, "yes", 7),
        Constraint('unique_genes', SAME_AS, 0, 7),
        Constraint('duplicate_probe', SAME_AS, 0, 7),
        Constraint('group_fc', SAME_AS, 0, 7),
        Constraint("cluster_alg", MUST_BE, "hierarchical", 7),

        # Heatmap (normalized).
        Constraint('logged', MUST_BE, "yes", 8),
        Constraint("contents", SAME_AS, 0, 8),
        Constraint("preprocess", SAME_AS, 0, 8),
        Constraint('quantile_norm', SAME_AS, 1, 8),
        Constraint('combat_norm', SAME_AS, 1, 8),
        Constraint('dwd_norm', SAME_AS, 1, 8),
        Constraint('shiftscale_norm', SAME_AS, 1, 8),
        Constraint('bfrm_norm', SAME_AS, 1, 8),
        Constraint("gene_center", MUST_BE, "mean", 8),
        Constraint("gene_normalize", MUST_BE, "variance", 8),
        Constraint('annotate', SAME_AS, 0, 8),
        Constraint('relabel_sample', SAME_AS, 0, 8),
        Constraint('platform', SAME_AS, 0, 8),
        Constraint('num_features', SAME_AS, 0, 8),
        Constraint("filter_and_threshold", MUST_BE, "yes", 8),
        Constraint('unique_genes', SAME_AS, 0, 8),
        Constraint('duplicate_probe', SAME_AS, 0, 8),
        Constraint('group_fc', SAME_AS, 0, 8),
        Constraint("cluster_alg", MUST_BE, "hierarchical", 8),
        
        Consequence("logged", SAME_AS_CONSTRAINT, 0),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT, 0),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT, 0),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("gene_center", SAME_AS_CONSTRAINT, 1),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT, 1),
        Consequence("gene_order", SAME_AS_CONSTRAINT, 0),
        Consequence("annotate", SAME_AS_CONSTRAINT, 0),
        Consequence("relabel_sample", SAME_AS_CONSTRAINT, 0),
        Consequence("platform", SAME_AS_CONSTRAINT, 0),
        Consequence("num_features", SAME_AS_CONSTRAINT, 0),
        Consequence("filter_and_threshold", SAME_AS_CONSTRAINT, 0),
        Consequence("unique_genes", SAME_AS_CONSTRAINT, 0),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT, 0),
        Consequence("group_fc", SAME_AS_CONSTRAINT, 0),

        help="Do comprehensive preprocessing for a gene expression matrix.",
        ),
    
    ModuleNode(
        'do_complete_preprocessing_illumina',
        [
            SignalFile.SignalFile,           #  0   no normalization
            SignalFile.SignalFile,           #  1   normalized
            EV.SignalDistributionBoxplot,    #  2
            EV.ActbPlot,                     #  3   no normalization
            EV.ActbPlot,                     #  4   normalized
            EV.PCAPlot,                      #  5   no normalization
            EV.PCAPlot,                      #  6   normalized
            EV.Heatmap,                      #  7   no normalization
            EV.Heatmap,                      #  8   normalized
            SignalFile.IlluminaControlFile,  #  9
            EV.HousekeepingPlot,             # 10
            EV.BiotinPlot,                   # 11
            EV.IlluHybridizationProbePlot,   # 12
            ],
        CompleteExpressionPreprocessing,

        # SignalFile (not normalized)
        Constraint("logged", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("preprocess", MUST_BE, "illumina", 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("filter_missing_values", CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"], 0),
        Constraint('quantile_norm', MUST_BE, "no", 0),
        Constraint('combat_norm', MUST_BE, "no", 0),
        Constraint('dwd_norm', MUST_BE, "no", 0),
        Constraint('shiftscale_norm', MUST_BE, "no", 0),
        Constraint('bfrm_norm', MUST_BE, "no", 0),
        Constraint('gene_center', MUST_BE, "no", 0),
        Constraint('gene_normalize', MUST_BE, "no", 0),
        Constraint("gene_order", CAN_BE_ANY_OF, BDT.GENE_ORDER, 0),
        Constraint('annotate', CAN_BE_ANY_OF, YESNO, 0),
        Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("platform", CAN_BE_ANY_OF, ["no", "yes", 'u133A'], 0),
        Constraint('num_features', CAN_BE_ANY_OF, YESNO, 0),
        Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            'unique_genes', CAN_BE_ANY_OF,
            ['no', 'average_genes', 'high_var', 'first_gene'], 0),
        Constraint(
            'duplicate_probe', CAN_BE_ANY_OF,
            ["no", "closest_probe", "high_var_probe"], 0),
        Constraint('group_fc', CAN_BE_ANY_OF, YESNO, 0),

        # SignalFile (normalized)
        Constraint('logged', SAME_AS, 0, 1),
        Constraint("contents", SAME_AS, 0, 1),
        Constraint("preprocess", SAME_AS, 0, 1),
        Constraint("filter_missing_values", SAME_AS, 0, 1),
        Constraint("missing_algorithm", SAME_AS, 0, 1),
        Constraint('quantile_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('combat_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('dwd_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('shiftscale_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('bfrm_norm', CAN_BE_ANY_OF, YESNO, 1),
        Constraint('gene_center', CAN_BE_ANY_OF, ['median', 'mean', 'no'], 1),
        Constraint(
            'gene_normalize', CAN_BE_ANY_OF,
            ['variance', 'sum_of_squares', 'no'], 1),
        Constraint("gene_order", SAME_AS, 0, 1),
        Constraint("annotate", SAME_AS, 0, 1),
        Constraint("relabel_sample", SAME_AS, 0, 1),
        Constraint("platform", SAME_AS, 0, 1),
        Constraint('num_features', SAME_AS, 0, 1),
        Constraint("filter_and_threshold", SAME_AS, 0, 1),
        Constraint('unique_genes', SAME_AS, 0, 1),
        Constraint('duplicate_probe', SAME_AS, 0, 1),
        Constraint('group_fc', SAME_AS, 0, 1),

        # SignalDistributionBoxplot
        Constraint("logged", MUST_BE, "yes", 2),
        Constraint("contents", SAME_AS, 0, 2),
        Constraint("preprocess", SAME_AS, 0, 2),
        Constraint('relabel_sample', SAME_AS, 0, 2),

        # ActbPlot (No normalization)
        Constraint('logged', MUST_BE, "yes", 3),
        Constraint("preprocess", SAME_AS, 0, 3),
        Constraint("contents", SAME_AS, 0, 3),
        Constraint('quantile_norm', SAME_AS, 0, 3),
        Constraint('combat_norm', SAME_AS, 0, 3),
        Constraint('dwd_norm', SAME_AS, 0, 3),
        Constraint('shiftscale_norm', SAME_AS, 0, 3),
        Constraint('bfrm_norm', SAME_AS, 0, 3),
        Constraint('gene_center', SAME_AS, 0, 3),
        Constraint('gene_normalize', SAME_AS, 0, 3),
        Constraint('annotate', MUST_BE, "yes", 3),
        Constraint('relabel_sample', SAME_AS, 0, 3),
        Constraint('platform', SAME_AS, 0, 3),
        Constraint('num_features', SAME_AS, 0, 3),
        Constraint('filter_and_threshold', MUST_BE, "no", 3),
        Constraint('unique_genes', MUST_BE, "no", 3),
        Constraint('duplicate_probe', MUST_BE, "no", 3),
        Constraint('group_fc', MUST_BE, "no", 3),

        # ActbPlot (Normalized)
        Constraint('logged', MUST_BE, "yes", 4),
        Constraint("contents", SAME_AS, 0, 4),
        Constraint("preprocess", SAME_AS, 0, 4),
        Constraint('quantile_norm', SAME_AS, 1, 4),
        Constraint('combat_norm', SAME_AS, 1, 4),
        Constraint('dwd_norm', SAME_AS, 1, 4),
        Constraint('shiftscale_norm', SAME_AS, 1, 4),
        Constraint('bfrm_norm', SAME_AS, 1, 4),
        Constraint('gene_center', SAME_AS, 1, 4),
        Constraint('gene_normalize', SAME_AS, 1, 4),
        Constraint('annotate', MUST_BE, "yes", 4),
        Constraint('relabel_sample', SAME_AS, 0, 4),
        Constraint('platform', SAME_AS, 0, 4),
        Constraint('num_features', SAME_AS, 0, 4),
        Constraint('filter_and_threshold', MUST_BE, "no", 4),
        Constraint('unique_genes', SAME_AS, 0, 4),
        Constraint('duplicate_probe', SAME_AS, 0, 4),
        Constraint('group_fc', SAME_AS, 0, 4),

        # PcaPlot (No normalization)
        Constraint("logged", MUST_BE, "yes", 5),
        Constraint("contents", SAME_AS, 0, 5),
        Constraint("preprocess", SAME_AS, 0, 5),
        Constraint('quantile_norm', SAME_AS, 0, 5),
        Constraint('combat_norm', SAME_AS, 0, 5),
        Constraint('dwd_norm', SAME_AS, 0, 5),
        Constraint('shiftscale_norm', SAME_AS, 0, 5),
        Constraint('bfrm_norm', SAME_AS, 0, 5),
        Constraint('gene_center', SAME_AS, 0, 5),
        Constraint('gene_normalize', SAME_AS, 0, 5),
        Constraint('annotate', SAME_AS, 0, 5),
        Constraint('relabel_sample', SAME_AS, 0, 5),
        Constraint('platform', SAME_AS, 0, 5),
        Constraint('num_features', SAME_AS, 0, 5),
        Constraint('filter_and_threshold', MUST_BE, "yes", 5),
        Constraint('unique_genes', SAME_AS, 0, 5),
        Constraint('duplicate_probe', SAME_AS, 0, 5),
        Constraint('group_fc', SAME_AS, 0, 5),

        # PCAPlot (normalized).
        Constraint('logged', MUST_BE, "yes", 6),
        Constraint("contents", SAME_AS, 0, 6),
        Constraint("preprocess", SAME_AS, 0, 6),
        Constraint('quantile_norm', SAME_AS, 1, 6),
        Constraint('combat_norm', SAME_AS, 1, 6),
        Constraint('dwd_norm', SAME_AS, 1, 6),
        Constraint('shiftscale_norm', SAME_AS, 1, 6),
        Constraint('bfrm_norm', SAME_AS, 1, 6),
        Constraint('gene_center', SAME_AS, 1, 6),
        Constraint('gene_normalize', SAME_AS, 1, 6),
        Constraint('annotate', SAME_AS, 0, 6),
        Constraint('relabel_sample', SAME_AS, 0, 6),
        Constraint('platform', SAME_AS, 0, 6),
        Constraint('num_features', SAME_AS, 0, 6),
        Constraint('filter_and_threshold', MUST_BE, "yes", 6),
        Constraint('unique_genes', SAME_AS, 0, 6),
        Constraint('duplicate_probe', SAME_AS, 0, 6),
        Constraint('group_fc', SAME_AS, 0, 6),

        # Heatmap (no normalization).
        Constraint('logged', MUST_BE, "yes", 7),
        Constraint("contents", SAME_AS, 0, 7),
        Constraint("preprocess", SAME_AS, 0, 7),
        Constraint('quantile_norm', SAME_AS, 0, 7),
        Constraint('combat_norm', SAME_AS, 0, 7),
        Constraint('dwd_norm', SAME_AS, 0, 7),
        Constraint('shiftscale_norm', SAME_AS, 0, 7),
        Constraint('bfrm_norm', SAME_AS, 0, 7),
        Constraint("gene_center", MUST_BE, "mean", 7),
        Constraint("gene_normalize", MUST_BE, "variance", 7),
        Constraint('annotate', SAME_AS, 0, 7),
        Constraint('relabel_sample', SAME_AS, 0, 7),
        Constraint('platform', SAME_AS, 0, 7),
        Constraint('num_features', SAME_AS, 0, 7),
        Constraint("filter_and_threshold", MUST_BE, "yes", 7),
        Constraint('unique_genes', SAME_AS, 0, 7),
        Constraint('duplicate_probe', SAME_AS, 0, 7),
        Constraint('group_fc', SAME_AS, 0, 7),
        Constraint("cluster_alg", MUST_BE, "hierarchical", 7),

        # Heatmap (normalized).
        Constraint('logged', MUST_BE, "yes", 8),
        Constraint("contents", SAME_AS, 0, 8),
        Constraint("preprocess", SAME_AS, 0, 8),
        Constraint('quantile_norm', SAME_AS, 1, 8),
        Constraint('combat_norm', SAME_AS, 1, 8),
        Constraint('dwd_norm', SAME_AS, 1, 8),
        Constraint('shiftscale_norm', SAME_AS, 1, 8),
        Constraint('bfrm_norm', SAME_AS, 1, 8),
        Constraint("gene_center", MUST_BE, "mean", 8),
        Constraint("gene_normalize", MUST_BE, "variance", 8),
        Constraint('annotate', SAME_AS, 0, 8),
        Constraint('relabel_sample', SAME_AS, 0, 8),
        Constraint('platform', SAME_AS, 0, 8),
        Constraint('num_features', SAME_AS, 0, 8),
        Constraint("filter_and_threshold", MUST_BE, "yes", 8),
        Constraint('unique_genes', SAME_AS, 0, 8),
        Constraint('duplicate_probe', SAME_AS, 0, 8),
        Constraint('group_fc', SAME_AS, 0, 8),
        Constraint("cluster_alg", MUST_BE, "hierarchical", 8),

        # IlluminaControlFile
        Constraint("contents", SAME_AS, 0, 9),
        
        # HousekeepingPlot
        Constraint("contents", SAME_AS, 0, 10),
        Constraint('relabel_sample', SAME_AS, 0, 10),

        # BiotinPlot
        Constraint("contents", SAME_AS, 0, 11),
        Constraint('relabel_sample', SAME_AS, 0, 11),

        # IlluHybridizationProbePlot
        Constraint("contents", SAME_AS, 0, 12),
        

        Consequence("logged", SAME_AS_CONSTRAINT, 0),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT, 0),
        Consequence("filter_missing_values", SAME_AS_CONSTRAINT, 0),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 1),
        Consequence("gene_center", SAME_AS_CONSTRAINT, 1),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT, 1),
        Consequence("gene_order", SAME_AS_CONSTRAINT, 0),
        Consequence("annotate", SAME_AS_CONSTRAINT, 0),
        Consequence("relabel_sample", SAME_AS_CONSTRAINT, 0),
        Consequence("platform", SAME_AS_CONSTRAINT, 0),
        Consequence("num_features", SAME_AS_CONSTRAINT, 0),
        Consequence("filter_and_threshold", SAME_AS_CONSTRAINT, 0),
        Consequence("unique_genes", SAME_AS_CONSTRAINT, 0),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT, 0),
        Consequence("group_fc", SAME_AS_CONSTRAINT, 0),

        help="Do comprehensive preprocessing for illumina",
        ),
    
    ModuleNode(
        'make_preprocessing_report',
        CompleteExpressionPreprocessing, ExpressionPreprocessingReport,
        ),

    ModuleNode(
        "make_batch_effect_report",
        [
            SignalFile.SignalFile,
            EV.PCAPlot,
         
            SignalFile.SignalFile,
            EV.PCAPlot,
         
            SignalFile.SignalFile,
            EV.PCAPlot,
         
            SignalFile.SignalFile,
            EV.PCAPlot,
         
            SignalFile.SignalFile,
            EV.PCAPlot,
         
            SignalFile.SignalFile,
            EV.PCAPlot
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
        "convert_complete_expression_preprocessing_preprocess",
        CompleteExpressionPreprocessing, CompleteExpressionPreprocessing,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SET_TO, "any"),
        Constraint("annotate", MUST_BE, "yes"),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        ),
    
    ModuleNode(
        'convert_preprocessing_report_preprocess',
        ExpressionPreprocessingReport, ExpressionPreprocessingReport,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SET_TO, "any"),
        help="convert preprocess from others to any in PreprocessingReport",
        ),
    ]
