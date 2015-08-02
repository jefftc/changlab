#ArrayPlatforms
from Betsy.bie3 import *
import GeneExpProcessing
import BasicDataTypes as BDT

ILLU_MANIFEST = [
    'HumanHT-12_V3_0_R2_11283641_A.txt',
    'HumanHT-12_V4_0_R2_15002873_B.txt',
    'HumanHT-12_V3_0_R3_11283641_A.txt',
    'HumanHT-12_V4_0_R1_15002873_B.txt',
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
    'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
    'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
    'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
    'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
    'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
    'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
    'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
    'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
    'ilmn_RatRef_12_V1_0_R5_11222119_A.chip',
    ]

AgilentFiles = DataType(
    "AgilentFiles",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="A folder of agilent files.")

CELFiles = DataType(
    "CELFiles",
    AttributeDef(
        "version", ["unknown", "cc", "v3_v4"], "unknown", "v3_v4",
        help="cel file version"),
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="A folder of cel files")

ControlFile = DataType(
    "ControlFile",
    AttributeDef(
        'preprocess', ["illumina"], "illumina", "illumina",
        help="preprocess for ControlFile"),
    AttributeDef(
        'missing_values', ["unknown", "no", "yes"], "no", "no",
        help="missing values yes or not"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill", help="missing algorithm"),
    AttributeDef("logged", ["no"], "no", "no", help="logged yes or not"),
    AttributeDef('format', ["gct"], "gct", "gct", help="file format"),
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="The control file in the ILLUFolder")

GPRFiles = DataType(
    "GPRFiles",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="A folder of GPR files.")

IDATFiles = DataType(
    "IDATFiles",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="A folder of IDAFiles.")

ILLUFolder = DataType(
    "ILLUFolder",
    AttributeDef(
        "illu_manifest", ILLU_MANIFEST, 'HumanHT-12_V4_0_R2_15002873_B.txt',
        'HumanHT-12_V4_0_R2_15002873_B.txt', help="illumina manifest"),
    AttributeDef(
        'illu_chip', ILLU_CHIP, 'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
        'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
        help="illumina chip type"),
    AttributeDef(
        'illu_bg_mode', ['false', 'true'], "false", "false",
        help="illumina background mode"),
    AttributeDef(
        'illu_coll_mode', ['none', 'max', 'median'], "none", "none",
        help="illumina coll mode"),
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="A folder generated from preprocess_illumina," \
    "it contains SignalFile_Postprocess and ControlFile.")

ActbPlot = DataType(
    'ActbPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help='contents'),
    AttributeDef(
        "preprocess", GeneExpProcessing.PREPROCESS, 'unknown', 'unknown',
        help="preprocess method"),
    help="Actb plot file")

Hyb_barPlot = DataType(
    'Hyb_barPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help='contents'),
    help="Hyb bar plot file")

ControlPlot = DataType(
    'ControlPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help='contents'),
    AttributeDef(
        "preprocess", GeneExpProcessing.PREPROCESS, 'unknown', 'unknown',
        help="preprocess method"),
    help="control plot file")

BiotinPlot = DataType(
    'BiotinPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help='contents'),
    help="Biotin plot file")

HousekeepingPlot = DataType(
    'HousekeepingPlot',
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help='contents'),
    help="Housekeeping plot file")

list_files = [
    AgilentFiles,
    CELFiles,
    ControlFile,
    GPRFiles,
    IDATFiles,
    ILLUFolder,
    ActbPlot,
    Hyb_barPlot,
    ControlPlot,
    BiotinPlot,
    HousekeepingPlot,
    ]

all_modules = [
    Module(
        "extract_matrix_data", BDT.ExpressionFiles,
        GeneExpProcessing._SignalFile_Postprocess,
        Constraint("filetype", MUST_BE, 'matrix'),
        Consequence("preprocess", SET_TO, "unknown"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        help="extract SignalFile_Postprocess files from Expression Files"),
    
    #CELFiles
    Module(
        "extract_CEL_files", BDT.ExpressionFiles, CELFiles,
        Constraint("filetype", MUST_BE, 'cel'),
        Consequence("version", SET_TO, "unknown"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract CEL files from Expression Files folder"),
    Module(
        "detect_CEL_version", CELFiles, CELFiles,
        Constraint("version", MUST_BE, "unknown"),
        Consequence("version", BASED_ON_DATA, ["cc", "v3_v4"]),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="detect the version of cel files"),
    Module(
        "convert_CEL_to_v3", CELFiles, CELFiles,
        Constraint("version", MUST_BE, "cc"),
        Consequence("version", SET_TO, "v3_v4"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="convert cel files to v3 version"),
    Module(
        "preprocess_rma", CELFiles, GeneExpProcessing._SignalFile_Postprocess,
        Constraint("version", MUST_BE, 'v3_v4'),
        Consequence("logged", SET_TO, "yes"),
        Consequence("preprocess", SET_TO, "rma"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence('predataset', SET_TO, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help=
        "preprocess CELFiles with rma method,generate SignalFile_Postprocess"),
    Module(
        "preprocess_mas5", CELFiles, GeneExpProcessing._SignalFile_Postprocess,
        Constraint("version", MUST_BE, 'v3_v4'),
        Consequence("logged", SET_TO, "no"),
        Consequence("preprocess", SET_TO, "mas5"),
        Consequence('predataset', SET_TO, "no"),
        Consequence("format", SET_TO, "jeffs"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="preprocess CELFiles with mas5 method,generate "
        "SignalFile_Postprocess",
        ),
    # IDATFiles
    Module(
        "extract_illumina_idat_files", BDT.ExpressionFiles, IDATFiles,
        Constraint("filetype", MUST_BE, 'idat'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract idat files from Expression File folder"),
    Module(
        "preprocess_illumina", IDATFiles, ILLUFolder,
        Consequence('illu_manifest', SET_TO_ONE_OF, ILLU_MANIFEST),
        Consequence('illu_chip', SET_TO_ONE_OF, ILLU_CHIP),
        Consequence('illu_bg_mode', SET_TO_ONE_OF, ["false", "true"]),
        Consequence(
            'illu_coll_mode', SET_TO_ONE_OF, ["none", "max", "median"]),
        
        OptionDef("illu_clm", '', help="illumina clm"),
        OptionDef("illu_custom_chip", '', help="illumina custom chip name"),
        OptionDef(
            "illu_custom_manifest", '', help='illumina custrom manifest file'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="preprocess idat files,generate SignalFile_Postprocess"),
    
    Module(
        "get_illumina_control", ILLUFolder, ControlFile,
        Consequence('preprocess', SET_TO, "illumina"),
        Consequence("format", SET_TO, "gct"),
        Consequence("logged", SET_TO, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract illumina ControlFile from ILLUFolder"),
    Module(
        "get_illumina_signal", ILLUFolder,
        GeneExpProcessing._SignalFile_Postprocess,
        Consequence('preprocess', SET_TO, "illumina"),
        Consequence('format', SET_TO, "gct"),
        Consequence('logged', SET_TO, "no"),
        Consequence('predataset', SET_TO, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract the SignalFile_Postprocess from ILLUFolder"),
    # AgilentFiles
    Module(
        "extract_agilent_files", BDT.ExpressionFiles, AgilentFiles,
        Constraint("filetype", MUST_BE, 'agilent'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract agilent files from ExpressionFiles"),
    Module(
        "preprocess_agilent", AgilentFiles,
        GeneExpProcessing._SignalFile_Postprocess,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence('logged', SET_TO, "unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence('preprocess', SET_TO, "agilent"),
        Consequence('format', SET_TO, "tdf"),
        help="preprocess agilent, generate SignalFile_Postprocess"),
    # GPRFiles
    Module(
        "extract_gpr_files", BDT.ExpressionFiles, GPRFiles,
        Constraint("filetype", MUST_BE, 'gpr'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="extract gpr files from ExpressionFiles"),
    Module(
        "normalize_with_loess", GPRFiles,
        GeneExpProcessing._SignalFile_Postprocess,
        Consequence("format", SET_TO, "tdf"),
        Consequence("logged", SET_TO, "unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence("preprocess", SET_TO, "loess"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="normalize GPRFiles,generate SignalFile_Postprocess"),
    Module(
        'plot_actb_line', GeneExpProcessing._SignalFile_Impute, ActbPlot,
        Constraint("preprocess", CAN_BE_ANY_OF, GeneExpProcessing.PREPROCESS1),
        Constraint("missing_values", MUST_BE, 'no'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        help="plot actb line"),
    Module(
        'plot_affy_affx_line', GeneExpProcessing.SignalFile, ControlPlot,
        Constraint("preprocess", CAN_BE_ANY_OF, GeneExpProcessing.PREPROCESS),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot affy affx line"),
    Module(
        'plot_illu_hyb_bar', ControlFile, Hyb_barPlot,
        Constraint("preprocess", MUST_BE, "illumina"),
        Constraint("format", MUST_BE, "gct"),
        Constraint("logged", MUST_BE, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot illumina affx line"),
    Module(
        'plot_illu_biotin_line', ControlFile, BiotinPlot,
        Constraint("preprocess", MUST_BE, 'illumina'),
        Constraint("format", MUST_BE, "gct"),
        Constraint("logged", MUST_BE, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot illumina biotin line"),
    Module(
        'plot_illu_housekeeping_line', ControlFile, HousekeepingPlot,
        Constraint("preprocess", MUST_BE, 'illumina'),
        Constraint("format", MUST_BE, "gct"),
        Constraint("logged", MUST_BE, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot illumina housekeeping line"),
    ]