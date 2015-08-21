from Betsy.bie3 import *
import BasicDataTypes as BDT

# TODO: Need code to extract _SignalFile_Postprocess from series
# matrix file.

GEOSeries = DataType(
    "GEOSeries",
    AttributeDef(
        "contents", BDT.CONTENTS,
        'unspecified', 'unspecified', help="contents"),
    help="GEOID to download from the GEO database")

GEOFamilySoftFile = DataType(
    "GEOFamilySoftFile",
    AttributeDef(
        "contents", BDT.CONTENTS,
        'unspecified', 'unspecified',help="contents"),
    help="Family soft file from the GEO database.")

GEOSeriesMatrixFile = DataType(
    "GEOSeriesMatrixFile",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="GEO matrix series file download from the GEO database")

GEOSampleMetadata = DataType(
    "GEOSampleMetadata",
    help="Metadata for the samples in a GEO data set.")

list_files = [
    GEOSeries,
    GEOFamilySoftFile,
    GEOSeriesMatrixFile,
    GEOSampleMetadata,
    ]

all_modules = [
    ModuleNode(
        "download_geo_seriesmatrix", GEOSeries, GEOSeriesMatrixFile,
        OptionDef("GSEID"),
        OptionDef("GPLID", default=""),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="Download SeriesMatrix file from GEO."
        ),
    ModuleNode(
        'download_GEO_family_soft', GEOSeries, GEOFamilySoftFile,
        OptionDef("GSEID", help='GSEID for download family_soft file'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="download geo family soft file"),

    ModuleNode(
        "download_geo_supplement", GEOSeries, BDT.ExpressionFiles,
        OptionDef("GSEID"),
        OptionDef("GPLID", default=""),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("filetype", SET_TO, 'unknown'),
        help="Download the supplemental expression files from GEO.",
        ),

    #ModuleNode(
    #    "acquire_geo_expression_files",
    #    [BDT.ExpressionFiles, GEOSeriesMatrixFile],
    #    BDT.ExpressionFiles,
    #    Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    #    Constraint("contents", SAME_AS, 0, 1),
    #    Consequence("contents", SAME_AS_CONSTRAINT, 0),
    #    Constraint(
    #        "filetype", CAN_BE_ANY_OF,
    #        ["series_matrix", "cel", "gpr", "idat", "agilent"]),
    #    help="Download from GEO the expression files of a specific type.",
    #    ),
    
    ModuleNode(
        "identify_type_of_expression_files",
        BDT.ExpressionFiles, BDT.ExpressionFiles,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("filetype", MUST_BE, "unknown", 0),
        Consequence(
            "filetype", BASED_ON_DATA, ['cel', 'gpr', 'idat', 'agilent']),
        help="Identify the type of expression files in a folder.",
        ),

    ModuleNode(
        "get_geo_sample_metadata", GEOSeriesMatrixFile, GEOSampleMetadata,
        OptionDef("set_NA_to", "NA", help='Convert "NA" to another value.'),
        help="Get the metadata for the samples for a GEO data set.",
        ),

    ModuleNode(
        'convert_family_soft_to_rename', GEOFamilySoftFile, BDT.RenameFile,
        OptionDef("GSEID", help='GSEID for download family_soft file'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("labels_from", SET_TO_ONE_OF, ["title","description"]),
        help="convert famliy soft file to RenameFile"),
    ]
