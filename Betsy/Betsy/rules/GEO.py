from Betsy.bie3 import *
import BasicDataTypes as BDT

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
    help="GEO fmaily soft file download from the GEO database")

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
    Module(
        "download_geo_supplement", GEOSeries, BDT.ExpressionFiles,
        OptionDef("GSEID"),
        OptionDef("GPLID", ""),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("filetype", SET_TO, 'unknown'),
        help="Download the supplemental expression files from GEO.",
        ),
     Module(
        "check_geo_file_type", [BDT.ExpressionFiles, GEOSeriesMatrixFile],
        BDT.ExpressionFiles,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("filetype", MUST_BE, 'unknown', 0),
        Consequence(
            "filetype", BASED_ON_DATA,
            ['matrix', 'cel', 'gpr', 'idat', 'agilent']),
        help="check the file type downloaded from geo database"
        ),
    Module(
        "download_geo_seriesmatrix", GEOSeries, GEOSeriesMatrixFile,
        OptionDef("GSEID"),
        OptionDef("GPLID", ""),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="Download SeriesMatrix file from GEO."
        ),
    Module(
        "get_geo_sample_metadata", GEOSeriesMatrixFile, GEOSampleMetadata,
        help="Get the metadata for the samples for a GEO data set.",
        ),

    Module(
        'download_GEO_family_soft', GEOSeries, GEOFamilySoftFile,
        OptionDef("GSEID", help='GSEID for download family_soft file'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="download geo family soft file"),
    Module(
        'convert_family_soft_to_rename', GEOFamilySoftFile, BDT.RenameFile,
        OptionDef("GSEID", help='GSEID for download family_soft file'),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("labels_from", SET_TO_ONE_OF, ["title","description"]),
        help="convert famliy soft file to RenameFile"),
    ]
