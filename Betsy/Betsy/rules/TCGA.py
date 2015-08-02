from Betsy.bie3 import *
import BasicDataTypes as BDT
import GeneExpProcessing as GEP

TCGAID = DataType(
    "TCGAID",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    help="TCGA ID to download from TCGA database")

TCGAFile = DataType(
    "TCGAFile",
    AttributeDef(
        "contents", BDT.CONTENTS, 'unspecified', 'unspecified',
        help="contents"),
    AttributeDef(
        "preprocess", ['RSEM_genes', 'RSEM_exons',
                       'humanmethylation450', 'mirnaseq',
                       'rppa', 'clinical', 'agilent', 'affymetrix'],
        'RSEM_genes', 'RSEM_genes',help="TCGA data type"),
    AttributeDef(
        "tumor_only", ['yes', 'no'],'no', 'no',
        help="select tumor sample only"),
    help="TCGA file download from TCGA database")

list_files = [
    TCGAID,
    TCGAFile,
    ]

all_modules = [
    Module(
        'download_tcga', TCGAID, TCGAFile,
        OptionDef("disease", help="tcga disease type"),
        OptionDef("date", "", help="date for tcga disease"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence(
            "preprocess", SET_TO_ONE_OF,
            ['RSEM_genes', 'RSEM_exons', 'humanmethylation450', 'mirnaseq',
             'rppa', 'clinical']),
        Consequence("tumor_only", SET_TO, 'no'),
        help="download data from tcga website according to TCGAID"),
    
    Module(
        'download_tcga_agilent', TCGAID, TCGAFile,
        OptionDef("disease", help="tcga disease type"),
        OptionDef("date", "", help="date for tcga disease"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SET_TO, 'agilent'),
        Consequence("tumor_only", SET_TO, 'no'),
        help="download agilent data from tcga website according to TCGAID"),
    
    Module(
        'download_tcga_affymetrix', TCGAID, TCGAFile,
        OptionDef("disease", help="tcga disease type"),
        OptionDef("date", "", help="date for tcga disease"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SET_TO, 'affymetrix'),
        Consequence("tumor_only", SET_TO, 'no'),
        help="download affymetrix data from tcga website according to TCGAID"),

    Module(
        'select_tumor_only', TCGAFile, TCGAFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("tumor_only", MUST_BE, 'no'),
        Consequence("tumor_only", SET_TO, 'yes'),
        help="select the tumor sample only in the TCGAFile"),
    
    Module(
        'preprocess_tcga', TCGAFile, GEP._SignalFile_Postprocess,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("tumor_only", MUST_BE, 'yes'),
        Constraint(
            "preprocess", CAN_BE_ANY_OF,
            ['RSEM_genes', 'RSEM_exons', 'humanmethylation450', 'mirnaseq',
             'rppa', 'clinical', 'affymetrix', 'agilent']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence('logged', SET_TO, "unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence('preprocess', SAME_AS_CONSTRAINT),
        Consequence('format', SET_TO, "tdf"),
        help="preprocess tcga file, generate to SignalFile_Postprocess"),
    ]