from Betsy.bie3 import *
import BasicDataTypes as BDT
import SignalFile

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
        #"preprocess", ['RSEM_genes', 'RSEM_exons',
        #               'humanmethylation450', 'mirnaseq',
        #               'rppa', 'clinical', 'agilent', 'affymetrix'],
        # XXX exons
        #"preprocess", ['tpm', 
        #               'humanmethylation450', 'mirnaseq',
        #               'rppa', 'clinical', 'agilent', 'affymetrix'],
        # Probably should be renamed something else.
        "preprocess", ['tpm'], 'tpm', 'tpm', help="TCGA data type"),
    AttributeDef(
        "tumor_only", ['yes', 'no'],'no', 'no',
        help="select tumor sample only"),
    help="TCGA file download from TCGA database")

all_data_types = [
    TCGAID,
    TCGAFile,
    ]

all_modules = [
    ModuleNode(
        'download_tcga', TCGAID, TCGAFile,
        OptionDef("disease", help="Which TCGA disease, e.g. BRCA"),
        OptionDef("date", "", help="Date of data set, e.g. 20140715"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence(
            "preprocess", SET_TO_ONE_OF, ["tpm"],
            # XXX RSEM_exons
            #['tpm', 'humanmethylation450', 'mirnaseq', 'rppa', 'clinical']),
            #['RSEM_genes', 'RSEM_exons', 'humanmethylation450', 'mirnaseq',
            # 'rppa', 'clinical']),
            ),
        Consequence("tumor_only", SET_TO, 'no'),
        help="download data from tcga website according to TCGAID"),
    
    ## ModuleNode(
    ##     'download_tcga_agilent', TCGAID, TCGAFile,
    ##     OptionDef("disease", help="tcga disease type"),
    ##     OptionDef("date", "", help="date for tcga disease"),
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence("preprocess", SET_TO, 'agilent'),
    ##     Consequence("tumor_only", SET_TO, 'no'),
    ##     help="download agilent data from tcga website according to TCGAID"),
    
    ## ModuleNode(
    ##     'download_tcga_affymetrix', TCGAID, TCGAFile,
    ##     OptionDef("disease", help="tcga disease type"),
    ##     OptionDef("date", "", help="date for tcga disease"),
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence("preprocess", SET_TO, 'affymetrix'),
    ##     Consequence("tumor_only", SET_TO, 'no'),
    ##     help="download affymetrix data from tcga website according to TCGAID"),

    ModuleNode(
        'select_tumor_only', TCGAFile, TCGAFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("tumor_only", MUST_BE, 'no'),
        Consequence("tumor_only", SET_TO, 'yes'),
        help="select the tumor sample only in the TCGAFile"),
    
    ModuleNode(
        'preprocess_tcga', TCGAFile, SignalFile.UnprocessedSignalFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("tumor_only", MUST_BE, 'yes'),
        Constraint(
            "preprocess", CAN_BE_ANY_OF, ["tpm"],
            #['tpm', 'humanmethylation450', 'mirnaseq',
            # 'rppa', 'clinical', 'affymetrix', 'agilent']
            ),
            #['RSEM_genes', 'RSEM_exons', 'humanmethylation450', 'mirnaseq',
            # 'rppa', 'clinical', 'affymetrix', 'agilent']),
        Consequence('preprocess', SAME_AS_CONSTRAINT),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence('logged', SET_TO, "unknown"),
        #Consequence('filter_and_threshold', SET_TO, "no"),
        Consequence('format', SET_TO, "tdf"),
        help="preprocess tcga file, generate to SignalFile_Postprocess",
        ),
    ]
