from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS

FastQCFolder = DataType(
    "FastQCFolder",
    help="Folder that holds FastQC results.",
    )

RNASeQCFile = DataType(
    "RNASeQCFile",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    help="File contains sample group infomation"
    )

all_data_types = [
    FastQCFolder,
    RNASeQCFile,
    ]

all_modules = [
    ModuleNode(
        "run_fastqc",
        NGS.FastqFolder, FastQCFolder,
        Constraint("compressed", MUST_BE, "no"),
        help="Run FastQC on a file of fastq files",
        ),
    ModuleNode(
        "run_RNA_SeQC",
        NGS.BamFolder, RNASeQCFile,
        OptionDef(
            "RNA_ref", help="ref file for RNA_SeQC"),
        OptionDef(
            "RNA_gtf", help="gtf file for RNA_SeQC"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
        Constraint("duplicates_marked",MUST_BE, "yes"),
        Constraint("indexed",MUST_BE, "yes"),
        Constraint("sorted", MUST_BE, "yes"),
        Constraint("sample_type", MUST_BE, "RNA"),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="run RNA-SeQC"),
    ]
