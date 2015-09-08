# SampleGroupFile
# FastqFile
# FastqFolder
# SamFile
# SamFolder
# BamFile
# BamFolder
#
# SaiFile
# VcfFile

from Betsy.bie3 import *
import BasicDataTypes as BDT


SampleGroupFile = DataType(
    "SampleGroupFile",
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified", help="contents"),
    help="File contains sample group infomation"
    )

FastqFile = DataType(
    "FastqFile",
    AttributeDef(
        "read", ["single", "pair", "pair1", "pair2"],
        "single", "single", help="single or pair read"),
    AttributeDef(
        "ref", ["hg18", "hg19", "mm9", "dm3"], "hg19", "hg19",
        help="ref species"),
    AttributeDef(
        "contents", BDT.CONTENTS,
        "unspecified", "unspecified", help="contents"),
    help="Fastq file"
    )

FastqFolder = DataType(
    "FastqFolder",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help=""),
    AttributeDef(
        "compressed", ["yes", "no", "unknown"], "unknown", "yes",
        help="Whether the files are compressed (gz, bz2, xz)."),
    help="RNA seq Fastq folder"
    )

SamFile = DataType(
    "SamFile",
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified",
                 help="contents"),
    AttributeDef("sorted", ["yes", "no"], "no", "no",
                 help="sorted or not"),
    AttributeDef("duplicates_marked", ["yes", "no"], "no", "no",
                 help="mark duplicate or not"),
    AttributeDef("recalibration", ["yes", "no"], "no", "no",
                 help="recalibration or not"),
    AttributeDef("has_header", ["yes", "no"], "no", "no",
                 help="fix header or not"),
    AttributeDef("read", ["single", "pair"],
                 "single", "single",
                 help="single or pair read"),
    AttributeDef("ref", ["hg18", "hg19", "mm9", "dm3"],
                 "hg19", "hg19",
                 help="ref species"),
    help="Sam file"
    )

SamFolder = DataType(
    "SamFolder",
    AttributeDef("ref", ["human", "mouse", "hg18", "hg19"], "human", "human",
                 help="ref species"),
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified", help="contents"),
    AttributeDef("sample_type", ["RNA", "DNA"],
                 "RNA", "RNA", help="RNA or DNA type"),
    help="RNA seq Sam folder"
    )

BamFile = DataType(
    "BamFile",
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified",
                 help="contents"),
    AttributeDef("sorted", ["yes", "no"], "no", "no",
                 help="sorted or not"),
    AttributeDef("duplicates_marked", ["yes", "no"], "no", "no",
                 help="mark duplicate or not"),
    AttributeDef("recalibration", ["yes", "no"], "no", "no",
                 help="recalibration or not"),
    AttributeDef("has_header", ["yes", "no"], "no", "no",
                 help="fix header or not"),
    AttributeDef("read", ["single", "pair"],
                 "single", "single",
                 help="single or pair read"),
    AttributeDef("ref", ["hg18", "hg19", "mm9", "dm3"],
                 "hg19", "hg19",
                 help="ref species"),
    help="Bam file"
    )

BamFolder = DataType(
    "BamFolder",
    AttributeDef("ref", ["human", "mouse", "hg18", "hg19"], "human", "human",
                 help="ref species"),
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified", help="contents"),
    AttributeDef("sample_type", ["RNA", "DNA"],
                 "RNA", "RNA", help="RNA or DNA type"),
    AttributeDef("duplicates_marked", ["yes", "no"], "no", "no",
                 help="mark duplicate or not"),
    AttributeDef("sorted", ["yes", "no", "unknown"], "unknown", "unknown",
                 help="sorted or not"),
    AttributeDef("indexed", ["yes", "no"], "no", "no",
                 help="indexed or not"),
    help="RNA seq Bam folder"
    )

SaiFile = DataType(
    "SaiFile",
    AttributeDef("read", ["single", "pair", "pair1", "pair2"],
                 "single", "single", help="single or pair read"),
    AttributeDef("ref", ["hg18", "hg19", "mm9", "dm3"],
                 "hg19", "hg19", help="ref species"),
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified", help="contents"),
    help="Sai file"
    )

VcfFile = DataType(
    "VcfFile",
    AttributeDef("contents", BDT.CONTENTS,
                 "unspecified", "unspecified",
                 help="contents"),
    AttributeDef("recalibration", ["yes", "no"], "no", "no",
                 help="recalibration or not"),
    AttributeDef("read", ["single", "pair"],
                 "single", "single",
                 help="single or pair read"),
    AttributeDef("ref", ["hg18", "hg19", "mm9", "dm3"],
                 "hg19", "hg19",
                 help="ref species"),
    AttributeDef("vcf_filter", ["yes", "no"], "no", "no",
                 help="filter VcfFile or not"),
    AttributeDef("reheader", ["standard", "bcftool"],
                 "standard", "standard",
                 help="method to convert to VcfFile"),
    AttributeDef("vcf_annotate", ["yes", "no"], "no", "no",
                 help="annotate VcfFile or not"),
    help="Vcf file"
    )

all_data_types = [
    SampleGroupFile,
    FastqFile,
    FastqFolder,
    SamFile,
    SamFolder,
    BamFile,
    BamFolder,
    SaiFile,
    VcfFile,
    ]

all_modules = [
    ModuleNode(
        "is_fastq_folder_compressed",
        FastqFolder, FastqFolder,
        Constraint("compressed", MUST_BE, "unknown"),
        Consequence("compressed", BASED_ON_DATA, ["no", "yes"]),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "uncompress_fastq_folder",
        FastqFolder, FastqFolder,
        Constraint("compressed", MUST_BE, "yes"),
        Consequence("compressed", SET_TO, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    

    ## ModuleNode(
    ##     "is_fastq_folder",
    ##     RNASeqFile, RNASeqFile,
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("format_type", MUST_BE, "unknown"),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence(
    ##         "format_type", BASED_ON_DATA, ["not_fastqfolder", "fastqfolder"]),
    ##     help=("extract rna files with different format")
    ##     ),
    ## ModuleNode(
    ##     "is_sam_folder",
    ##     RNASeqFile, RNASeqFile,
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("format_type", MUST_BE, "not_fastqfolder"),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence(
    ##         "format_type", BASED_ON_DATA, ["not_samfolder", "samfolder"]),
    ##     help=("extract rna files with different format")
    ##     ),
    ## ModuleNode(
    ##     "is_bam_folder",
    ##     RNASeqFile, RNASeqFile,
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("format_type", MUST_BE, "not_samfolder"),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence(
    ##         "format_type", BASED_ON_DATA, ["not_bamfolder", "bamfolder"]),
    ##     help=("extract rna files with different format")
    ##     ),
    ModuleNode(
        "convert_sam_to_bam",
        SamFolder, BamFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        #Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
        help="Convert SAM to BAM files.",
        ),
    ]
