#RNASeq
from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
import GeneExpProcessing

# TODO: Should move out modules that aren"t specifically for RNA-Seq.
# E.g. is_sam_folder.

#RNASeqFile = DataType(
#    "RNASeqFile",
#    AttributeDef("contents", BDT.CONTENTS,
#                 "unspecified", "unspecified", help="contents"),
#    AttributeDef(
#        "format_type",
#        ["unknown", "not_fastqfolder", "not_samfolder", "not_bamfolder",
#         "samfolder", "bamfolder", "fastqfolder"],
#        "unknown", "unknown", help="format type"),
#    help="RNA Seq File"
#    )

all_data_types = [
#    RNASeqFile,
    ]
all_modules = [
#    ModuleNode(
#        # What does this do?
#        "extract_rna_files_sam",
#        RNASeqFile, NGS.SamFolder,
#        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
#        Constraint("format_type",MUST_BE, "samfolder"),
#        Consequence("contents", SAME_AS_CONSTRAINT),
#        Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
#        help=("extract rna files with different format")
#        ),
#    ModuleNode(
#        # What does this do?
#        "extract_rna_files_bam",
#        RNASeqFile, NGS.BamFolder,
#        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
#        Constraint("format_type",MUST_BE, "bamfolder"),
#        Consequence("contents", SAME_AS_CONSTRAINT),
#        Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
#        help=("extract rna files with bam format")
#        ),
#    ModuleNode(
#        # What does this do?
#        "extract_rna_files_fastq",
#        RNASeqFile, NGS.FastqFolder,
#        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
#        Constraint("format_type",MUST_BE, "fastqfolder"),
#        Consequence("contents", SAME_AS_CONSTRAINT),
#        help=("extract rna files with fa or fastq format")
#        ),
    ModuleNode(
        "normalize_with_rsem",
        [NGS.BamFolder, NGS.SampleGroupFile],
        GeneExpProcessing.UnprocessedSignalFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS,0),
        Constraint("contents", SAME_AS,0,1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        Consequence("preprocess", SET_TO, "RSEM"),
        Consequence("logged", SET_TO, "unknown"),
        Consequence("predataset", SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        help=("process BamFolder, generate SignalFile_Postprocess with preprocess rsem"),
        ),
    ]
