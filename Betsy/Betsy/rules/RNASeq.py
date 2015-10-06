#RNASeq
from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
import GeneExpProcessing as GXP

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

RSEMIndexedGenome = DataType(
    "RSEMIndexedGenome",
    help="Indexed for rsem.",
    )
RSEMResults = DataType(
    "RSEMResults",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    help="Results from an rsem-calculate-expression analysis.",
    )

all_data_types = [
    RSEMIndexedGenome,
    RSEMResults,
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
        "index_rsem_reference",
        NGS.ReferenceGenome, RSEMIndexedGenome,
        OptionDef(
            "assembly", default="genome",
            help="Optional name for the genome assembly, e.g. hg19",
            ),
        OptionDef(
            "gtf_file", 
            help="Gene annotations in GTF format.",
            ),
        ),
    
    ModuleNode(
        "normalize_with_rsem",
        [NGS.FastqFolder, NGS.SampleGroupFile, RSEMIndexedGenome],
        #GXP.UnprocessedSignalFile,
        RSEMResults,
        Constraint(
            "orientation", CAN_BE_ANY_OF, NGS.ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        help="Use RSEM to estimate TPM or FPKM.",
        ),

    ModuleNode(
        "extract_rsem_signal", RSEMResults, GXP.UnprocessedSignalFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        Consequence("preprocess", SET_TO_ONE_OF, ["tpm", "fpkm"]),
        Consequence("logged", SET_TO, "no"),
        # What is this for?
        #Consequence("predataset", SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        ),
    ]
