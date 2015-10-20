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

RSEMReferenceGenome = DataType(
    "RSEMReferenceGenome",
    AttributeDef(
        "bowtie1_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    AttributeDef(
        "bowtie2_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    AttributeDef(
        "rsem_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    help="Indexed for rsem.",
    )
#FullyIndexedRSEMReferenceGenome = DataType(
#    "FullyIndexedReferenceGenome",
#    help="Used only for indexing a new reference genome."
#    )

RSEMResults = DataType(
    "RSEMResults",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    help="Results from an rsem-calculate-expression analysis.",
    )

all_data_types = [
    RSEMReferenceGenome,
    #FullyIndexedRSEMReferenceGenome,
    RSEMResults,
    ]
all_modules = [
    ModuleNode(
        "is_rsemreference_rsem_indexed",
        RSEMReferenceGenome, RSEMReferenceGenome,
        Constraint("rsem_indexed", MUST_BE, "unknown"),
        Consequence("rsem_indexed", BASED_ON_DATA, ["no", "yes"]),
        ),
    ModuleNode(
        "is_rsemreference_bowtie1_indexed",
        RSEMReferenceGenome, RSEMReferenceGenome,
        Constraint("bowtie1_indexed", MUST_BE, "unknown"),
        Consequence("bowtie1_indexed", BASED_ON_DATA, ["no", "yes"]),
        ),
    ModuleNode(
        "is_rsemreference_bowtie2_indexed",
        RSEMReferenceGenome, RSEMReferenceGenome,
        Constraint("bowtie2_indexed", MUST_BE, "unknown"),
        Consequence("bowtie2_indexed", BASED_ON_DATA, ["no", "yes"]),
        ),
    ModuleNode(
        "index_reference_rsem",
        RSEMReferenceGenome, RSEMReferenceGenome,
        OptionDef(
            "gtf_file", 
            help="Gene annotations in GTF format.",
            ),
        Constraint("rsem_indexed", MUST_BE, "no"),
        Consequence("rsem_indexed", SET_TO, "yes"),
        Constraint("bowtie1_indexed", MUST_BE, "no"),
        Consequence("bowtie1_indexed", SET_TO, "yes"),
        Constraint("bowtie2_indexed", MUST_BE, "no"),
        Consequence("bowtie2_indexed", SET_TO, "yes"),
        ),
    
    ModuleNode(
        "normalize_with_rsem",
        [NGS.FastqFolder, NGS.SampleGroupFile, RSEMReferenceGenome],
        #GXP.UnprocessedSignalFile,
        RSEMResults,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint(
            "orientation", CAN_BE_ANY_OF, NGS.ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("rsem_indexed", MUST_BE, "yes", 2),
        
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
