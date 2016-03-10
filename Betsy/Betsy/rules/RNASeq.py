#RNASeq

# TophatAlignmentFolder
# TophatAlignmentSummary

# align_with_tophat                Tophat
# extract_tophat_bamfolder

from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
import GeneExpProcessing as GXP

YESNO = BDT.YESNO  # for convenience



TophatAlignmentSummary = DataType(
    "TophatAlignmentSummary",
    help="Summarizes the alignment from tophat (.xls file).",
    )

TophatAlignmentFolder = DataType(
    "TophatAlignmentFolder",
    #AttributeDef(
    #    "contents", BDT.CONTENTS,
    #    "unspecified", "unspecified", help="contents"),
    help="A folder that contains alignments from Tophat.",
    )

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

RSEMResults = DataType(
    "RSEMResults",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    help="A folder of results from an rsem-calculate-expression analysis.",
    )

STARReferenceGenome = DataType(
    "STARReferenceGenome",
    help="Indexed for STAR.",
    )

STARAlignmentFolder = DataType(
    "STARAlignmentFolder",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    AttributeDef(
        "is_subset", YESNO, "no", "no",
        help="Generate a small sample file with a limited number of reads."
        ),
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="Results from a STAR alignment.  Includes SAM files and other stuff.",
    )

HTSeqCountResults = DataType(
    "HTSeqCountResults",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    help="Results from HTSeq-Count.",
    )

HTSeqCountSummary = DataType(
    "HTSeqCountSummary",
    help="Contains a summary of the results from htseq-count as a "
    "tab-delimited text file.",
    )



all_data_types = [
    TophatAlignmentFolder,
    TophatAlignmentSummary,
    RSEMReferenceGenome,
    #FullyIndexedRSEMReferenceGenome,
    RSEMResults,
    STARReferenceGenome,
    STARAlignmentFolder,
    HTSeqCountResults,
    HTSeqCountSummary,
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
        [RSEMReferenceGenome, NGS.GTFGeneModel], RSEMReferenceGenome,
        Constraint("rsem_indexed", MUST_BE, "no", 0),
        Consequence("rsem_indexed", SET_TO, "yes"),
        Consequence("bowtie1_indexed", SET_TO, "yes"),
        Consequence("bowtie2_indexed", SET_TO, "yes"),
        ),
    
    ModuleNode(
        "normalize_with_rsem",
        [NGS.FastqFolder, NGS.SampleGroupFile, NGS.ReadStrandedness,
         RSEMReferenceGenome],
        RSEMResults,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("rsem_indexed", MUST_BE, "yes", 3),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        help="Use RSEM to estimate TPM or FPKM.",
        ),

    ModuleNode(
        "extract_rsem_signal", RSEMResults, GXP.UnprocessedSignalFile,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("preprocess", SET_TO_ONE_OF, ["tpm", "fpkm"]),
        Consequence("logged", SET_TO, "no"),
        # What is this for?
        #Consequence("predataset", SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        ),

    ModuleNode(
        "align_with_tophat",
        [NGS.FastqFolder, NGS.SampleGroupFile, NGS.ReadStrandedness,
         NGS.ReferenceGenome, NGS.GTFGeneModel],
        TophatAlignmentFolder,
        OptionDef(
            "tophat_transcriptome_fa", default="",
            help="FASTA file for Transcriptome, e.g. "
            "hg19.knownGene.fa.  Must be indexed by tophat."
            "Either this or tophat_gtf_file must be given."
            ),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 3),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        help="Align to a reference genome with tophat.",
        ),
    ModuleNode(
        "summarize_tophat_alignment",
        TophatAlignmentFolder, TophatAlignmentSummary,
        help="Summarize the alignment, e.g. number of reads aligned.",
        ),
    ModuleNode(
        "extract_tophat_bamfolder",
        TophatAlignmentFolder, NGS.BamFolder,

        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "tophat"),
        help="Pull the BAM files out of the Tophat result folders."
        ),

    ModuleNode(
        "index_reference_star",
        [NGS.ReferenceGenome, NGS.GTFGeneModel], STARReferenceGenome,
        ),

    ModuleNode(
        "align_with_star",
        [NGS.FastqFolder, NGS.SampleGroupFile, NGS.ReadStrandedness,
         STARReferenceGenome],
        STARAlignmentFolder,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 2),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        #Constraint(
        #    "orientation", CAN_BE_ANY_OF, NGS.ORIENTATION_NOT_UNKNOWN, 0),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        help="Align to a reference genome with star.  "
        "Running with too many processors will kill the machine.  8 is "
        "pretty safe."
        ),

    ModuleNode(
        "extract_star_samfolder",
        STARAlignmentFolder, NGS.SamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "star"),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        help="Pull out the SAM files from the STAR results folder.",
        ),

    ModuleNode(
        "count_with_htseq_count",
        [NGS.BamFolder, NGS.SampleGroupFile, NGS.GTFGeneModel, 
         NGS.ReadStrandedness],
        HTSeqCountResults,
        #OptionDef(
        #    "gtf_file", 
        #    help="Gene annotations in GTF format.",
        #    ),
        OptionDef(
            "htseq_count_mode", default="union",
            help="union, intersection-strict, or intersection-nonempty.  "
            "See htseq-count documentation.  union is recommended.",
            ),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 3),
        # Get error when using coordinate sorting:
        # Maximum alignment buffer size exceeded while pairing SAM alignments.
        #Constraint("sorted", CAN_BE_ANY_OF, ["name", "coordinate"], 0),
        Constraint("sorted", MUST_BE, "name", 0),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        #Constraint("aligner", MUST_BE, "bwa_backtrack", 0),
        #Constraint(
        #    "orientation", CAN_BE_ANY_OF,
        #    ["single", "paired", "paired_fr", "paired_rf"], 0),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        help="Use RSEM to estimate TPM or FPKM.  name sorting the BAM file "
        "is better.  Otherwise, may run into error related to buffer size.",
        ),

    ModuleNode(
        "extract_htseq_count_signal",
        HTSeqCountResults, GXP.UnprocessedSignalFile,
        OptionDef(
            "ignore_htseq_count_errors", default="no",
            help='Whether to ignore errors in the files.  '
            'Should be "yes" or "no".',
            ),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        Consequence("preprocess", SET_TO, "counts"),
        Consequence("logged", SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        ),

    ModuleNode(
        "summarize_htseq_count",
        HTSeqCountResults, HTSeqCountSummary,
        ),

    ModuleNode(
        "convert_counts_to_cpm",
        GXP.UnprocessedSignalFile, GXP.UnprocessedSignalFile,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),

        Constraint("preprocess", MUST_BE, "counts"),
        Consequence("preprocess", SET_TO, "cpm"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("format", SET_TO, "tdf"),
        ),
    
    ]
