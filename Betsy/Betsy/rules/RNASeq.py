# DataTypes:
# TophatAlignmentFolder      Tophat
# TophatAlignmentSummary     
# RSEMReferenceGenome        RSEM
# RSEMResults
# STARReferenceGenome        STAR
# STARAlignmentFolder
# HTSeqCountResults          HTSeqCount
# HTSeqCountSummary
#
# RNASeqUnprocessedSignalFile  Contains attributes from RNA-Seq analysis.
# 
# CompleteRNASeqAnalysis     A complete RNA-Seq analysis.
#                            FastQCSummary, FastQCFolder
#                            RSeQCResults
#                            BamFolder
#                            UnprocessedSignalFile (TPM)
#                            UnprocessedSignalFile (TPM, isoforms)
#                            NumAlignedReads (Star alignment)
#                            UnprocessedSignalFile (counts)
#                            HTSeqCountSummary
#
#
# Modules:
# index_reference_rsem               rsem
# normalize_with_rsem
# extract_rsem_signal
# 
# align_with_tophat                  Tophat
# extract_tophat_bamfolder
# summarize_tophat_alignment
# 
# index_reference_star               Star
# align_with_star
# extract_star_samfolder
# 
# count_with_htseq_count             HTSeq-Count
# summarize_htseq_count
# extract_htseq_count_signal
# convert_counts_to_cpm

from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
import SignalFile
import NGSQC as QC

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
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the reads in here have had adapters trimmed."),
    help="A folder that contains alignments from Tophat.",
    )

RSEMReferenceGenome = DataType(
    "RSEMReferenceGenome",
    #AttributeDef(
    #    "bowtie1_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    #AttributeDef(
    #    "bowtie2_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    #AttributeDef(
    #    "rsem_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    help="Indexed for rsem.",
    )

RSEMResults = DataType(
    "RSEMResults",
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the reads in here have had adapters trimmed."),
    AttributeDef(
        "align_to", ["genome", "transcriptome"], "genome", "genome",
        help="Align to the genome or transcriptome."),
    AttributeDef(
        "aligner", ["star", "bowtie1"], "star", "star"),
    #AttributeDef(
    #    "contents", BDT.CONTENTS, "unspecified", "unspecified"),
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
        "aligner", ["star", "star_unstranded"], "star", "star"),
    AttributeDef(
        "is_subset", YESNO, "no", "no",
        help="Generate a small sample file with a limited number of reads."
        ),
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the reads in here have had adapters trimmed."),
    help="Results from a STAR alignment.  Includes SAM files and other stuff.",
    )

HTSeqCountResults = DataType(
    "HTSeqCountResults",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef("aligner", NGS.RNA_ALIGNERS, "star", "star"),
    help="Results from HTSeq-Count.",
    )

HTSeqCountSummary = DataType(
    "HTSeqCountSummary",
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef("aligner", NGS.RNA_ALIGNERS, "star", "star"),
    help="Contains a summary of the results from htseq-count as a "
    "tab-delimited text file.",
    )

RNASeqUnprocessedSignalFile = DataType(
    "RNASeqUnprocessedSignalFile",

    AttributeDef(
        "format", ["unknown", "tdf", "pcl", "gct", "res", "jeffs"],
        "unknown", "tdf", help="file format"),
    AttributeDef(
        # Should limit these options to RNA-Seq specific ones.
        "preprocess", BDT.PREPROCESS, "unknown", "unknown",
        help="preprocess method"),
    AttributeDef(
        "logged", ["unknown", "no", "yes"], "unknown", "yes",
        help="logged or not"),
    #AttributeDef(
    #    "contents", BDT.CONTENTS, "unspecified", "unspecified",
    #    help="contents"),
    
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef("aligner", NGS.RNA_ALIGNERS, "star", "star"),

    AttributeDef(
        "expression_of", ["gene", "isoform"], "gene", "gene",
        help="Whether to get expression for overall genes or "
        "individual isoforms."),
    
    help="Contains attributes for an RNA-Seq analysis.",
    )
    

CompleteRNASeqAnalysis = DataType(
    "CompleteRNASeqAnalysis",
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef("aligner", NGS.RNA_ALIGNERS, "star", "star"),
    help="Contains the results from a standard RNA-Seq experiment.",
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
    
    RNASeqUnprocessedSignalFile,
    CompleteRNASeqAnalysis,
    ]


all_modules = [
    ModuleNode(
        "index_reference_rsem",
        [NGS.ReferenceGenome, NGS.GTFGeneModel], RSEMReferenceGenome,
        #Constraint("rsem_indexed", MUST_BE, "no", 0),
        #Consequence("rsem_indexed", SET_TO, "yes"),
        #Consequence("bowtie1_indexed", SET_TO, "yes"),
        #Consequence("bowtie2_indexed", SET_TO, "yes"),
        ),
    #ModuleNode(
    #    "is_rsemreference_rsem_indexed",
    #    RSEMReferenceGenome, RSEMReferenceGenome,
    #    Constraint("rsem_indexed", MUST_BE, "unknown"),
    #    Consequence("rsem_indexed", BASED_ON_DATA, ["no", "yes"]),
    #    ),
    #ModuleNode(
    #    "is_rsemreference_bowtie1_indexed",
    #    RSEMReferenceGenome, RSEMReferenceGenome,
    #    Constraint("bowtie1_indexed", MUST_BE, "unknown"),
    #    Consequence("bowtie1_indexed", BASED_ON_DATA, ["no", "yes"]),
    #    ),
    #ModuleNode(
    #    "is_rsemreference_bowtie2_indexed",
    #    RSEMReferenceGenome, RSEMReferenceGenome,
    #    Constraint("bowtie2_indexed", MUST_BE, "unknown"),
    #    Consequence("bowtie2_indexed", BASED_ON_DATA, ["no", "yes"]),
    #    ),
    
    ModuleNode(
        "normalize_with_rsem",
        [NGS.FastqFolder, NGS.SampleGroupFile, NGS.ReadStrandedness,
         RSEMReferenceGenome],
        RSEMResults,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Consequence("align_to", SET_TO_ONE_OF, ["genome", "transcriptome"]),
        # bowtie not implemented.
        Consequence("aligner", SET_TO, "star"),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        help="Use RSEM to estimate TPM or FPKM.",
        ),

    ModuleNode(
        "extract_rsem_signal", RSEMResults, RNASeqUnprocessedSignalFile,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        #OptionDef(
        #    "genes_or_isoforms", default="genes",
        #    help='Get the expression value for "genes" (DEFAULT) or '
        #    '"isoforms".'),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SET_TO_ONE_OF, ["tpm", "fpkm"]),
        Consequence("logged", SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        Constraint("aligner", MUST_BE, "star"),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("align_to", MUST_BE, "genome"),
        Consequence("expression_of", SET_TO_ONE_OF, ["gene", "isoform"]),
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
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 3),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        help="Align to a reference genome with tophat.",
        ),
    ModuleNode(
        "extract_tophat_bamfolder",
        TophatAlignmentFolder, NGS.BamFolder,

        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "tophat"),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        help="Pull the BAM files out of the Tophat result folders."
        ),
    ModuleNode(
        "summarize_tophat_alignment",
        TophatAlignmentFolder, TophatAlignmentSummary,
        help="Summarize the alignment, e.g. number of reads aligned.",
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
        OptionDef(
            "two_pass", default="yes", help="Whether to do 2-pass mapping.  "
            "Helpful for variant calling, novel junction discovery.  "
            '"yes" or "no"'),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        #Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 2),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "star"),
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
        "align_with_star_unstranded",
        [NGS.FastqFolder, NGS.SampleGroupFile, STARReferenceGenome],
        STARAlignmentFolder,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        #Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "star_unstranded"),
        
        help="Align to a reference genome with star, unstranded.  "
        "This is only used internally to determine strandedness of reads."
        ),
    ModuleNode(
        "extract_star_bamfolder",
        STARAlignmentFolder, NGS.BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("aligner", CAN_BE_ANY_OF, ["star", "star_unstranded"]),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        #Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
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
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("adapters_trimmed", SAME_AS, 0, 3),

        # Optimization.  Don't allow this.
        #Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        #Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 3),
        # Get error when using coordinate sorting:
        # Maximum alignment buffer size exceeded while pairing SAM alignments.
        #Constraint("sorted", CAN_BE_ANY_OF, ["name", "coordinate"], 0),
        Constraint("sorted", MUST_BE, "name", 0),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.RNA_ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
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
        "summarize_htseq_count",
        HTSeqCountResults, HTSeqCountSummary,
        Constraint("aligner", CAN_BE_ANY_OF, NGS.RNA_ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "extract_htseq_count_signal",
        HTSeqCountResults, RNASeqUnprocessedSignalFile,
        OptionDef(
            "ignore_htseq_count_errors", default="no",
            help='Whether to ignore errors in the files.  '
            'Should be "yes" or "no".',
            ),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.RNA_ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SET_TO, "counts"),
        Consequence("logged", SET_TO, "no"),
        Consequence("format", SET_TO, "tdf"),
        Consequence("expression_of", SET_TO, "gene"),
        ),
    ModuleNode(
        "convert_counts_to_cpm",
        RNASeqUnprocessedSignalFile, RNASeqUnprocessedSignalFile,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),

        Constraint("preprocess", MUST_BE, "counts"),
        Consequence("preprocess", SET_TO, "cpm"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("format", SET_TO, "tdf"),
        ),
    ModuleNode(
        "do_std_rna_seq_analysis",
        [NGS.BamFolder,                      # 0
         QC.FastQCSummary, QC.FastQCFolder,  # 1-2  no trimming
         QC.FastQCSummary, QC.FastQCFolder,  # 3-4  trimmed
         QC.RSeQCResults,                    # 5
         RNASeqUnprocessedSignalFile,        # 6    gene expression
         RNASeqUnprocessedSignalFile,        # 7    isoform expression
         NGS.NumAlignedReads,                # 8
         RNASeqUnprocessedSignalFile,        # 9
         HTSeqCountSummary,                  # 10
         ],
        CompleteRNASeqAnalysis,
        # BamFolder
        Constraint("aligner", CAN_BE_ANY_OF, NGS.RNA_ALIGNERS, 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        # FastQCSummary, FastQCFolder
        Constraint("adapters_trimmed", MUST_BE, "no", 1),
        Constraint("adapters_trimmed", MUST_BE, "no", 2),
        Constraint("adapters_trimmed", MUST_BE, "yes", 3),
        Constraint("adapters_trimmed", MUST_BE, "yes", 4),
        # RSeQCResults
        Constraint("aligner", SAME_AS, 0, 5),
        Constraint("adapters_trimmed", SAME_AS, 0, 5),
        # RNASeqUnprocessedSignalFile
        Constraint("aligner", SAME_AS, 0, 6),
        Constraint("adapters_trimmed", SAME_AS, 0, 6),
        Constraint("logged", MUST_BE, "no", 6),
        Constraint("preprocess", MUST_BE, "tpm", 6),
        Constraint("expression_of", MUST_BE, "gene", 6),
        # RNASeqUnprocessedSignalFile
        Constraint("aligner", SAME_AS, 0, 7),
        Constraint("adapters_trimmed", SAME_AS, 0, 7),
        Constraint("logged", MUST_BE, "no", 7),
        Constraint("preprocess", MUST_BE, "tpm", 7),
        Constraint("expression_of", MUST_BE, "isoform", 7),
        # NumAlignedReads
        Constraint("aligner", SAME_AS, 0, 8),
        Constraint("adapters_trimmed", SAME_AS, 0, 8),
        # RNASeqUnprocessedSignalFile
        Constraint("aligner", SAME_AS, 0, 9),
        Constraint("adapters_trimmed", SAME_AS, 0, 9),
        Constraint("logged", MUST_BE, "no", 9),
        Constraint("preprocess", MUST_BE, "counts", 9),
        # HTSeqCountSummary
        Constraint("aligner", SAME_AS, 0, 10),
        Constraint("adapters_trimmed", SAME_AS, 0, 10),

        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT, 0),
        ),

    ModuleNode(
        "rnasequnprocessedsignalfile_to_unprocessedsignalfile",
        RNASeqUnprocessedSignalFile, SignalFile.UnprocessedSignalFile,

        Constraint(
            "format", CAN_BE_ANY_OF,
            ["unknown", "tdf", "pcl", "gct", "res", "jeffs"]),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("logged", CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence("logged", SAME_AS_CONSTRAINT),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),

        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Constraint("aligner", CAN_BE_ANY_OF, ["star", "tophat"]),
        ),
    ]
