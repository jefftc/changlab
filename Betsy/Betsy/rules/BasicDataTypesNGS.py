# Data Types:
# ReferenceGenome
# FullyIndexedReferenceGenome
# SampleGroupFile
# GTFGeneModel
# ReadOrientation            # Contains orientation of reads
# ReadStrandedness           # Contains strandedness of reads
# 
# FastqFolder
# FastaFolder
# SamFolder
# BamFolder
# SaiFolder                  # from BWA Backtrack
#
# RealignTargetFolder
# RecalibrationReport
#
# Bowtie1AlignmentSummary
# Bowtie2AlignmentSummary
# TrimmomaticSummary
# 
# NumAlignedReads
# DepthOfCoverage          # Average sequencing depth.
# AlignmentCIGARFolder     # Summarizes each alignment.  CIGAR, edit distance
# PerfectAlignmentSummary  # How many alignments have no mismatches
#
#
# Modules:
# is_fastq_folder_compressed
# uncompress_fastq_folder
# merge_reads
# trim_adapters_trimmomatic
# summarize_trimmomatic_trimming
#
# is_reference_dict_added          dictionary
# add_dict_to_reference
# is_reference_samtools_indexed    Samtools
# index_reference_samtools
# is_reference_bowtie1_indexed     Bowtie1
# index_reference_bowtie1
# align_with_bowtie1
# summarize_bowtie1_alignment
# is_reference_bowtie2_indexed     Bowtie2
# index_reference_bowtie2
# align_with_bowtie2
# summarize_bowtie2_alignment
# is_reference_bwa_indexed         BWA
# index_reference_bwa              
# align_with_bwa_aln
# align_with_bwa_mem
# convert_sai_to_sam_folder
#
# infer_read_orientation
# infer_read_strandedness
# 
# summarize_aligned_reads
# calculate_coverage_old
# 
# convert_sam_to_bam_folder
# index_bam_folder
# sort_bam_folder_by_coordinate
# sort_bam_folder_by_name
# sort_bam_folder_by_contig
# 
# add_read_groups_to_bam_folder       GATK variant calling pipeline
# mark_duplicates_bam_folder
# split_n_trim_bam_folder
# create_realign_targets
# realign_indels_bam_folder
# recalibrate_base_quality_score


# OBSOLETE:
# Bowtie1IndexedGenome       # Move to file for alignment.
# Bowtie2IndexedGenome
# BWAIndexedGenome


from Betsy.bie3 import *
import BasicDataTypes as BDT

YESNO = BDT.YESNO  # for convenience

ALIGNERS = [
    "unknown", "bowtie1", "bowtie2", "bwa_backtrack", "bwa_mem", "tophat",
    "star"]
##REFERENCE_GENOMES = ["hg18", "hg19", "mm9", "dm3"]

COMPRESSION = ["unknown", "no", "gz", "bz2", "xz"]
COMPRESSION_NOT_UNKNOWN = [x for x in COMPRESSION if x != "unknown"]

SORT_ORDERS = ["no", "coordinate", "name", "contig"]
SORTED = [x for x in SORT_ORDERS if x != "no"]
# coordinate  samtools sort
# name        by read name samtools sort -n
# contig      match contig ordering of reference genome (Picard)

## ORIENTATION = [
##     "unknown", "single", "paired", "paired_fr", "paired_rf", "paired_ff",
##     ]
## ORIENTATION_NOT_UNKNOWN = [x for x in ORIENTATION if x != "unknown"]

## STRANDED = [
##     "unknown", "unstranded", "firststrand", "secondstrand",
##     ]
## STRANDED_NOT_UNKNOWN = [x for x in STRANDED if x != "unknown"]


GTFGeneModel = DataType(
    "GTFGeneModel",
    #OptionDef(
    #    "assembly", default="genome",
    #    help="Optional name for the genome assembly, e.g. hg19",
    #    ),
    help="GTF file that contains the gene model.",
    )


ReadOrientation = DataType(
    # How read pairs oriented towards each other.
    # orientation =
    #  single
    #  paired_fr  FR  -> <-
    #  paired_rf  RF  <- ->
    #  paired_ff  FF  -> ->
    "ReadOrientation",
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="Contains information about orientation and strandedness of reads "
    "(.json)."
    )

ReadStrandedness = DataType(
    # Whether the first read (/1) aligns to the 5' end of the gene.
    # stranded = 
    #   unstranded     Can go on either strand.
    #   firststrand    First read (/1) matches revcomp (stranded=reverse)
    #   secondstrand   First read matches the transcript.

    "ReadStrandedness",
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="Contains information about strandedness of reads (.json)."
    )


ReferenceGenome = DataType(
    "ReferenceGenome",
    # Should be AttributeDef.
    #OptionDef(
    #    "assembly", default="genome",
    #    help="Optional name for the genome assembly, e.g. hg19",
    #    ),
    AttributeDef(
        "dict_added", ["unknown", "no", "yes"], "unknown", "unknown"),
    AttributeDef(
        "samtools_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    AttributeDef(
        "bowtie1_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    AttributeDef(
        "bowtie2_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    AttributeDef(
        "bwa_indexed", ["unknown", "no", "yes"], "unknown", "unknown"),
    help="Should be FASTA file with reference genome.",
    )

FullyIndexedReferenceGenome = DataType(
    "FullyIndexedReferenceGenome",
    help="Used only for indexing a new reference genome.",
    )

SampleGroupFile = DataType(
    "SampleGroupFile",
    #AttributeDef(
    #    "contents", BDT.CONTENTS,
    #    "unspecified", "unspecified", help="contents"),
    # Needed to prevent cycles.
    #AttributeDef(
    #    "mouse_reads_subtracted", ["yes", "no"], "no", "no",
    #    help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="File contains sample group infomation"
    )

## FastqFile = DataType(
##     "FastqFile",
##     AttributeDef(
##         "read", ["single", "pair", "pair1", "pair2"],
##         "single", "single", help="single or pair read"),
##     #AttributeDef(
##     #    "ref", REFERENCE_GENOMES, "hg19", "hg19",
##     #    help="ref species"),
##     AttributeDef(
##         "contents", BDT.CONTENTS,
##         "unspecified", "unspecified", help="contents"),
##     )

FastqFolder = DataType(
    "FastqFolder",
    #AttributeDef(
    #    "contents", BDT.CONTENTS, "unspecified", "unspecified",
    #    help=""),
    AttributeDef(
        "compressed", COMPRESSION, "unknown", "no",
        help="Whether the files are compressed (gz, bz2, xz)."),
    AttributeDef(
        "reads_merged", ["yes", "no"], "no", "yes",
        help="Whether reads for a sample are merged into one file."),
    AttributeDef(
        "adapters_trimmed", ["yes", "no"], "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    #AttributeDef(
    #    "orientation", ORIENTATION, "unknown", "unknown",
    #    help="Either single-end reads, paired-end reads with orientation.  "
    #    "See the bowtie manual for a description of the orientation.",
    #    ),
    AttributeDef(
        "is_subset", YESNO, "no", "no",
        help="Generate a small sample file with a limited number of reads."
        ),
    help="A folder containing FASTQ files."
    )


FastaFolder = DataType(
    "FastaFolder",
    AttributeDef(
        "compressed", COMPRESSION, "unknown", "no",
        help="Whether the files are compressed (gz, bz2, xz)."),
    AttributeDef(
        "reads_merged", ["yes", "no"], "no", "yes",
        help="Whether reads for a sample are merged into one file."),
    #AttributeDef(
    #    "adapters_trimmed", ["yes", "no"], "no", "no",
    #    help="Whether the adapters are trimmed."),
    help="A folder containing FASTA files."
    )


SAM_ATTRIBUTES = [
    #AttributeDef(
    #    "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    #AttributeDef(
    #    "orientation", ORIENTATION, "unknown", "unknown",
    #    help="Either single-end reads, paired-end reads with orientation.  "
    #    "See the bowtie manual.",
    #    ),
    #AttributeDef(
    #    "stranded", STRANDED, "unknown", "unknown",
    #    help="Whether these reads are stranded.  See the tophat manual.",
    #    ),
    AttributeDef(
        "is_subset", YESNO, "no", "no",
        help="Generate a small sample file with a limited number of reads."
        ),

    #AttributeDef(
    #    "ref", REFERENCE_GENOMES, "hg19", "hg19",
    #    help="ref species"),
    #AttributeDef(
    #    "sample_type", ["RNA", "DNA"],
    #    "RNA", "RNA", help="RNA or DNA type"),
    ]

BAM_ATTRIBUTES = SAM_ATTRIBUTES + [
    AttributeDef(
        "indexed", ["yes", "no"], "no", "no",
        ),
    AttributeDef("sorted", SORT_ORDERS, "no", "coordinate"),
    AttributeDef(
        "duplicates_marked", ["yes", "no"], "no", "no",
        help="mark duplicate or not"),
    AttributeDef(
        "split_n_trim", ["yes", "no"], "no", "no",
        help="Whether SplitNCigarReads has been applied."),
    AttributeDef(
        "base_quality_recalibrated", ["yes", "no"], "no", "no",
        help="base quality score recalibration"),
    AttributeDef(
        "indel_realigned", ["yes", "no"], "no", "no",
        help="realigned or not"),
    #AttributeDef(
    #    "has_header", ["yes", "no"], "no", "yes",
    #    help="fix header or not"),
    AttributeDef(
        "has_read_groups", ["yes", "no"], "no", "no",
        help="Whether the files contain read groups.",
        ),
    AttributeDef(
        "has_md_tags", ["yes", "no"], "no", "no",
        help="Whether the files have MD tags.",
        ),
    ]


## SamFile = DataType("SamFile", *SAM_ATTRIBUTES)
## BamFile = DataType("BamFile", *BAM_ATTRIBUTES)

SamFolder = DataType(
    "SamFolder", *SAM_ATTRIBUTES, **{"help":"A folder containing SAM files."})

BamFolder = DataType(
    "BamFolder", *BAM_ATTRIBUTES, **{"help":"A folder containing BAM files."})

RealignTargetFolder = DataType(
    "RealignTargetFolder",
    AttributeDef(
        "duplicates_marked", ["yes", "no"], "no", "no",
        help="Whether the duplicates are marked."),
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    )

RecalibrationReport = DataType(
    "RecalibrationReport",
    AttributeDef(
        "duplicates_marked", ["yes", "no"], "no", "no",
        help="Whether the duplicates are marked."),
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    )

SaiFolder = DataType(
    "SaiFolder",
    #AttributeDef(
    #    "orientation", ORIENTATION, "unknown", "unknown",
    #    help="Either single-end reads, paired-end reads with orientation.  "
    #    "See the bowtie manual.",
    #    ),
    AttributeDef(
        "is_subset", YESNO, "no", "no",
        help="Generate a small sample file with a limited number of reads."
        ),
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    #AttributeDef(
    #    "contents", BDT.CONTENTS,
    #    "unspecified", "unspecified", help="contents"),
    #AttributeDef(
    #    "orientation", ORIENTATION, "unknown", "unknown",
    #    help="Either single-end reads, paired-end reads with orientation.  "
    #    "See the bowtie manual for a description of the orientation.",
    #    ),
    #AttributeDef(
    #    "read", ["single", "pair", "pair1", "pair2"],
    #    "single", "single", help="single or pair read"),
    #AttributeDef(
    #    "ref", REFERENCE_GENOMES,
    #    "hg19", "hg19", help="ref species"),
    help=".sai file generated by BWA aln."
    )

Bowtie1AlignmentSummary = DataType(
    "Bowtie1AlignmentSummary",
    help="Summarizes the alignment from bowtie1 (.xls file).",
    )

Bowtie2AlignmentSummary = DataType(
    "Bowtie2AlignmentSummary",
    help="Summarizes the alignment from bowtie2 (.xls file).",
    )

NumAlignedReads = DataType(
    "NumAlignedReads",
    help="Count the number of aligned reads (.xls file).",
    )

AlignmentCIGARFolder = DataType(
    "AlignmentCIGARFolder",
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="For each alignment, show the CIGAR, MD, NM, and NH data "
    "(folder of .txt file).",
    )

PerfectAlignmentSummary = DataType(
    "PerfectAlignmentSummary",
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="Summarizes the fraction of perfect alignments (.txt file).",
    )

DepthOfCoverage = DataType(
    "DepthOfCoverage",
    help="Count the average coverage.  Saves a directory of files.",
    )

TrimmomaticSummary = DataType(
    "TrimmomaticSummary",
    help="Summarizes the results from trimmomatic in an Excel .xls file.",
    )

all_data_types = [
    ReadOrientation,
    ReadStrandedness,
    ReferenceGenome,
    FullyIndexedReferenceGenome,
    SampleGroupFile,
    FastqFolder,
    SamFolder,
    BamFolder,
    SaiFolder,
    GTFGeneModel,

    RealignTargetFolder,
    RecalibrationReport,
    
    Bowtie1AlignmentSummary,
    Bowtie2AlignmentSummary,
    TrimmomaticSummary,
    
    NumAlignedReads,
    DepthOfCoverage,
    AlignmentCIGARFolder,
    PerfectAlignmentSummary,
    ]

all_modules = [
    ModuleNode(
        "is_fastq_folder_compressed",
        FastqFolder, FastqFolder,
        Constraint("compressed", MUST_BE, "unknown"),
        Consequence("compressed", BASED_ON_DATA, COMPRESSION_NOT_UNKNOWN),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "uncompress_fastq_folder",
        FastqFolder, FastqFolder,
        Constraint("compressed", CAN_BE_ANY_OF, ["gz", "bz2", "xz"]),
        Consequence("compressed", SET_TO, "no"),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "merge_reads",
        [FastqFolder, SampleGroupFile], FastqFolder,

        Constraint("reads_merged", MUST_BE, "no", 0),
        Consequence("reads_merged", SET_TO, "yes"),
        Constraint("compressed", CAN_BE_ANY_OF, ["no", "gz", "bz2", "xz"], 0),
        Consequence("compressed", SET_TO, "no"),
        #Constraint("compressed", MUST_BE, "no", 0),
        #Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        # Bug: why does this cause the RSEM pipeline to not work?
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION, 1),
        ),
    ## ModuleNode(
    ##     "check_orientation",
    ##     [FastqFolder, SampleGroupFile, ReferenceGenome], FastqFolder,
    ##     #[SampleGroupFile, BamFolder, ReferenceGenome], SampleGroupFile,
    ##     Constraint("orientation", MUST_BE, "unknown", 0),
    ##     Consequence("orientation", BASED_ON_DATA, ORIENTATION_NOT_UNKNOWN),
    ##     Constraint("compressed", MUST_BE, "no", 0),
    ##     Consequence("compressed", SAME_AS_CONSTRAINT),
    ##     Constraint("reads_merged", MUST_BE, "yes", 0),
    ##     Consequence("reads_merged", SAME_AS_CONSTRAINT),
    ##     Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
    ##     Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
    ##     Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
    ##     Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
    ##     ),
    ModuleNode(
        "trim_adapters_trimmomatic",
        [FastqFolder, SampleGroupFile], FastqFolder,
        OptionDef(
            "adapters_fasta", 
            help="Full path to FASTA file containing adapters to trim.",
            ),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", MUST_BE, "no", 0),
        Consequence("adapters_trimmed", SET_TO, "yes"),
        Constraint("compressed", MUST_BE, "no", 0),
        Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Consequence("reads_merged", SAME_AS_CONSTRAINT),
        help="Trim adapters.  Files should be merged.",
        ),
    ModuleNode(
        "make_sample_fastq_folder",
        FastqFolder, FastqFolder,
        OptionDef(
            # Use 200k because that's what infer_experiments.py uses.
            "num_samples", default=200000, help="Number of reads to sample.",
            ),
        Constraint("is_subset", MUST_BE, "no"),
        Consequence("is_subset", SET_TO, "yes"),
        Constraint("compressed", MUST_BE, "no"),
        Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("reads_merged", MUST_BE, "yes"),
        Consequence("reads_merged", SAME_AS_CONSTRAINT),
        help="Make a small sample Fastq folder.",
        ),
    ModuleNode(
        "summarize_trimmomatic_trimming",
        FastqFolder, TrimmomaticSummary,
        Constraint("compressed", MUST_BE, "no"),
        Constraint("adapters_trimmed", MUST_BE, "yes"),
        Constraint("reads_merged", MUST_BE, "yes"),
        help="Summarize the trimmomatic results.",
        # This rule isn't quite right.  Actually requires the "log"
        # files saved by trimmomatic, not the FASTQ files.  So this
        # might fail if the user doesn't provide the log files along
        # with the FASTQ files.
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
        "index_reference_complete",
        [ReferenceGenome, ReferenceGenome, ReferenceGenome, ReferenceGenome,
         ReferenceGenome],
        FullyIndexedReferenceGenome,
        Constraint("dict_added", MUST_BE, "yes", 0),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Constraint("bowtie1_indexed", MUST_BE, "yes", 2),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 3),
        Constraint("bwa_indexed", MUST_BE, "yes", 4),
        ## rsem index includes bowtie1 and bowtie2.
        #Constraint("rsem_indexed", MUST_BE, "yes", 3),
        help="Do all known indexing on a reference genome.",
        ),
    ModuleNode(
        "is_reference_dict_added",
        ReferenceGenome, ReferenceGenome,
        Constraint("dict_added", MUST_BE, "unknown"),
        Consequence("dict_added", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "add_dict_to_reference",
        ReferenceGenome, ReferenceGenome,
        Constraint("dict_added", MUST_BE, "no"),
        Consequence("dict_added", SET_TO, "yes"),
        help="CreateSequenceDictionary.jar",
        ),
    ModuleNode(
        "is_reference_samtools_indexed",
        ReferenceGenome, ReferenceGenome,
        Constraint("samtools_indexed", MUST_BE, "unknown"),
        Consequence("samtools_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_samtools",
        ReferenceGenome, ReferenceGenome,
        Constraint("samtools_indexed", MUST_BE, "no"),
        Consequence("samtools_indexed", SET_TO, "yes"),
        help="samtools faidx",
        ),

    ModuleNode(
        "is_reference_bowtie1_indexed",
        ReferenceGenome, ReferenceGenome,
        Constraint("bowtie1_indexed", MUST_BE, "unknown"),
        Consequence("bowtie1_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_bowtie1",
        ReferenceGenome, ReferenceGenome,
        #OptionDef(
        #    "assembly", default="genome",
        #    help="Optional name for the genome assembly, e.g. hg19",
        #    ),
        Constraint("bowtie1_indexed", MUST_BE, "no"),
        Consequence("bowtie1_indexed", SET_TO, "yes"),
        ),
    ModuleNode(
        "align_with_bowtie1",
        [FastqFolder, SampleGroupFile, ReadOrientation, ReferenceGenome],
        SamFolder,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("bowtie1_indexed", MUST_BE, "yes", 3),
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 0),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("is_subset", MUST_BE, "no", 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "bowtie1"),
        help="Align to a reference genome with bowtie.",
        ),
    ModuleNode(
        "summarize_bowtie1_alignment",
        SamFolder, Bowtie1AlignmentSummary,
        Constraint("aligner", MUST_BE, "bowtie1"),
        help="Summarize the alignment, e.g. number of reads aligned.",
        # This rule isn't quite right.  Actually requires the "log"
        # files saved by bowtie, not the SAM files.  So this might
        # fail if the user doesn't provide the log files along with
        # the SAM files.
        ),
    
    ModuleNode(
        "is_reference_bowtie2_indexed",
        ReferenceGenome, ReferenceGenome,
        Constraint("bowtie2_indexed", MUST_BE, "unknown"),
        Consequence("bowtie2_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_bowtie2",
        ReferenceGenome, ReferenceGenome,
        #OptionDef(
        #    "assembly", default="genome",
        #    help="Optional name for the genome assembly, e.g. hg19",
        #    ),
        Constraint("bowtie2_indexed", MUST_BE, "no"),
        Consequence("bowtie2_indexed", SET_TO, "yes"),
        ),
    ModuleNode(
        "align_with_bowtie2",
        [FastqFolder, SampleGroupFile, ReadOrientation, ReferenceGenome],
        BamFolder,
        #OptionDef(
        #    "orientation", default="fr",
        #    help="Which orientation.  See bowtie2 manual.  "
        #    "Values: fr, rf, ff.",
        #    ),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 3),
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 0),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("is_subset", MUST_BE, "no", 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "bowtie2"),
        #Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
        help="Align to a reference genome with bowtie 2.",
        ),
    ModuleNode(
        "summarize_bowtie2_alignment",
        BamFolder, Bowtie2AlignmentSummary,
        Constraint("aligner", MUST_BE, "bowtie2"),
        Constraint("sorted", MUST_BE, "no"),
        help="Summarize the alignment, e.g. number of reads aligned.",
        # This rule isn't quite right.  Actually requires the "log"
        # files saved by bowtie, not the SAM files.  So this might
        # fail if the user doesn't provide the log files along with
        # the SAM files.
        ),
        
    ModuleNode(
        "is_reference_bwa_indexed",
        ReferenceGenome, ReferenceGenome,
        Constraint("bwa_indexed", MUST_BE, "unknown"),
        Consequence("bwa_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_bwa",
        ReferenceGenome, ReferenceGenome,
        #OptionDef(
        #    "assembly", default="genome",
        #    help="Optional name for the genome assembly, e.g. hg19",
        #    ),
        Constraint("bwa_indexed", MUST_BE, "no"),
        Consequence("bwa_indexed", SET_TO, "yes"),
        ),
    # TODO: Choose bwa aln or bwa mem automatically.
    ModuleNode(
        "align_with_bwa_aln",
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        SaiFolder,
        
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        #Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 0),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("is_subset", MUST_BE, "no", 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),

        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Consequence("aligner", SET_TO, "bwa"),
        #Consequence("sorted", SET_TO, "no"),
        #Consequence("duplicates_marked", SET_TO, "no"),
        #Consequence("recalibration", SET_TO, "no"),
        #Consequence("has_header", SET_TO, "no"),
        #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19", "mm9", "dm3"]),
        #Consequence("ref", SAME_AS_CONSTRAINT),
        #help="generate algiment in SaiFile to SamFile"
        help="bwa aln for reads < 70 bp.  Generally not used anymore.",
        ),
    ModuleNode(
        "align_with_bwa_mem",
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        BamFolder,

        # bwa mem figures out the orientation itself.
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        #Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 0),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("is_subset", MUST_BE, "no", 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "bwa_mem"),
        help="bwa mem for reads >= 70 bp.  Now preferred",
        ),
    ModuleNode(
        "convert_sai_to_sam_folder",
        [FastqFolder, SaiFolder, ReadOrientation, SampleGroupFile,
         ReferenceGenome], SamFolder,

        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("bwa_indexed", MUST_BE, "yes", 4),
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 0),
        #Constraint("orientation", SAME_AS, 0, 1),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("is_subset", MUST_BE, "no", 0),
        Constraint("is_subset", SAME_AS, 0, 1),
        Consequence("is_subset", SAME_AS_CONSTRAINT),

        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        #Constraint("contents", SAME_AS, 0, 3),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "bwa_backtrack"),
        help="Convert bwa's .sai alignments into .sam format.",
        ),

    ModuleNode(
        "count_aligned_reads",
        [FastqFolder, SampleGroupFile, BamFolder], NumAlignedReads,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 2),
        Constraint("sorted", CAN_BE_ANY_OF, SORT_ORDERS, 2),
        Constraint("indexed", MUST_BE, "yes", 2),
        help="Summarize the alignment, e.g. number of reads aligned.",
        ),

    ModuleNode(
        "summarize_alignment_cigar",
        BamFolder, AlignmentCIGARFolder,
        Constraint("mouse_reads_subtracted", MUST_BE, "no"),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("has_md_tags", MUST_BE, "yes"),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        Constraint("sorted", CAN_BE_ANY_OF, SORT_ORDERS),
        Constraint("indexed", CAN_BE_ANY_OF, YESNO),
        help="Summarize the number of matches for each alignment.",
        ),

    ModuleNode(
        "summarize_perfect_alignments",
        [FastqFolder, SampleGroupFile, AlignmentCIGARFolder],
        PerfectAlignmentSummary,
        OptionDef(
            "num_mismatches", default=0, help="Remove all reads with <= "
            "this number of mismatches against the mouse genome.",
            ),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 2),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        help="Summarize the fraction of reads with perfect alignments.",
        ),

    ModuleNode(
        "calculate_coverage",
        [BamFolder, ReferenceGenome], DepthOfCoverage,
        OptionDef(
            "ignore_coverage_below", default="1",
            help="If given, will ignore all regions of the genome with a "
            "coverage below this value (e.g. 1) when calculating mean "
            "coverage.  Provides better estimate of WES.",
            ),
        OptionDef(
            "features_bed", default="",
            help="A BED6 file for the regions of the genome to include.  "
            "e.g. for looking at coverage of exons.  Should have:  "
            "chrom chromStart (0-based) chromEnd name score strand.",
            ),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        help="Summarize the distribution of coverage across the genome.",
        ),
    ModuleNode(
        "convert_sam_to_bam_folder",
        SamFolder, BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        
        Consequence("has_read_groups", SET_TO, "no"),
        Consequence("sorted", SET_TO, "no"),
        Consequence("duplicates_marked", SET_TO, "no"),
        Consequence("base_quality_recalibrated", SET_TO, "no"),
        Consequence("indel_realigned", SET_TO, "no"),
        Consequence("has_md_tags", SET_TO, "no"),
        #help="Convert SAM to BAM files.",
        ),
    ModuleNode(
        "index_bam_folder",
        BamFolder, BamFolder,
        Constraint("indexed", MUST_BE, "no"),
        Constraint("sorted", MUST_BE, "coordinate"),
        Consequence("indexed", SET_TO, "yes"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        ),
    # Sorting.
    # coordinate -> name -> contig.
    # Sorting by contig must be last, because RNA-SeQC needs it.
    # Don't allow sorting by contig -> coordinate.  This prevents a cycle.
    ModuleNode(
        "sort_bam_folder_by_coordinate",
        BamFolder, BamFolder,
        Constraint("indexed", MUST_BE, "no", 0),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        Constraint("sorted", CAN_BE_ANY_OF, ["no"], 0),
        Consequence("sorted", SET_TO, "coordinate"),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("duplicates_marked", MUST_BE, "no"),
        #Constraint("recalibrated", MUST_BE, "no"),
        #Constraint("has_header", MUST_BE, "no"),
        #Consequence("has_header", SAME_AS_CONSTRAINT),
        #Consequence("recalibrated", SAME_AS_CONSTRAINT),
        #Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "sort_bam_folder_by_name",
        BamFolder, BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", CAN_BE_ANY_OF, ["no", "coordinate"]),
        Consequence("sorted", SET_TO, "name"),
        Constraint("indexed", MUST_BE, "no"),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "sort_bam_folder_by_contig",
        [BamFolder, ReferenceGenome], BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", CAN_BE_ANY_OF, ["no", "name", "coordinate"]),
        Constraint("indexed", MUST_BE, "no"),
        Consequence("sorted", SET_TO, "contig"),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "add_read_groups_to_bam_folder",
        BamFolder, BamFolder,
        Constraint("has_read_groups", MUST_BE, "no"),
        Consequence("has_read_groups", SET_TO, "yes"),
        Constraint("sorted", MUST_BE, "coordinate"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("indexed", MUST_BE, "no"),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "add_md_tags_to_bam_folder",
        [BamFolder, ReferenceGenome], BamFolder,
        Constraint("has_md_tags", MUST_BE, "no"),
        Consequence("has_md_tags", SET_TO, "yes"),
        # Needs to be sorted or will run impossibly slow.
        Constraint("sorted", MUST_BE, "coordinate"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        #Constraint("indexed", MUST_BE, "yes"),
        Consequence("indexed", SET_TO, "no"),
        ),
    ModuleNode(
        "mark_duplicates_bam_folder",
        BamFolder, BamFolder,
        Constraint("duplicates_marked", MUST_BE, "no"),
        Consequence("duplicates_marked", SET_TO, "yes"),
        Constraint("sorted", MUST_BE, "coordinate"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("has_read_groups", MUST_BE, "yes"),
        Consequence("has_read_groups", SAME_AS_CONSTRAINT),
        Consequence("indexed", SET_TO, "no"),
        Consequence("indel_realigned", SET_TO, "no"),
        Consequence("base_quality_recalibrated", SET_TO, "no"),
        help="mark duplicates in SamFile"
        ),
    ModuleNode(
        "split_n_trim_bam_folder",
        [BamFolder, ReferenceGenome], BamFolder,
        Constraint("split_n_trim", MUST_BE, "no", 0),
        Consequence("split_n_trim", SET_TO, "yes"),
        Constraint("indexed", MUST_BE, "yes", 0),
        Consequence("indexed", SET_TO, "no"),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Consequence("sorted", SET_TO, "no"),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Consequence("has_read_groups", SAME_AS_CONSTRAINT),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("indel_realigned", SET_TO, "no"),
        Consequence("base_quality_recalibrated", SET_TO, "no"),
        help="Run SplitNCigarReads (for variant calling in RNA-Seq)",
        ),
    ModuleNode(
        "create_realign_targets",
        [BamFolder, ReferenceGenome], RealignTargetFolder,
        OptionDef(
            "filter_reads_with_N_cigar", default="no",
            help="Do for RNA-Seq data.",
            ),
        OptionDef(
            "realign_known_sites1", 
            help="For filtering out known SNPs.",
            ),
        OptionDef(
            "realign_known_sites2", default="",
            help="(OPTIONAL) Second database for filtering out known SNPs.",
            ),
        OptionDef(
            "realign_known_sites3", default="",
            help="(OPTIONAL).",
            ),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, ["yes", "no"], 0),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        #Constraint("base_recalibrated", MUST_BE, "no", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        help="Find the intervals to target for realignment "
        "(RealignerTargetCreator).",
        ),
    ModuleNode(
        "realign_indels_bam_folder",
        [BamFolder, ReferenceGenome, RealignTargetFolder], BamFolder,
        #DefaultAttributesFrom(0),
        OptionDef(
            "realign_known_sites1", 
            help="For filtering out known SNPs.",
            ),
        OptionDef(
            "realign_known_sites2", default="",
            help="(OPTIONAL) Second database for filtering out known SNPs.",
            ),
        OptionDef(
            "realign_known_sites3", default="",
            help="(OPTIONAL).",
            ),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Consequence("has_read_groups", SAME_AS_CONSTRAINT),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, ["yes", "no"], 0),
        Constraint("duplicates_marked", SAME_AS, 0, 2),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        #Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("indel_realigned", MUST_BE, "no", 0),
        Consequence("indel_realigned", SET_TO, "yes"),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Constraint("aligner", SAME_AS, 0, 2),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Consequence("indexed", SET_TO, "no"),
        Consequence("base_quality_recalibrated", SET_TO, "no"),
        help="Realign indels (IndelRealigner)."
        ),
    ModuleNode(
        "make_base_recalibration_report",
        [BamFolder, ReferenceGenome], RecalibrationReport,
        OptionDef(
            "recal_known_sites1", 
            help="For filtering out known SNPs.",
            ),
        OptionDef(
            "recal_known_sites2", default="",
            help="(OPTIONAL) Second database for filtering out known SNPs.",
            ),
        OptionDef(
            "recal_known_sites3", default="",
            help="(OPTIONAL).",
            ),
        Constraint("base_quality_recalibrated", MUST_BE, "no", 0),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        #Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, ["yes", "no"], 0),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Calculate the statistics for base recalibration "
        "(BaseRecalibrator).",
        ),
    ModuleNode(
        "recalibrate_base_quality_score",
        [BamFolder, ReferenceGenome, RecalibrationReport], BamFolder,
        
        Constraint("base_quality_recalibrated", MUST_BE, "no", 0),
        Consequence("base_quality_recalibrated", SET_TO, "yes"),
        #Constraint("sorted", MUST_BE, "contig", 0),
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Consequence("has_read_groups", SAME_AS_CONSTRAINT),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Consequence("indel_realigned", SAME_AS_CONSTRAINT),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, ["yes", "no"], 0),
        Constraint("duplicates_marked", SAME_AS, 0, 2),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Constraint("aligner", SAME_AS, 0, 2),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        
        Consequence("indexed", SET_TO, "no"),
        help="PrintReads",
        ),

    ModuleNode(
        "subtract_mouse_reads",
        [FastqFolder, SampleGroupFile, AlignmentCIGARFolder],
        FastqFolder,
        OptionDef(
            "num_mismatches", default=0, help="Remove all reads with <= "
            "this number of mismatches against the mouse genome.",
            ),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 2),
        Consequence("mouse_reads_subtracted", SET_TO, "yes"),
        Constraint("compressed", MUST_BE, "no", 0),
        Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Consequence("reads_merged", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "infer_read_orientation",
        [FastqFolder, SampleGroupFile, ReferenceGenome], ReadOrientation,
        Constraint("is_subset", MUST_BE, "yes", 0),
        Constraint("compressed", MUST_BE, "no", 0),
        #Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        #Consequence("reads_merged", SAME_AS_CONSTRAINT),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        #Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "infer_read_strandedness",
        [BamFolder, GTFGeneModel], ReadStrandedness,
        # Needs to be done with a non-RNA-Seq aligner, otherwise will
        # cause a cycle.
        #Constraint(
        #    "aligner", CAN_BE_ANY_OF,
        #    ["bowtie1", "bowtie2", "bwa_backtrack", "bwa_mem"], 0),
        Constraint("aligner", MUST_BE, "bowtie2", 0),
        Constraint("mouse_reads_subtracted", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("is_subset", MUST_BE, "yes", 0),
        ),

    ModuleNode(
        "extract_fasta_folder",
        BamFolder, FastaFolder,
        Consequence("compressed", SET_TO, "no"),
        Consequence("reads_merged", SET_TO, "yes"),
        ),

    ## ModuleNode(
    ##     "align_sample_fastq_folder",
    ##     [FastqFolder, SampleGroupFile, ReferenceGenome], SamFolder,
    ##     Constraint("compressed", MUST_BE, "no", 0),
    ##     Constraint("reads_merged", MUST_BE, "yes", 0),
    ##     Constraint("is_subset", MUST_BE, "yes", 0),
    ##     Consequence("is_subset", SAME_AS_CONSTRAINT),
    ##     Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
    ##     Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
    ##     Consequence("aligner", SET_TO, "bowtie2"),
    ##     help="Align to a reference genome with bowtie 2.",
    ##     ),
    
    ]
