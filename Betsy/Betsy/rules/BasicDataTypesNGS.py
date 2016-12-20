# TODO: clean up this file.
# Data Types:
# ReferenceGenome
# FullyIndexedReferenceGenome
# GTFGeneModel
#
# FastqFolder
# FastaFolder
# SamFolder
# BamFolder
# SaiFolder                  # from BWA Backtrack
#
# ReadOrientation            # Contains orientation of reads
# ReadStrandedness           # Contains strandedness of reads
#
# SampleGroupFile
#
# RealignTargetFolder
# RecalibrationReport
#
# TrimmomaticSummary
# NumAlignedReads
# AlignmentSummaryMetrics
# Bowtie1AlignmentSummary
# Bowtie2AlignmentSummary
# PerfectAlignmentSummary  # How many alignments have no mismatches
# AlignmentCIGARFolder     # Summarizes each alignment.  CIGAR, edit distance
# DepthOfCoverage          # Average sequencing depth.
# InsertSizeMetrics
# ReadDuplicationRate
#
# Modules:
# is_fastq_folder_compressed
# uncompress_fastq_folder
# merge_reads
# make_sample_fastq_folder
# trim_adapters_trimmomatic
# summarize_trimmomatic_trimming
#
# subtract_mouse_reads
# extract_fasta_folder
# 
# index_reference_complete
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
# count_aligned_reads
# summarize_alignment_cigar
# summarize_perfect_alignments
# calculate_coverage
# 
# infer_read_orientation
# infer_read_strandedness
#
# convert_sam_to_bam_folder
# index_bam_folder
# sort_bam_folder_by_coordinate
# sort_bam_folder_by_name
# sort_bam_folder_by_contig
#
# add_read_groups_to_bam_folder       GATK variant calling pipeline
# add_md_tags_to_bam_folder
# mark_duplicates_bam_folder
# split_n_trim_bam_folder
# create_realign_targets
# realign_indels_bam_folder
# make_base_recalibration_report
# recalibrate_base_quality_score
# 
# calculate_insert_size_metrics
# calculate_alignment_summary_metrics
# count_read_duplication_rate



from Betsy.bie3 import *
import BasicDataTypes as BDT

YESNO = BDT.YESNO  # for convenience

DNA_ALIGNERS = ["bowtie1", "bowtie2", "bwa_backtrack", "bwa_mem"]
RNA_ALIGNERS = ["tophat", "star", "star_unstranded"]
ALIGNERS = ["unknown"] + DNA_ALIGNERS + RNA_ALIGNERS

# star_unstranded is to figure out the strand (without causing cycles).
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


ReferenceGenome = DataType(
    "ReferenceGenome",
    # Should be AttributeDef.
    #OptionDef(
    #    "assembly", default="genome",
    #    help="Optional name for the genome assembly, e.g. hg19",
    #    ),

    # Process attributes in this order.
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


GTFGeneModel = DataType(
    "GTFGeneModel",
    #OptionDef(
    #    "assembly", default="genome",
    #    help="Optional name for the genome assembly, e.g. hg19",
    #    ),
    help="GTF file that contains the gene model.",
    )


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
        "mouse_reads_subtracted", YESNO, "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the reads in here have had adapters trimmed."),
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
    ]

BAM_ATTRIBUTES = SAM_ATTRIBUTES + [
    AttributeDef("indexed", YESNO, "no", "no"),
    AttributeDef("sorted", SORT_ORDERS, "no", "coordinate"),
    AttributeDef(
        "duplicates_marked", YESNO, "no", "no",
        help="mark duplicate or not"),
    AttributeDef(
        "duplicates_filtered", YESNO, "no", "no",
        help="Whether duplicates are filtered out of the BAM file."),
    AttributeDef(
        "split_n_trim", YESNO, "no", "no",
        help="Whether SplitNCigarReads has been applied."),
    AttributeDef(
        "base_quality_recalibrated", YESNO, "no", "no",
        help="base quality score recalibration"),
    AttributeDef(
        "indel_realigned", YESNO, "no", "no",
        help="realigned or not"),
    #AttributeDef(
    #    "has_header", YESNO, "no", "yes",
    #    help="fix header or not"),
    AttributeDef(
        "has_read_groups", YESNO, "no", "no",
        help="Whether the files contain read groups.",
        ),
    AttributeDef(
        "has_md_tags", YESNO, "no", "no",
        help="Whether the files have MD tags.",
        ),
    ]


SamFolder = DataType(
    "SamFolder",
    *SAM_ATTRIBUTES,
    help="A folder containing SAM files."
    )


BamFolder = DataType(
    "BamFolder",
    *BAM_ATTRIBUTES,
    help="A folder containing BAM files."
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
        "mouse_reads_subtracted", YESNO, "no", "no",
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
        "adapters_trimmed", ["yes", "no"], "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="Contains information about strandedness of reads (.json)."
    )


SampleGroupFile = DataType(
    "SampleGroupFile",
    #AttributeDef(
    #    "contents", BDT.CONTENTS,
    #    "unspecified", "unspecified", help="contents"),
    # Needed to prevent cycles.
    AttributeDef(
        "mouse_reads_subtracted", ["yes", "no"], "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="File that associates FastQ files with sample names and pairing "
    "information.  It can be either an Excel file or a tab-delimited text "
    "file with at least three columns.  These three columns should contain "
    'headers "Filename", "Sample", and "Pair".  '
    "The Filename column contains the name of the fastq file (without "
    "any path information.  The Sample column contains the name of each "
    "sample.  There may be duplicates in this column if multiple FastQ "
    "files contain data for the same sample.  The Pair column should be "
    'either "1" or "2" for mate pairs for paired-end sequencing, '
    "or blank for single end sequencing."
    )


RealignTargetFolder = DataType(
    "RealignTargetFolder",
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    AttributeDef(
        "duplicates_marked", YESNO, "no", "no",
        help="Whether the duplicates are marked."),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the reads in here have had adapters trimmed."),
    AttributeDef(
        "split_n_trim", YESNO, "no", "no",
        help="Whether SplitNCigarReads has been applied."),
    AttributeDef(
        "sorted", SORT_ORDERS, "no", "no",
        help="Whether SplitNCigarReads has been applied."),
    )


RecalibrationReport = DataType(
    "RecalibrationReport",
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    AttributeDef(
        "duplicates_marked", YESNO, "no", "no",
        help="Whether the duplicates are marked."),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the reads in here have had adapters trimmed."),
    AttributeDef(
        "split_n_trim", YESNO, "no", "no",
        help="Whether SplitNCigarReads has been applied."),
    AttributeDef(
        "sorted", SORT_ORDERS, "no", "no",
        help="Whether SplitNCigarReads has been applied."),
    )

TrimmomaticSummary = DataType(
    "TrimmomaticSummary",
    help="Summarizes the results from trimmomatic in an Excel .xls file.",
    )

NumAlignedReads = DataType(
    "NumAlignedReads",
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    help="Count the number of aligned reads (.xls file).",
    )

AlignmentSummaryMetrics = DataType(
    "AlignmentSummaryMetrics",
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    help="Picard CollectAlignmentSummaryMetrics.  Saves a directory of files.",
    )

Bowtie1AlignmentSummary = DataType(
    "Bowtie1AlignmentSummary",
    help="Summarizes the alignment from bowtie1 (.xls file).",
    )

Bowtie2AlignmentSummary = DataType(
    "Bowtie2AlignmentSummary",
    help="Summarizes the alignment from bowtie2 (.xls file).",
    )

PerfectAlignmentSummary = DataType(
    "PerfectAlignmentSummary",
    AttributeDef(
        "mouse_reads_subtracted", YESNO, "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    help="Summarizes the fraction of perfect alignments (.txt file).",
    )

AlignmentCIGARFolder = DataType(
    "AlignmentCIGARFolder",
    AttributeDef(
        "mouse_reads_subtracted", YESNO, "no", "no",
        help="For subtracting mouse reads from PDX models of FastqFolder"),
    AttributeDef(
        "adapters_trimmed", YESNO, "no", "no",
        help="Whether the adapters are trimmed."),
    help="For each alignment, show the CIGAR, MD, NM, and NH data "
    "(folder of .txt file).",
    )

DepthOfCoverage = DataType(
    "DepthOfCoverage",
    AttributeDef(
        "duplicates_filtered", YESNO, "no", "no",
        help="Whether duplicates are filtered out."),
    help="Count the average coverage.  Saves a directory of files.",
    )

InsertSizeMetrics = DataType(
    "InsertSizeMetrics",
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),
    #AttributeDef(
    #    "indel_realigned", ["yes", "no"], "no", "no",
    #    help="realigned or not"),
    help="Look at insert sizes of DNA-Seq.  Saves a directory of files.",
    )

ReadDuplicationRate = DataType(
    "ReadDuplicationRate",
    help="Summarizes the rate of duplicates in a BAM file.",
    )


all_data_types = [
    ReferenceGenome,
    FullyIndexedReferenceGenome,
    GTFGeneModel,
    
    FastqFolder,
    FastaFolder
    SamFolder,
    BamFolder,
    SaiFolder,
    
    ReadOrientation,
    ReadStrandedness,
    
    SampleGroupFile,

    RealignTargetFolder,
    RecalibrationReport,
    
    TrimmomaticSummary,
    NumAlignedReads,
    AlignmentSummaryMetrics,
    Bowtie1AlignmentSummary,
    Bowtie2AlignmentSummary,
    PerfectAlignmentSummary,
    AlignmentCIGARFolder,
    DepthOfCoverage,
    InsertSizeMetrics,
    ReadDuplicationRate,
    ]

all_modules = [
    ModuleNode(
        "is_fastq_folder_compressed",
        FastqFolder, FastqFolder,
        Constraint("compressed", MUST_BE, "unknown"),
        Consequence("compressed", BASED_ON_DATA, COMPRESSION_NOT_UNKNOWN),

        # Optimization: Don't check when is_subset=yes.  is_subset is
        # created by the system, and will never be in a compressed
        # state.
        Constraint("is_subset", MUST_BE, "no"),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "uncompress_fastq_folder",
        FastqFolder, FastqFolder,
        Constraint("compressed", CAN_BE_ANY_OF, ["gz", "bz2", "xz"]),
        Consequence("compressed", SET_TO, "no"),
        
        # Optimization: Don't uncompress when is_subset=yes.
        # is_subset is created by the system, and will never be in a
        # compressed state.
        Constraint("is_subset", MUST_BE, "no"),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "merge_reads",
        [FastqFolder, SampleGroupFile], FastqFolder,

        Constraint("reads_merged", MUST_BE, "no", 0),
        Consequence("reads_merged", SET_TO, "yes"),
        Constraint("compressed", CAN_BE_ANY_OF, ["no", "gz", "bz2", "xz"], 0),
        Consequence("compressed", SET_TO, "no"),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        # Bug: why does this cause the RSEM pipeline to not work?
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION, 1),

        # Optimization: Don't merge when is_subset=yes.  is_subset is
        # created by the system, and will never be in a compressed
        # state.
        Constraint("is_subset", MUST_BE, "no"),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "make_sample_fastq_folder",
        FastqFolder, FastqFolder,
        OptionDef(
            # Use 200k because that's what infer_experiments.py uses.
            # Actually, 1 million is better, because 200k on
            # RNA-Seq/bowtie2 may result in too few good reads.
            #"num_samples", default=1000000, help="Number of reads to sample.",
            "num_samples", default=200000, help="Number of reads to sample.",
            ),
        Constraint("is_subset", MUST_BE, "no"),
        Consequence("is_subset", SET_TO, "yes"),
        Constraint("compressed", MUST_BE, "no"),
        Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("reads_merged", MUST_BE, "yes"),
        Consequence("reads_merged", SAME_AS_CONSTRAINT),
        help="Make a small sample Fastq folder.",
        ),
    ModuleNode(
        "trim_adapters_trimmomatic",
        [FastqFolder, SampleGroupFile], FastqFolder,
        OptionDef(
            "adapters_fasta", 
            help="Full path to FASTA file containing adapters to trim.",
            ),
        Constraint("adapters_trimmed", MUST_BE, "no", 0),
        Consequence("adapters_trimmed", SET_TO, "yes"),
        Constraint("compressed", MUST_BE, "no", 0),
        Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Consequence("reads_merged", SAME_AS_CONSTRAINT),
        help="Trim adapters.  Files should be merged.",
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
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "extract_fasta_folder",
        BamFolder, FastaFolder,
        Consequence("compressed", SET_TO, "no"),
        Consequence("reads_merged", SET_TO, "yes"),
        ),

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
        #Constraint("dict_added", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("dict_added", SAME_AS_CONSTRAINT),
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
        #Constraint("dict_added", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("dict_added", SAME_AS_CONSTRAINT),
        #Constraint("samtools_indexed", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("samtools_indexed", SAME_AS_CONSTRAINT),
        Constraint("bowtie1_indexed", MUST_BE, "unknown"),
        Consequence("bowtie1_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_bowtie1",
        ReferenceGenome, ReferenceGenome,
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
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("bowtie1_indexed", MUST_BE, "yes", 3),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
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
        #Constraint("dict_added", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("dict_added", SAME_AS_CONSTRAINT),
        #Constraint("samtools_indexed", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("samtools_indexed", SAME_AS_CONSTRAINT),
        #Constraint("bowtie1_indexed", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("bowtie1_indexed", SAME_AS_CONSTRAINT),
        Constraint("bowtie2_indexed", MUST_BE, "unknown"),
        Consequence("bowtie2_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_bowtie2",
        ReferenceGenome, ReferenceGenome,
        Constraint("bowtie2_indexed", MUST_BE, "no"),
        Consequence("bowtie2_indexed", SET_TO, "yes"),
        ),
    ModuleNode(
        "align_with_bowtie2",
        [FastqFolder, SampleGroupFile, ReadOrientation, ReferenceGenome],
        BamFolder,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 3),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 2),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "bowtie2"),
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
        #Constraint("dict_added", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("dict_added", SAME_AS_CONSTRAINT),
        #Constraint("samtools_indexed", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("samtools_indexed", SAME_AS_CONSTRAINT),
        #Constraint("bowtie1_indexed", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("bowtie1_indexed", SAME_AS_CONSTRAINT),
        #Constraint("bowtie2_indexed", CAN_BE_ANY_OF, ["yes", "no"]),
        #Consequence("bowtie2_indexed", SAME_AS_CONSTRAINT),
        Constraint("bwa_indexed", MUST_BE, "unknown"),
        Consequence("bwa_indexed", BASED_ON_DATA, YESNO),
        ),
    ModuleNode(
        "index_reference_bwa",
        ReferenceGenome, ReferenceGenome,
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
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        help="bwa aln for reads < 70 bp.  Generally not used anymore.",
        ),
    ModuleNode(
        "align_with_bwa_mem",
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        BamFolder,

        # bwa mem figures out the orientation itself.
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
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
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("bwa_indexed", MUST_BE, "yes", 4),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("is_subset", SAME_AS, 0, 1),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Constraint("mouse_reads_subtracted", SAME_AS, 0, 1),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "bwa_backtrack"),
        help="Convert bwa's .sai alignments into .sam format.",
        ),
    ModuleNode(
        "count_aligned_reads",
        [FastqFolder, SampleGroupFile, BamFolder], NumAlignedReads,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT, 2),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 2),
        Consequence("aligner", SAME_AS_CONSTRAINT, 2),
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
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
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
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
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
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_filtered", SAME_AS_CONSTRAINT),
        help="Summarize the distribution of coverage across the genome.",
        ),
    ModuleNode(
        "convert_sam_to_bam_folder",
        SamFolder, BamFolder,
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no"),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("is_subset", CAN_BE_ANY_OF, YESNO),
        Consequence("is_subset", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
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
        Consequence("indexed", SET_TO, "no"),
        ),
    ModuleNode(
        "mark_duplicates_bam_folder",
        BamFolder, BamFolder,
        Constraint("duplicates_marked", MUST_BE, "no"),
        Consequence("duplicates_marked", SET_TO, "yes"),
        Constraint("duplicates_filtered", MUST_BE, "no"),
        Consequence("duplicates_filtered", SAME_AS_CONSTRAINT),
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
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_filtered", SAME_AS_CONSTRAINT),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("sorted", MUST_BE, "coordinate"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("split_n_trim", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("split_n_trim", SAME_AS_CONSTRAINT),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
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
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("duplicates_marked", SAME_AS, 0, 2),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_filtered", SAME_AS_CONSTRAINT),
        Constraint("indel_realigned", MUST_BE, "no", 0),
        Consequence("indel_realigned", SET_TO, "yes"),
        # create_realign_targets requires sort order to be "coordinate"
        Constraint("sorted", MUST_BE, "coordinate", 2),
        Constraint("sorted", SAME_AS, 2, 0),
        Consequence("sorted", SAME_AS_CONSTRAINT, 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Constraint("aligner", SAME_AS, 0, 2),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Consequence("indexed", SET_TO, "no"),
        Consequence("base_quality_recalibrated", SET_TO, "no"),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("split_n_trim", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("split_n_trim", SAME_AS, 0, 2),
        Consequence("split_n_trim", SAME_AS_CONSTRAINT),
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
        # BamFolder attributes should match those in
        # recalibrate_base_quality_scores.  Otherwise, the BamFolder
        # used for mkaing the recalibration report may be different
        # from the one whose scores are recalibrated.
        
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_filtered", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("split_n_trim", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("split_n_trim", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        help="Calculate the statistics for base recalibration "
        "(BaseRecalibrator).",
        ),
    ModuleNode(
        "recalibrate_base_quality_score",
        [BamFolder, ReferenceGenome, RecalibrationReport], BamFolder,
        
        Constraint("base_quality_recalibrated", MUST_BE, "no", 0),
        Consequence("base_quality_recalibrated", SET_TO, "yes"),
        
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("adapters_trimmed", SAME_AS, 0, 2),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        
        Constraint("split_n_trim", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("split_n_trim", SAME_AS, 0, 2),
        Consequence("split_n_trim", SAME_AS_CONSTRAINT),

        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("sorted", SAME_AS, 0, 2),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Consequence("has_read_groups", SAME_AS_CONSTRAINT),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Consequence("indel_realigned", SAME_AS_CONSTRAINT),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("duplicates_marked", SAME_AS, 0, 2),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("duplicates_filtered", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS, 0),
        Constraint("aligner", SAME_AS, 0, 2),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        
        Consequence("indexed", SET_TO, "no"),
        help="PrintReads",
        ),

    ModuleNode(
        "infer_read_orientation",
        [FastqFolder, SampleGroupFile, ReferenceGenome], ReadOrientation,
        Constraint("is_subset", MUST_BE, "yes", 0),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "infer_read_strandedness",
        [BamFolder, GTFGeneModel], ReadStrandedness,
        # ReadStrandedness is only used for RNA aligners.  So use an
        # RNA-Seq aligner.  Have to use a special one that doesn't
        # require ReadStrandedness, or it will cause a cycle.
        Constraint("aligner", MUST_BE, "star_unstranded", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),
        Constraint("mouse_reads_subtracted", MUST_BE, "no", 0),
        Consequence("mouse_reads_subtracted", SAME_AS_CONSTRAINT),
        Constraint("is_subset", MUST_BE, "yes", 0),
        ),

    ModuleNode(
        "calculate_insert_size_metrics",
        BamFolder, InsertSizeMetrics,
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("indexed", MUST_BE, "yes"),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO),
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO),
        #Constraint("indel_realigned", CAN_BE_ANY_OF, YESNO),
        #Consequence("indel_realigned", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        ),
    
    ModuleNode(
        "calculate_alignment_summary_metrics",
        [BamFolder, ReferenceGenome], AlignmentSummaryMetrics,
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("indexed", MUST_BE, "yes"),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO),
        Constraint("duplicates_filtered", CAN_BE_ANY_OF, YESNO),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO),
        #Constraint("indel_realigned", CAN_BE_ANY_OF, YESNO),
        #Consequence("indel_realigned", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "count_read_duplication_rate",
        BamFolder, ReadDuplicationRate,
        Constraint("indexed", CAN_BE_ANY_OF, YESNO),
        Constraint("sorted", CAN_BE_ANY_OF, SORT_ORDERS),
        Constraint("duplicates_marked", MUST_BE, "yes"),
        Constraint("split_n_trim", CAN_BE_ANY_OF, YESNO),
        Constraint("base_quality_recalibrated", CAN_BE_ANY_OF, YESNO),
        Constraint("indel_realigned", CAN_BE_ANY_OF, YESNO),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO),
        Constraint("has_md_tags", CAN_BE_ANY_OF, YESNO),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        ),
    ]
