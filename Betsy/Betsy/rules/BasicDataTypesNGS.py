# Data Types:
# ReferenceGenome
# FullyIndexedReferenceGenome
# SampleGroupFile
# 
# FastqFolder
# SamFolder
# BamFolder
# SaiFolder                  # from BWA Backtrack
#
# RealignTargetFolder
# RecalibrationReport
#
# TophatAlignmentFolder
# Bowtie1AlignmentSummary
# Bowtie2AlignmentSummary
# TophatAlignmentSummary
# AlignedReadsSummary
#
# CoverageSummary
# TrimmomaticSummary
#
#
# Modules:
# is_fastq_folder_compressed
# uncompress_fastq_folder
# merge_reads
# check_single_or_paired_orientation
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
# align_with_tophat                Tophat
# extract_tophat_bamfolder
# is_reference_bwa_indexed         BWA
# index_reference_bwa              
# align_with_bwa_aln
# align_with_bwa_mem
# convert_sai_to_sam_folder
# 
# summarize_aligned_reads
# calculate_coverage
# 
# convert_sam_to_bam_folder
# index_bam_folder
# sort_bam_folder_by_coordinate
# sort_bam_folder_by_name
# sort_bam_folder_by_contig
# 
# add_read_groups_to_bam_folder
# mark_duplicates_bam_folder
# create_realign_targets
# realign_indels_bam_folder
# recalibrate_base_quality_score


# OBSOLETE:
# Bowtie1IndexedGenome       # Move to file for alignment.
# Bowtie2IndexedGenome
# BWAIndexedGenome


from Betsy.bie3 import *
import BasicDataTypes as BDT

ALIGNERS = [
    "unknown", "bowtie1", "bowtie2", "bwa_backtrack", "bwa_mem", "tophat"]
##REFERENCE_GENOMES = ["hg18", "hg19", "mm9", "dm3"]

COMPRESSION = ["unknown", "no", "gz", "bz2", "xz"]
COMPRESSION_NOT_UNKNOWN = [x for x in COMPRESSION if x != "unknown"]

SORT_ORDERS = ["no", "coordinate", "name", "contig"]
SORTED = [x for x in SORT_ORDERS if x != "no"]
# coordinate  samtools sort
# name        by read name samtools sort -n
# contig      match contig ordering of reference genome (Picard)

ORIENTATION = [
    "unknown", "single", "paired", "paired_fr", "paired_rf", "paired_ff"]
ORIENTATION_NOT_UNKNOWN = [x for x in ORIENTATION if x != "unknown"]


ReferenceGenome = DataType(
    "ReferenceGenome",
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
    help="Used only for indexing a new reference genome."
    )

SampleGroupFile = DataType(
    "SampleGroupFile",
    AttributeDef(
        "contents", BDT.CONTENTS,
        "unspecified", "unspecified", help="contents"),
    AttributeDef(
        "orientation", ORIENTATION, "unknown", "unknown",
        help="Either single-end reads, paired-end reads with orientation.  "
        "See the bowtie manual for a description of the orientation.",
        ),
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
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help=""),
    AttributeDef(
        "compressed", COMPRESSION, "unknown", "no",
        help="Whether the files are compressed (gz, bz2, xz)."),
    AttributeDef(
        "adapters_trimmed", ["yes", "no"], "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef(
        "reads_merged", ["yes", "no"], "no", "yes",
        help="Whether reads for a sample are merged into one file."),
    #AttributeDef(
    #    "orientation", ORIENTATION, "unknown", "unknown",
    #    help="Either single-end reads, paired-end reads with orientation.  "
    #    "See the bowtie manual for a description of the orientation.",
    #    ),
    help="A folder containing FASTQ files."
    )


SAM_ATTRIBUTES = [
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified"),
    AttributeDef(
        "aligner", ALIGNERS, "unknown", "bowtie2",
        help="Alignment algorithm."),

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
        "base_recalibrated", ["yes", "no"], "no", "no",
        help="recalibrated or not"),
    AttributeDef(
        "indel_realigned", ["yes", "no"], "no", "no",
        help="realigned or not"),
    #AttributeDef(
    #    "has_header", ["yes", "no"], "no", "yes",
    #    help="fix header or not"),
    AttributeDef(
        "has_read_groups", ["yes", "no"], "no", "no",
        help="Whether the file contains read groups.",
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
    )

RecalibrationReport = DataType(
    "RecalibrationReport",
    )

SaiFolder = DataType(
    "SaiFolder",
    AttributeDef(
        "contents", BDT.CONTENTS,
        "unspecified", "unspecified", help="contents"),
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

TophatAlignmentSummary = DataType(
    "TophatAlignmentSummary",
    help="Summarizes the alignment from tophat (.xls file).",
    )

AlignedReadsSummary = DataType(
    "AlignedReadsSummary",
    help="Summarizes the number of aligned reads (.xls file).",
    )

CoverageSummary = DataType(
    "CoverageSummary",
    help="Summarizes the coverage for an alignment.",
    )

TrimmomaticSummary = DataType(
    "TrimmomaticSummary",
    help="Summarizes the results from trimmomatic in an Excel .xls file.",
    )

TophatAlignmentFolder = DataType(
    "TophatAlignmentFolder",
    AttributeDef(
        "contents", BDT.CONTENTS,
        "unspecified", "unspecified", help="contents"),
    help="A folder that contains alignments from Tophat.",
    )

all_data_types = [
    ReferenceGenome,
    FullyIndexedReferenceGenome,
    SampleGroupFile,
    FastqFolder,
    SamFolder,
    BamFolder,
    SaiFolder,

    RealignTargetFolder,
    RecalibrationReport,
    
    Bowtie1AlignmentSummary,
    Bowtie2AlignmentSummary,
    TophatAlignmentFolder,
    TophatAlignmentSummary,
    
    AlignedReadsSummary,
    CoverageSummary,
    TrimmomaticSummary,
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
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT),

        #Constraint("compressed", CAN_BE_ANY_OF, ["no", "gz", "bz2", "xz"], 0),
        #Consequence("compressed", SET_TO, "no"),
        # Don't deal with compression here.
        Constraint("compressed", MUST_BE, "no", 0),
        Consequence("compressed", SAME_AS_CONSTRAINT),

        Constraint("adapters_trimmed", CAN_BE_ANY_OF, ["no", "yes"], 0),
        Consequence("adapters_trimmed", SAME_AS_CONSTRAINT),

        Constraint("reads_merged", MUST_BE, "no", 0),
        Consequence("reads_merged", SET_TO, "yes"),

        # Bug: why does this cause the RSEM pipeline to not work?
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION, 1),
        ),
    ModuleNode(
        "check_single_or_paired_orientation",
        [SampleGroupFile, FastqFolder, ReferenceGenome], SampleGroupFile,
        Constraint("orientation", MUST_BE, "unknown", 0),
        Consequence("orientation", BASED_ON_DATA, ORIENTATION_NOT_UNKNOWN),
        Constraint("compressed", MUST_BE, "no", 1),
        Constraint("reads_merged", MUST_BE, "yes", 1),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
        ),
    ModuleNode(
        "trim_adapters_trimmomatic",
        [FastqFolder, SampleGroupFile], FastqFolder,
        OptionDef(
            "adapters_fasta", 
            help="Full path to FASTA file containing adapters to trim.",
            ),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT),

        Constraint("compressed", MUST_BE, "no", 0),
        Consequence("compressed", SAME_AS_CONSTRAINT),
        Constraint("adapters_trimmed", MUST_BE, "no", 0),
        Consequence("adapters_trimmed", SET_TO, "yes"),
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
        Consequence("dict_added", BASED_ON_DATA, ["no", "yes"]),
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
        Consequence("samtools_indexed", BASED_ON_DATA, ["no", "yes"]),
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
        Consequence("bowtie1_indexed", BASED_ON_DATA, ["no", "yes"]),
        ),
    ModuleNode(
        "index_reference_bowtie1",
        ReferenceGenome, ReferenceGenome,
        #OptionDef(
        #    "assembly", default="genome",o
        #    help="Optional name for the genome assembly, e.g. hg19",
        #    ),
        Constraint("bowtie1_indexed", MUST_BE, "no"),
        Consequence("bowtie1_indexed", SET_TO, "yes"),
        ),
    ModuleNode(
        "align_with_bowtie1",
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        SamFolder,
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, ["no", "yes"], 0),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("bowtie1_indexed", MUST_BE, "yes", 2),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
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
        Consequence("bowtie2_indexed", BASED_ON_DATA, ["no", "yes"]),
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
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        SamFolder,
        #OptionDef(
        #    "orientation", default="fr",
        #    help="Which orientation.  See bowtie2 manual.  "
        #    "Values: fr, rf, ff.",
        #    ),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, ["no", "yes"], 0),
        #Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        Consequence("aligner", SET_TO, "bowtie2"),
        #Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
        help="Align to a reference genome with bowtie 2.",
        ),
    ModuleNode(
        "summarize_bowtie2_alignment",
        SamFolder, Bowtie2AlignmentSummary,
        Constraint("aligner", MUST_BE, "bowtie2"),
        help="Summarize the alignment, e.g. number of reads aligned.",
        # This rule isn't quite right.  Actually requires the "log"
        # files saved by bowtie, not the SAM files.  So this might
        # fail if the user doesn't provide the log files along with
        # the SAM files.
        ),
        
    ModuleNode(
        "align_with_tophat",
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        TophatAlignmentFolder,
        OptionDef(
            "tophat_gtf_file", default="",
            help="GTF file containing the gene information.",
            ),
        OptionDef(
            "tophat_transcriptome_fa", default="",
            help="FASTA file for Transcriptome, e.g. "
            "hg19.knownGene.fa.  Must be indexed by tophat."
            "Either this or tophat_gtf_file must be given."
            ),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, ["no", "yes"], 0),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("bowtie2_indexed", MUST_BE, "yes", 2),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        
        help="Align to a reference genome with tophat.",
        ),
    ModuleNode(
        "summarize_tophat_alignment",
        TophatAlignmentFolder, TophatAlignmentSummary,
        help="Summarize the alignment, e.g. number of reads aligned.",
        ),
    ModuleNode(
        "extract_tophat_bamfolder",
        TophatAlignmentFolder, BamFolder,

        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "tophat"),
        help="Pull the BAM files out of the Tophat result folders."
        ),
    ModuleNode(
        "is_reference_bwa_indexed",
        ReferenceGenome, ReferenceGenome,
        Constraint("bwa_indexed", MUST_BE, "unknown"),
        Consequence("bwa_indexed", BASED_ON_DATA, ["no", "yes"]),
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
        #[FastqFolder, SaiFile], SamFolder,
        [FastqFolder, SampleGroupFile, ReferenceGenome],
        SaiFolder,
        
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        #Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, ["no", "yes"], 0),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT),

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
        SamFolder,

        # bwa mem figures out the orientation itself.
        #Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        #Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", CAN_BE_ANY_OF, ["no", "yes"], 0),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT),

        Consequence("aligner", SET_TO, "bwa_mem"),
        help="bwa mem for reads >= 70 bp.  Now preferred",
        ),
    ModuleNode(
        "convert_sai_to_sam_folder",
        [FastqFolder, SaiFolder, ReferenceGenome, SampleGroupFile], SamFolder,

        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("adapters_trimmed", MUST_BE, "yes", 0),
        Constraint("bwa_indexed", MUST_BE, "yes", 2),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 3),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Constraint("contents", SAME_AS, 0, 3),
        Consequence("contents", SAME_AS_CONSTRAINT),
        
        Consequence("aligner", SET_TO, "bwa_backtrack"),
        help="Convert bwa's .sai alignments into .sam format.",
        ),

    ModuleNode(
        "summarize_aligned_reads",
        BamFolder, AlignedReadsSummary,
        Constraint("sorted", CAN_BE_ANY_OF, SORT_ORDERS),
        Constraint("indexed", MUST_BE, "yes"),
        help="Summarize the alignment, e.g. number of reads aligned.",
        ),

    ModuleNode(
        "calculate_coverage",
        [BamFolder, ReferenceGenome], CoverageSummary,
        OptionDef(
            "ignore_coverage_below", default="",
            help="If given, will ignore all regions of the genome with a "
            "coverage below this value (e.g. 1) when calculating mean "
            "coverage.  Provides better estimate of WES.",
            ),
        OptionDef(
            "features_bed", default="",
            help="A bed file for the regions of the genome to include.  "
            "e.g. for looking at coverage of exons.",
            ),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        help="Calculate the coverage for an alignment.",
        ),
    ModuleNode(
        "convert_sam_to_bam_folder",
        SamFolder, BamFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
      
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Consequence("has_read_groups", SET_TO, "no"),
        Consequence("sorted", SET_TO, "no"),
        Consequence("duplicates_marked", SET_TO, "no"),
        Consequence("base_recalibrated", SET_TO, "no"),
        Consequence("indel_realigned", SET_TO, "no"),
        #help="Convert SAM to BAM files.",
        ),
    ModuleNode(
        "index_bam_folder",
        BamFolder, BamFolder,
        Constraint("indexed", MUST_BE, "no"),
        Constraint("sorted", CAN_BE_ANY_OF, SORTED),
        Consequence("indexed", SET_TO, "yes"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        help="index bam folder",
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
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Consequence("contents", SAME_AS_CONSTRAINT),
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
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", CAN_BE_ANY_OF, ["no", "coordinate"]),
        Consequence("sorted", SET_TO, "name"),
        Constraint("indexed", MUST_BE, "no"),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "sort_bam_folder_by_contig",
        [BamFolder, ReferenceGenome], BamFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
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
        Consequence("base_recalibrated", SET_TO, "no"),
        help="mark duplicates in SamFile"
        ),
    ModuleNode(
        "create_realign_targets",
        [BamFolder, ReferenceGenome], RealignTargetFolder,
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
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        help="Find the intervals to target for realignment.",
        ),
    ModuleNode(
        "realign_indels_bam_folder",
        [BamFolder, ReferenceGenome, RealignTargetFolder], BamFolder,
        DefaultAttributesFrom(0),
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
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("indel_realigned", MUST_BE, "no", 0),
        Consequence("indel_realigned", SET_TO, "yes"),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("indexed", SET_TO, "no"),
        Consequence("base_recalibrated", SET_TO, "no"),
        help="Realign indels."
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
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        help="Calculate the statistics for base recalibration.",
        ),
    ModuleNode(
        "recalibrate_base_quality_score",
        [BamFolder, ReferenceGenome, RecalibrationReport], BamFolder,
        
        Constraint("base_recalibrated", MUST_BE, "no", 0),
        Consequence("base_recalibrated", SET_TO, "yes"),
        #Constraint("sorted", MUST_BE, "contig", 0),
        #Consequence("sorted", SAME_AS_CONSTRAINT),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Consequence("contents", SAME_AS_CONSTRAINT),
        
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Consequence("has_read_groups", SAME_AS_CONSTRAINT),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Consequence("indel_realigned", SAME_AS_CONSTRAINT),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        
        Consequence("indexed", SET_TO, "no"),
        ),
    ]
