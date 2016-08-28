# DataTypes:
# VCFFolder             # Folder of vcf files (one caller) of many samples.
# ManySampleVCFFile     # One VCF file with many samples.
# ManyCallerVCFFolders  # Multiple VCFFolders, one for each caller.
#
# SimpleVariantFile     # File with variants. many samples, many callers.
# SimpleVariantMatrix   # Summarizes variants in a matrix.
#
# Pipeline:
# 1.  ManyCallerVCFFolders               Calls for all variant callers.
# 2.  SimpleVariantFile.filtered=no      One file, merge from all folders.
# 3.  SimpleVariantFile.filtered=yes     Filter based on the FILTER VCF column.
# 4.  UnprocessedSimpleVariantMatrix     Compact Matrix format.
#
#     _SimpleVariantMatrix1    Annotation and filtering.
#     =====================
# 5.  SimpleVariantMatrix.annotated=yes  Annotated with Annovar.
#     annotate_simplevariantmatrix
#     No requirements.  Can be done at any time.
# 
# 6.  SimpleVariantMatrix.filtered=yes   Filter on reads, exons, etc.
#     filter_simplevariantmatrix_calls
#
#     _SimpleVariantMatrix2    Coverage
#     =====================
# 7.  SimpleVariantMatrix.with_coverage=yes
#     add_coverage_to_simplevariantmatrix
#
# 8.  SimpleVariantMatrix.with_rna_coverage=yes
#     add_rna_coverage_to_simplevariantmatrix
#     Must be: with_coverage="yes"
#
#     _SimpleVariantMatrix3    Adding external annotations.
#     =====================
# 8.  SimpleVariantMatrix.with_gxp=yes   Add gene expression from SignalFile.
#     add_gene_expression_to_simplevariantmatrix
# 
# 9.  SimpleVariantMatrix.with_cancer_genes=yes
#     add_cancer_genes_to_simplevariantmatrix
#
#     SimpleVariantMatrix      Final
#     ===================
#
#
# VCFRecalibrationReport            For GATK VariantRecalibrator.
# PileupSummary                     Used for VarScan calling.
# PositionsFile                     0-based coords
# NormalCancerFile
# IntervalListFile
# PositionSpecificDepthOfCoverage   Depth of coverage for specific position.
#                                   1-based coords
#
# NEED TO RENAME.  NOT VCF FOLDER
# AnnotatedVCFFolder      # No.  Not VCF files.  Need to clean this up.
# AnnotatedMultiVCFFile   # No.  Not VCF files.  Need to clean this up.
#
# FILTERS
# filter_simplevariantfile       SimpleVariantFile
#   remove_samples               Get rid of matched normal sample.
#   remove_radia_rna_samples     Get rid of RNA-Seq samples for Radia.
#   apply_filter                 Remove variants based on VCF FILTER column.
#   wgs_or_wes                   Needed for filtering MuSE calls.
# filter_simplevariantmatrix_calls   _SimpleVariantMatrix1
#   * Filters specific calls.
#   filter_by_min_alt_reads      At least this num of ALT reads. (No Strelka)
#   filter_by_min_total_reads
#   filter_by_min_vaf            VAF >= this value.
# filter_simplevariantmatrix_variants   SimpleVariantMatrix
#   * Filters variants.
#   min_callers_in_every_sample         Minimum callers in every sample.
#   min_coverage_in_every_sample
#   min_callers_in_any_sample           Minimum callers for >= 1 sample.
#   min_gene_expression_every_sample    Needs with_gxp=yes.
#   nonsynonymous_and_stopgain_only     Needs annotated=yes.
#   sift_polyphen_damaging
#
#
# Modules:
# summarize_variants_mpileup
# call_variants_mpileup
# call_variants_GATK
# call_variants_platypus
# call_consensus_varscan
# call_variants_varscan
# merge_variants_snp
#
# call_somatic_varscan
# call_variants_mutect
# call_variants_strelka
# call_variants_somaticsniper
# call_variants_jointsnvmix
# call_variants_muse
# merge_somatic_variants_snp
#
# make_vcf_recalibration_report_snp
# recalibrate_variants_snp
# filter_snps_only_multivcf
# annotate_with_annovar_vcffolder
# merge_vcf_folder
# annotate_multivcf_annovar
# extract_positions_from_multivcf_file
# backfill_vcf_folder


# Recalibrate variant scores with GATK.
# https://www.broadinstitute.org/gatk/guide/article?id=2805


#                 SNP   INDEL  NOTES
# gatk             Y      Y    Makes file with both.  Can filter out.
# mutect           Y      N
# varscan          Y      Y    Makes separate files.
# strelka          Y      Y    Makes separate files.
# somaticsniper    Y      N
# jointsnvmix      Y      Y    Makes file with both.  Can filter out.
# muse             Y      N
# radia            Y      Y    Makes file with both.  Can filter out.
# platypus         Y      Y    Makes file with both.  Can filter out.
        

# Calculate the consensus reads at specific positions.
# PositionsFile -> summarize_consensus_mpileup -> PileupSummary
#   -> call_consensus_varscan
# coordinates_from is either "unknown" or "simplevariantmatrix".
#
# Summarize the reads across the genome.
# summarize_reads_mpileup -> PileupSummary -> call_variants_varscan
#                                          -> call_somatic_varscan
# coordinate_from is "whole_genome".

    

from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS
import GeneExpProcessing as GXP

CALLERS = [
    "none", "mpileup", "gatk", "platypus", "varscan", "mutect", "strelka",
    "somaticsniper", "jointsnvmix", "muse", "radia"]
VARTYPES = ["all", "snp", "indel"]
#VARTYPE_NOT_CONSENSUS = [x for x in VARTYPES if x != "consensus"]
BACKFILLS = ["no", "yes", "consensus"]

COORDINATES_FROM = ["unknown", "simplevariantmatrix", "whole_genome"]

YESNO = BDT.YESNO  # for convenience



## BCFFolder = DataType(
##     "BCFFolder",
##     AttributeDef(
##         "contents", BDT.CONTENTS,
##         "unspecified", "unspecified", help="contents"),
##     AttributeDef(
##         "vartype", ["both", "snp", "indel"], "both", "snp",
##         help="What kind of variants are held in this file."),
##     AttributeDef(
##         "get_coverage", ["no", "yes"], "no", "no",
##         help="Whether the purpose is to get depth of coverage."),
##     help="Folder of .bcf files generated by samtools mpileup.",
##     )


VCFFolder = DataType(
    "VCFFolder",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "caller", CALLERS, "none", "mpileup",
        help="Which variant caller was used."),

    #AttributeDef(
    #    "mpileup_summary", ["no", "yes"], "no", "no",
    #    help="Whether this is just summary information from mpileup."),
    AttributeDef(
        "is_consensus", YESNO, "no", "no",
        help="Whether this is a consensus call."),
    AttributeDef(
        "vartype", VARTYPES, "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "vcf_recalibrated", YESNO, "no", "no",
        help="Whether quality scores are ready for filtering."),
    AttributeDef(
        "somatic", YESNO, "no", "no",
        help="Whether the variants here are from the somatic cancer genome "
        "(no germline)."),
    AttributeDef(
        "backfilled", BACKFILLS, "no", "no",
        help="Whether the mutations are backfilled."),
    AttributeDef(
        "coordinates_from", COORDINATES_FROM, "unknown", "unknown",
        help="Where do these coordinates come from."),
    AttributeDef(
        "aligner", NGS.ALIGNERS, "unknown", "unknown",
        help="What alignment algorithm used.  Helpful for determining "
        "whether this contains DNA or RNA data."),
    # For PileupSummary -> call_consensus_varscan ->
    # summarize_coverage_at_positions ->
    # PositionSpecificDepthOfCoverage
    AttributeDef(
        "filtered_calls", YESNO, "no", "no",
        help="Whether the calls are filtered."),
    )

VCFRecalibrationReport = DataType(
    "VCFRecalibrationReport",
    AttributeDef(
        "vartype", VARTYPES, "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "somatic", YESNO, "no", "no",
        help="Whether the variants here are from the somatic cancer genome "
        "(no germline)."),
    AttributeDef(
        "caller", CALLERS, "none", "mpileup",
        help="Which variant caller was used."),
    )

AnnotatedVCFFolder = DataType(
    "AnnotatedVCFFolder",
    )

ManySampleVCFFile = DataType(
    "ManySampleVCFFile",
    # Headers are:
    # #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [Samples...]
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "caller", CALLERS, "none", "mpileup",
        help="Which variant caller was used."),
    AttributeDef(
        "vartype", VARTYPES, "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "somatic", YESNO, "no", "no",
        help="Whether the variants here are from the somatic cancer genome "
        "(no germline)."),
    AttributeDef(
        "backfilled", ["no", "yes", "consensus"],
        # na         backfill is irrelevant for this pipeline
        # missing    not yet backfilled
        # filled     backfilled
        # consensus  contains information that can be used for backfill
        "no", "no",
        #"backfilled", YESNO, "no", "no",
        help="Whether the mutations are backfilled."),
    #AttributeDef(
    #    "backfilled", YESNO, "no", "no",
    #    help="Whether the mutations are backfilled."),
    )

AnnotatedMultiVCFFile = DataType(
    "AnnotatedMultiVCFFile",
    #AttributeDef(
    #    "backfilled", YESNO, "no", "no",
    #    help="Whether the mutations are backfilled."),
    )

PileupSummary = DataType(
    "PileupSummary",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "vartype", VARTYPES, "snp", "snp",
        help="What kind of variants are held in this file."),
    # Always consensus
    AttributeDef(
        "coordinates_from", COORDINATES_FROM, "unknown", "unknown",
        help="Where do these coordinates come from."),
    AttributeDef(
        "aligner", NGS.ALIGNERS, "unknown", "unknown",
        help="Which aligner used to generate this.  Helpful for "
        "distinguishing DNA and RNA data."),
    AttributeDef(
        "filtered_calls", YESNO, "no", "no",
        help="Whether the calls are filtered."),
    )

# For samtools.  Two columns (no header):
# <chrom>  <pos 0-based>  
PositionsFile = DataType(
    "PositionsFile",
    AttributeDef(
        "vartype", VARTYPES, "snp", "snp",
        help="What kind of variants are held in this file."),
    # Enables the modules to specify what positions to use.
    AttributeDef(
        "coordinates_from", COORDINATES_FROM, "unknown", "unknown",
        help="Where do these coordinates come from."),
    # Whether the positions are from a filtered or unfiltered
    # SimpleVariantMatrix2.
    AttributeDef(
        "filtered_calls", YESNO, "no", "no",
        help="Whether the calls are filtered."),
    )

NormalCancerFile = DataType(
    "NormalCancerFile",
    #AttributeDef(
    #    "mouse_reads_subtracted", ["yes", "no"], "no", "no",
    #    help="For subtracting mouse reads from PDX models of FastqFolder"),
    help="File contains correspondence between tumor and normal samples.  "
    "Should be a tab-delimited text file (or Excel file) containing two "
    'columns with headers "Normal" and "Cancer".  The "Normal" column '
    'contains the name of the normal sample.  The "Cancer" column '
    "contains the name of the cancer sample.  The corresponding BAM files "
    "should be named <sample>.bam."
    )

# To tell MuTect where to search.
IntervalListFile = DataType(
    "IntervalListFile",
    help="Genomic intervals in GATK format."
    )

PositionSpecificDepthOfCoverage = DataType(
    "PositionSpecificDepthOfCoverage",
    AttributeDef(
        "coordinates_from", COORDINATES_FROM, "unknown", "unknown",
        help="Where do these coordinates come from."),
    AttributeDef(
        "vartype", VARTYPES, "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "filtered_calls", YESNO, "no", "no",
        help="Whether the calls are filtered.  Requires annotated=yes"),
    AttributeDef(
        "aligner", NGS.ALIGNERS, "unknown", "unknown",
        help="What alignment algorithm used.  Helpful for determining "
        "whether this contains DNA or RNA data."),
    help="The amount of sequencing coverage at specific positions (0-based).",
    )


ManyCallerVCFFolders = DataType(
    "ManyCallerVCFFolders",
    AttributeDef(
        "vartype", ["snp", "indel"], "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "somatic", YESNO, "no", "no",
        help="Whether the variants here are from the somatic cancer genome "
        "(no germline)."),
    )

SimpleVariantFile = DataType(
    "SimpleVariantFile",
    AttributeDef(
        "vartype", ["snp", "indel"], "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "filtered", YESNO, "no", "no",
        help="Whether these variants are filtered."),
    help="1-based coordinates"
    )

UnprocessedSimpleVariantMatrix_ATTRIBUTES = [
    AttributeDef(
        "vartype", ["snp", "indel"], "snp", "snp",
        help="What kind of variants are held in this file."),
    ]

_SimpleVariantMatrix1_ATTRIBUTES = \
    UnprocessedSimpleVariantMatrix_ATTRIBUTES + [
    AttributeDef(
        "annotated", YESNO, "no", "no",
        help="Whether these variants are annotated with Annovar."),
    AttributeDef(
        "filtered_calls", YESNO, "no", "no",
        help="Whether the calls are filtered.  Requires annotated=yes"),
    ]

_SimpleVariantMatrix2_ATTRIBUTES = \
    _SimpleVariantMatrix1_ATTRIBUTES + [
    AttributeDef(
        "with_coverage", YESNO, "no", "no",
        help="Whether this file contains the coverage at each position."),
    AttributeDef(
        "with_rna_coverage", YESNO, "no", "no",
        help="Whether this file contains the coverage from RNA-Seq data."),
    ]

_SimpleVariantMatrix3_ATTRIBUTES = \
    _SimpleVariantMatrix2_ATTRIBUTES + [
    AttributeDef(
        "with_gxp", YESNO, "no", "no",
        help="Whether this file contains the expression of the genes."),
    AttributeDef(
        "with_cancer_genes", YESNO, "no", "no",
        help="Whether each gene is a known cancer gene."),
    ]
                                 
SimpleVariantMatrix_ATTRIBUTES = \
    _SimpleVariantMatrix3_ATTRIBUTES + [
    AttributeDef(
        "filtered_variants", YESNO, "no", "no",
        help="Whether the variants are filtered.  Requires XXX"),
    ]
    

UnprocessedSimpleVariantMatrix = DataType(
    "UnprocessedSimpleVariantMatrix",
    *UnprocessedSimpleVariantMatrix_ATTRIBUTES,
    **{
        "help" :
        "Contains information about variants.  Coordinates are 1-based.",
        }
    )

_SimpleVariantMatrix1 = DataType(
    "_SimpleVariantMatrix1",
    *_SimpleVariantMatrix1_ATTRIBUTES,
    **{"help" : "Should not be used by the user.",}
    )

_SimpleVariantMatrix2 = DataType(
    "_SimpleVariantMatrix2",
    *_SimpleVariantMatrix2_ATTRIBUTES,
    **{"help" : "Should not be used by the user.",}
    )

_SimpleVariantMatrix3 = DataType(
    "_SimpleVariantMatrix3",
    *_SimpleVariantMatrix3_ATTRIBUTES,
    **{"help" : "Should not be used by the user.",}
    )

SimpleVariantMatrix = DataType(
    "SimpleVariantMatrix",
    *SimpleVariantMatrix_ATTRIBUTES,
    **{
        "help" :
        "Contains information about variants.  Coordinates are 1-based.",
        }
    )

#SimpleCallMatrix = DataType(
#    "SimpleCallMatrix",
#    AttributeDef(
#        "vartype", ["snp", "indel"], "snp", "snp",
#        help="What kind of variants are held in this file."),
#    help="Contains variant calls.  Coordinates are 1-based."
#    )


all_data_types=[
    VCFFolder,
    ManySampleVCFFile,
    ManyCallerVCFFolders,
    SimpleVariantFile,
    UnprocessedSimpleVariantMatrix,
    _SimpleVariantMatrix1,
    _SimpleVariantMatrix2,
    _SimpleVariantMatrix3,
    SimpleVariantMatrix,
    #SimpleCallMatrix,

    VCFRecalibrationReport,
    PileupSummary,
    PositionsFile,
    NormalCancerFile,
    IntervalListFile,
    PositionSpecificDepthOfCoverage,

    # DEPRECATED
    AnnotatedVCFFolder,
    AnnotatedMultiVCFFile,
    ]

all_modules = [
    ModuleNode(
        "call_variants_mpileup",
        [NGS.BamFolder, NGS.ReferenceGenome], VCFFolder,
        #OptionDef("max_read_depth", default=100),
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("caller", SET_TO, "mpileup"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "all"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Use samtools mpileup to call variants."),
    ModuleNode(
        "summarize_consensus_mpileup",
        [NGS.BamFolder, NGS.ReferenceGenome, PositionsFile], PileupSummary,
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("has_read_groups", MUST_BE, "yes"),
        Constraint("indexed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Constraint("coordinates_from", CAN_BE_ANY_OF, COORDINATES_FROM, 2),
        Consequence("coordinates_from", SAME_AS_CONSTRAINT, 2),
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO, 2),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT, 2),

        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Consequence("caller", SET_TO, "mpileup"),
        #Consequence("mpileup_summary", SET_TO, "yes"),
        #Consequence("vartype", SET_TO_ONE_OF, ["all", "consensus"]),
        #Constraint("vartype", CAN_BE_ANY_OF, VARTYPE_NOT_CONSENSUS, 2),
        Constraint("vartype", CAN_BE_ANY_OF, VARTYPES, 2),
        #Constraint("vartype", SAME_AS, 0, 2),
        Consequence("vartype", SAME_AS_CONSTRAINT, 2),
        #Consequence("vartype", SET_TO, "consensus"),
        help="Use mpileup to summarize mapped reads."),
    ModuleNode(
        "summarize_reads_mpileup",
        [NGS.BamFolder, NGS.ReferenceGenome], PileupSummary,
        Constraint("sorted", MUST_BE, "coordinate", 0),
        #Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, ["yes", "no"], 0),
        Constraint("has_read_groups", MUST_BE, "yes"),
        Constraint("indexed", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        #Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"], 2),
        Consequence("vartype", SET_TO, "all"),
        Consequence("coordinates_from", SET_TO, "whole_genome"),
        help="Use mpileup to summarize mapped reads."),

    ModuleNode(
        "call_variants_GATK",
        [NGS.BamFolder, NGS.ReferenceGenome], VCFFolder,
        
        # Pipeline: read groups -> sort -> mark dups -> realign ->
        # recalibrate -> call variants
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Constraint(
            "base_quality_recalibrated", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        
        Consequence("caller", SET_TO, "gatk"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel", "all"]),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Use GATK HaplotypeCaller to call variants."),
    ModuleNode(
        "call_variants_platypus",
        [NGS.BamFolder, NGS.ReferenceGenome], VCFFolder,
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        
        Consequence("caller", SET_TO, "platypus"),
        #Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel", "all"]),
        Consequence("vartype", SET_TO, "all"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        help="Use Platypus to call variants."),

    ModuleNode(
        "filter_variants_platypus",
        VCFFolder, VCFFolder,

        Constraint("caller", MUST_BE, "platypus"),
        Consequence("caller", SAME_AS_CONSTRAINT),
        Constraint("vartype", MUST_BE, "all"),
        Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel"]),
        help="Filter the variants from Platypus."),

    ModuleNode(
        "call_consensus_varscan",
        PileupSummary, VCFFolder,
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("vartype", MUST_BE, "consensus"),
        Constraint("vartype", CAN_BE_ANY_OF, VARTYPES),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Consequence("is_consensus", SET_TO, "yes"),
        #Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("caller", SET_TO, "varscan"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("backfilled", SET_TO, "consensus"),
        Constraint("coordinates_from", CAN_BE_ANY_OF, COORDINATES_FROM),
        Consequence("coordinates_from", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        help="Use Varscan to generate consensus information."),

    ModuleNode(
        "call_variants_varscan",
        PileupSummary, VCFFolder,
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        #Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("coordinates_from", MUST_BE, "whole_genome"),
        Constraint("vartype", MUST_BE, "all"),
        Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel"]),
        #Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("caller", SET_TO, "varscan"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        # Since this PileupSummary doesn't come from calls, is not yet
        # filtered.
        Constraint("filtered_calls", MUST_BE, "no"),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        help="Use Varscan to call variants."),

    ModuleNode(
        "call_somatic_varscan",
        [PileupSummary, NormalCancerFile], VCFFolder,
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("coordinates_from", MUST_BE, "whole_genome"),
        Constraint("vartype", MUST_BE, "all", 0),
        Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel"]),
        #Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"], 0),
        #Consequence("vartype", SAME_AS_CONSTRAINT),
        #Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("caller", SET_TO, "varscan"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("somatic", SET_TO, "yes"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        # Since this PileupSummary doesn't come from calls, is not yet
        # filtered.
        Constraint("filtered_calls", MUST_BE, "no"),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        help="Use Varscan to call variants."),

    ModuleNode(
        "make_full_genome_intervals",
        NGS.ReferenceGenome, IntervalListFile,
        ),

    ModuleNode(
        "call_variants_mutect",
        [NGS.BamFolder, NormalCancerFile, NGS.ReferenceGenome,
         IntervalListFile], VCFFolder,
        OptionDef("mutect_cosmic_vcf"),
        OptionDef("mutect_dbsnp_vcf"),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("dict_added", MUST_BE, "yes", 2),
        Constraint("samtools_indexed", MUST_BE, "yes", 2),
        Consequence("caller", SET_TO, "mutect"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "snp"),
        Consequence("somatic", SET_TO, "yes"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Use MuTect to call variants.",
        ),

    ModuleNode(
        "call_variants_strelka",
        [NGS.BamFolder, NormalCancerFile, NGS.ReferenceGenome],
        VCFFolder,
        OptionDef(
            "strelka_skip_depth_filter", default="no",
            help='Set to "yes" for exome or other targeted sequencing.'),
        # Not sure if the BAM files need to be sorted or indexed.
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("aligner", MUST_BE, "bwa_mem", 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        Consequence("caller", SET_TO, "strelka"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel"]),
        Consequence("somatic", SET_TO, "yes"),
        help="Use Strelka to call variants.",
        ),

    ModuleNode(
        "call_variants_somaticsniper",
        [NGS.BamFolder, NormalCancerFile, NGS.ReferenceGenome],
        VCFFolder,
        # Not sure if the BAM files need to be sorted or indexed.
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("caller", SET_TO, "somaticsniper"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "snp"),
        Consequence("somatic", SET_TO, "yes"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Use SomaticSniper to call variants.",
        ),

    ModuleNode(
        "call_variants_jointsnvmix",
        [NGS.BamFolder, NormalCancerFile, NGS.ReferenceGenome],
        VCFFolder,
        # BAM files need to be indexed, duplicates_marked.
        # Reference must be samtools_indexed.
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("samtools_indexed", MUST_BE, "yes", 2),
        Consequence("caller", SET_TO, "jointsnvmix"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO_ONE_OF, ["snp", "indel", "all"]),
        Consequence("somatic", SET_TO, "yes"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Use JointSNVMix (museq) to call variants.",
        ),

    ModuleNode(
        "call_variants_muse",
        [NGS.BamFolder, NormalCancerFile, NGS.ReferenceGenome],
        VCFFolder,
        OptionDef("wgs_or_wes", help='Should be "wgs" or "wes".'),
        OptionDef(
            "muse_dbsnp_vcf", help="bgzip and indexed VCF file for dbsnp"),
        # Reference must be samtools_indexed.
        # BAM files need to be indexed, duplicates_marked.
        # Normal/Cancer should be re-aligned jointly.  How??
        # dbSNP VCF file must be bgzip compressed, tabix indexed.
        # Recommended aligners: bwa_mem.
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("samtools_indexed", MUST_BE, "yes", 2),
        Consequence("caller", SET_TO, "muse"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "snp"),
        Consequence("somatic", SET_TO, "yes"),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS, 0),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        help="Use MuSE to call variants.",
        ),

    ModuleNode(
        "call_variants_radia_with_rna",
        [NGS.BamFolder, NGS.BamFolder, NormalCancerFile, NGS.ReferenceGenome],
        VCFFolder,
        # Reference must be samtools_indexed.
        # BAM files need to be indexed.
        
        OptionDef(
            "radia_genome_assembly", help="Which genome assembly, e.g. hg19"),
        OptionDef(
            "snp_eff_genome",
            help="Which snpEff genome to use, e.g. GRCh37.75"),

        # DNA BamFolder
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("has_read_groups", CAN_BE_ANY_OF, YESNO, 0),
        Constraint(
            "aligner", CAN_BE_ANY_OF,
            # TODO: replace with constant: DNA_ALIGNERS
            ["bowtie1", "bowtie2", "bwa_backtrack", "bwa_mem"], 0),
        Consequence("aligner", SAME_AS_CONSTRAINT, 0),
        # RNA BamFolder
        Constraint("sorted", MUST_BE, "coordinate", 1),
        Constraint("indexed", MUST_BE, "yes", 1),
        Constraint(
            # TODO: replace with constant: RNA_ALIGNERS
            "aligner", CAN_BE_ANY_OF, ["tophat", "star"], 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 3),
        Consequence("caller", SET_TO, "radia"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "all"),
        Consequence("somatic", SET_TO, "yes"),
        help="Use Radia to call SNPs and INDELs using DNA and RNA.",
        ),

    ModuleNode(
        "select_radia_snps",
        VCFFolder, VCFFolder,
        Constraint("caller", MUST_BE, "radia"), 
        Consequence("caller", SAME_AS_CONSTRAINT),
        Constraint("vartype", MUST_BE, "all"),
        Consequence("vartype", SET_TO, "snp"),
        ),

    ModuleNode(
        "merge_somatic_variants_snp",
        [
            VCFFolder, # MuTect
            VCFFolder, # Varscan
            VCFFolder, # Strelka
            VCFFolder, # SomaticSniper
            VCFFolder, # JointSNVMix
            VCFFolder, # MuSE
            VCFFolder, # Radia
            ],
        ManyCallerVCFFolders,
        Consequence("somatic", SET_TO, "yes"),
        Constraint("caller", MUST_BE, "mutect", 0),
        Constraint("vartype", MUST_BE, "snp", 0),
        Constraint("somatic", MUST_BE, "yes", 0),
        Constraint("caller", MUST_BE, "varscan", 1),
        Constraint("vartype", MUST_BE, "snp", 1),
        Constraint("somatic", MUST_BE, "yes", 1),
        Constraint("caller", MUST_BE, "strelka", 2),
        Constraint("vartype", MUST_BE, "snp", 2),
        Constraint("somatic", MUST_BE, "yes", 2),
        Constraint("caller", MUST_BE, "somaticsniper", 3),
        Constraint("vartype", MUST_BE, "snp", 3),
        Constraint("somatic", MUST_BE, "yes", 3),
        Constraint("caller", MUST_BE, "jointsnvmix", 4),
        Constraint("vartype", MUST_BE, "snp", 4),
        Constraint("somatic", MUST_BE, "yes", 4),
        Constraint("caller", MUST_BE, "muse", 5),
        Constraint("vartype", MUST_BE, "snp", 5),
        Constraint("somatic", MUST_BE, "yes", 5),
        Constraint("caller", MUST_BE, "radia", 6),
        Constraint("vartype", MUST_BE, "snp", 6),
        Constraint("somatic", MUST_BE, "yes", 6),
        Consequence("vartype", SET_TO, "snp"),
        help="Call variants with all implemented somatic variant callers.",
        ),

    ModuleNode(
        "merge_variants_snp",
        [
            VCFFolder, # GATK
            VCFFolder, # Platypus
            VCFFolder, # Varscan
            ],
        ManyCallerVCFFolders,
        Consequence("vartype", SET_TO, "snp"),
        Consequence("somatic", SET_TO, "no"),
        Constraint("caller", MUST_BE, "gatk", 0),
        Constraint("vartype", MUST_BE, "snp", 0),
        Constraint("somatic", MUST_BE, "no", 0),
        Constraint("caller", MUST_BE, "platypus", 1),
        Constraint("vartype", MUST_BE, "snp", 1),
        Constraint("somatic", MUST_BE, "no", 1),
        Constraint("caller", MUST_BE, "varscan", 2),
        Constraint("vartype", MUST_BE, "snp", 2),
        Constraint("somatic", MUST_BE, "no", 2),
        help="Call variants with all implemented non-somatic variant callers.",
        ),

    ModuleNode(
        "merge_variants_indel",
        [
            VCFFolder, # GATK
            VCFFolder, # Platypus
            VCFFolder, # Varscan
            ],
        ManyCallerVCFFolders,
        Consequence("somatic", SET_TO, "no"),
        Consequence("vartype", SET_TO, "indel"),
        Constraint("caller", MUST_BE, "gatk", 0),
        Constraint("vartype", MUST_BE, "indel", 0),
        Constraint("somatic", MUST_BE, "no", 0),
        Constraint("caller", MUST_BE, "platypus", 1),
        Constraint("vartype", MUST_BE, "indel", 1),
        Constraint("somatic", MUST_BE, "no", 1),
        Constraint("caller", MUST_BE, "varscan", 2),
        Constraint("vartype", MUST_BE, "indel", 2),
        Constraint("somatic", MUST_BE, "no", 2),
        help="Call variants with all implemented non-somatic variant callers.",
        ),

    ModuleNode(
        "merge_manycallervcffolders",
        ManyCallerVCFFolders, SimpleVariantFile,
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "filter_simplevariantfile",
        SimpleVariantFile, SimpleVariantFile,
        OptionDef(
            "remove_samples", default="",
            help="Comma-separated list of samples to be removed."),
        OptionDef(
            "remove_radia_rna_samples", default="no",
            help="Radia generates calls from the RNA-Seq, with sample names: "
            "<sample>_RNA.  Remove these."),
        OptionDef(
            "apply_filter", default="yes",
            help='Whether to remove variants according to "Filter" column.'),
        OptionDef(
            "wgs_or_wes",
            help='Should be "wgs" or "wes".  Used for filtering MuSE calls.'),
        
        Constraint("filtered", MUST_BE, "no"),
        Consequence("filtered", SET_TO, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "convert_simplevariantfile_to_matrix",
        SimpleVariantFile, UnprocessedSimpleVariantMatrix,
        Constraint("filtered", MUST_BE, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "_convert_unprocessedsimplevariantmatrix_to_simplevariantmatrix1",
        UnprocessedSimpleVariantMatrix, _SimpleVariantMatrix1,
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        #Constraint("annotated", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("annotated", SAME_AS_CONSTRAINT),
        #Constraint("filtered", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("filtered", SAME_AS_CONSTRAINT),
        #Constraint("with_gxp", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_gxp", SAME_AS_CONSTRAINT),
        #Constraint("with_coverage", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_coverage", SAME_AS_CONSTRAINT),
        #Constraint("with_rna_coverage", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_rna_coverage", SAME_AS_CONSTRAINT),
        #Constraint("with_cancer_genes", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_cancer_genes", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "_convert_simplevariantmatrix1_to_simplevariantmatrix2",
        _SimpleVariantMatrix1, _SimpleVariantMatrix2,
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("annotated", CAN_BE_ANY_OF, YESNO),
        Consequence("annotated", SAME_AS_CONSTRAINT),
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        #Constraint("filtered_variants", CAN_BE_ANY_OF, YESNO),
        #Consequence("filtered_variants", SAME_AS_CONSTRAINT),
        #Constraint("with_gxp", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_gxp", SAME_AS_CONSTRAINT),
        #Constraint("with_rna_coverage", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_rna_coverage", SAME_AS_CONSTRAINT),
        #Constraint("with_coverage", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_coverage", SAME_AS_CONSTRAINT),
        #Constraint("with_cancer_genes", CAN_BE_ANY_OF, ["no", "yes"]),
        #Consequence("with_cancer_genes", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "_convert_simplevariantmatrix2_to_simplevariantmatrix3",
        _SimpleVariantMatrix2, _SimpleVariantMatrix3,
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("annotated", CAN_BE_ANY_OF, YESNO),
        Consequence("annotated", SAME_AS_CONSTRAINT),
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        #Constraint("filtered_variants", CAN_BE_ANY_OF, YESNO),
        #Consequence("filtered_variants", SAME_AS_CONSTRAINT),
        #Constraint("with_gxp", CAN_BE_ANY_OF, YESNO),
        #Consequence("with_gxp", SAME_AS_CONSTRAINT),
        Constraint("with_coverage", CAN_BE_ANY_OF, YESNO),
        Consequence("with_coverage", SAME_AS_CONSTRAINT),
        Constraint("with_rna_coverage", CAN_BE_ANY_OF, YESNO),
        Consequence("with_rna_coverage", SAME_AS_CONSTRAINT),
        #Constraint("with_cancer_genes", CAN_BE_ANY_OF, YESNO),
        #Consequence("with_cancer_genes", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "_convert_simplevariantmatrix3_to_simplevariantmatrix",
        _SimpleVariantMatrix3, SimpleVariantMatrix,
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("annotated", CAN_BE_ANY_OF, YESNO),
        Consequence("annotated", SAME_AS_CONSTRAINT),
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        #Constraint("filtered_variants", CAN_BE_ANY_OF, YESNO),
        #Consequence("filtered_variants", SAME_AS_CONSTRAINT),
        Constraint("with_gxp", CAN_BE_ANY_OF, YESNO),
        Consequence("with_gxp", SAME_AS_CONSTRAINT),
        Constraint("with_rna_coverage", CAN_BE_ANY_OF, YESNO),
        Consequence("with_rna_coverage", SAME_AS_CONSTRAINT),
        Constraint("with_coverage", CAN_BE_ANY_OF, YESNO),
        Consequence("with_coverage", SAME_AS_CONSTRAINT),
        Constraint("with_cancer_genes", CAN_BE_ANY_OF, YESNO),
        Consequence("with_cancer_genes", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "annotate_simplevariantmatrix",
        _SimpleVariantMatrix1, _SimpleVariantMatrix1,
        OptionDef("annovar_buildver", help="E.g. hg19.  See annovar docs."),
        Constraint("annotated", MUST_BE, "no"),
        Consequence("annotated", SET_TO, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "filter_simplevariantmatrix_calls",
        _SimpleVariantMatrix1, _SimpleVariantMatrix1,
        OptionDef(
            "filter_by_min_alt_reads", default="0",
            help="Remove calls with less than this number of ALT reads.  "
            "Does not work for Strelka because it doesn't report ALT reads."),
        OptionDef(
            "filter_by_min_total_reads", default="0",
            help="Remove calls with less than this number of reads."),
        OptionDef(
            "filter_by_min_vaf", default="0.0",
            help="Remove calls whose variant allele frequency is less than "
            "this."),
        #OptionDef(
        #    "filter_by_min_GQ", default="0",
        #    help="Remove variants with GQ less than this number.  "
        #    "SomaticSniper provides this.  "
        #    "JointSNVMix, Strelka don't use this.  "
        #    "MuTect and VarScan has it in the header, but doesn't actually "
        #    "generate any values."
        #    ),
        Constraint("filtered_calls", MUST_BE, "no"),
        Consequence("filtered_calls", SET_TO, "yes"),
        #Constraint("annotated", MUST_BE, "yes"),
        #Consequence("annotated", SAME_AS_CONSTRAINT),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "filter_simplevariantmatrix_variants",
        SimpleVariantMatrix, SimpleVariantMatrix,
        OptionDef(
            "min_callers_in_every_sample", default="",
            help="Remove variants that do not have this many callers in "
            "every sample."),
        OptionDef(
            "min_callers_in_any_sample", default="",
            help="Remove variants that do not have this many callers in "
            "any sample."),
        OptionDef(
            "min_gene_expression_in_every_sample", default="",
            help="Remove variants that do not have this many callers in "
            "every sample."),
        OptionDef(
            "min_coverage_in_every_sample", default="",
            help="Remove variants that do not have this coverage "
            "every sample."),
        OptionDef(
            "nonsynonymous_and_stopgain_only", default="no",
            help="Keep only non-synonymous and stopgain variants.  "
            '"yes" or "no".'),
        OptionDef(
            "sift_polyphen_damaging", default="no",
            help="Only if SIFT (D)amaging and Polyphen2 (P) or (D).  "
            '"yes" or "no".'),
        #OptionDef(
        #    "filter_by_min_total_reads", default="0",
        #    help="Remove variants with less than this number of reads."),
        #OptionDef(
        #    "nonsynonymous_and_stopgain_only", 
        #    help="Keep only non-synonymous and stopgain variants.  "
        #    '"yes" or "no".'),
        #OptionDef(
        #    "filter_by_min_GQ", default="0",
        #    help="Remove variants with GQ less than this number.  "
        #    "SomaticSniper provides this.  "
        #    "JointSNVMix, Strelka don't use this.  "
        #    "MuTect and VarScan has it in the header, but doesn't actually "
        #    "generate any values."
        #    ),
        Constraint("filtered_variants", MUST_BE, "no"),
        Consequence("filtered_variants", SET_TO, "yes"),
        #Constraint("annotated", MUST_BE, "yes"),
        #Consequence("annotated", SAME_AS_CONSTRAINT),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "extract_positions_from_simplevariantmatrix",
        _SimpleVariantMatrix2, PositionsFile,
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Consequence("coordinates_from", SET_TO, "simplevariantmatrix"),
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "add_coverage_to_simplevariantmatrix",
        [_SimpleVariantMatrix2, PositionSpecificDepthOfCoverage],
        _SimpleVariantMatrix2,
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("filtered_calls", SAME_AS, 0, 1),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        Constraint("with_coverage", MUST_BE, "no", 0),
        Consequence("with_coverage", SET_TO, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"], 0),
        Constraint("vartype", SAME_AS, 0, 1),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("coordinates_from", MUST_BE, "simplevariantmatrix", 1),
        Constraint(
            "aligner", CAN_BE_ANY_OF,
            # TODO: replace with constant: DNA_ALIGNERS
            ["bowtie1", "bowtie2", "bwa_backtrack", "bwa_mem"], 1),
        help="Add the coverage and baseline variant allele frequencies at "
        "each position."
        ),

    ModuleNode(
        "add_rna_coverage_to_simplevariantmatrix",
        [_SimpleVariantMatrix2, PositionSpecificDepthOfCoverage],
        _SimpleVariantMatrix2,
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("filtered_calls", SAME_AS, 0, 1),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        Constraint("with_rna_coverage", MUST_BE, "no", 0),
        Consequence("with_rna_coverage", SET_TO, "yes"),
        Constraint("with_coverage", MUST_BE, "yes", 0),
        Consequence("with_coverage", SAME_AS_CONSTRAINT),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"], 0),
        Constraint("vartype", SAME_AS, 0, 1),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint(
            # TODO: replace with constant: RNA_ALIGNERS
            "aligner", CAN_BE_ANY_OF, ["tophat", "star"], 1),
        Constraint("coordinates_from", MUST_BE, "simplevariantmatrix", 1),
        help="Add the coverage and baseline variant allele frequencies "
        "from RNA-Seq data at each position."
        ),

    ModuleNode(
        "add_gene_expression_to_simplevariantmatrix",
        [_SimpleVariantMatrix3, GXP.SignalFile], _SimpleVariantMatrix3,
        Constraint("with_gxp", MUST_BE, "no", 0),
        Consequence("with_gxp", SET_TO, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"], 0),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes", 1),
        Constraint("preprocess", MUST_BE, "tpm", 1),
        ),

    #ModuleNode(
    #    "add_misc_annotations_to_simplevariantmatrix",
    #    SimpleVariantMatrix, SimpleVariantMatrix,
    #    OptionDef(
    #        "cancer_genes_file",
    #        help='Has "Gene ID", "Gene Symbol", '
    #        "followed by cancer gene sets."),
    #    Constraint("with_misc_annotations", MUST_BE, "no"),
    #    Consequence("with_misc_annotations", SET_TO, "yes"),
    #    Constraint("annotated", MUST_BE, "yes"),
    #    Consequence("annotated", SAME_AS_CONSTRAINT),
    #    Constraint("filtered", MUST_BE, "yes"),
    #    Consequence("filtered", SAME_AS_CONSTRAINT),
    #    Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
    #    Consequence("vartype", SAME_AS_CONSTRAINT),
    #    ),
        

    ModuleNode(
        "add_cancer_genes_to_simplevariantmatrix",
        _SimpleVariantMatrix3, _SimpleVariantMatrix3,
        OptionDef(
            "cancer_genes_file",
            help='Has "Gene ID", "Gene Symbol", '
            "followed by cancer gene sets."),
        Constraint("with_cancer_genes", MUST_BE, "no"),
        Consequence("with_cancer_genes", SET_TO, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    #ModuleNode(
    #    "convert_simplevariantmatrix_to_simplecallmatrix",
    #    SimpleVariantMatrix, SimpleCallMatrix,
    #    OptionDef(
    #        "num_callers", 
    #        help="Minimum number of callers that need to call a position."),
    #    Constraint("filtered", MUST_BE, "yes"),
    #    Constraint("annotated", MUST_BE, "yes"),
    #    Constraint("vartype", CAN_BE_ANY_OF, ["snp", "indel"]),
    #    Consequence("vartype", SAME_AS_CONSTRAINT),
    #    ),

    
    ModuleNode(
        "merge_somatic_variants_indel",
        [
            VCFFolder, # Varscan
            VCFFolder, # Strelka
            VCFFolder, # JointSNVMix
            ],
        ManyCallerVCFFolders,
        Consequence("somatic", SET_TO, "yes"),
        Constraint("caller", MUST_BE, "varscan", 0),
        Constraint("vartype", MUST_BE, "indel", 0),
        Constraint("somatic", MUST_BE, "yes", 0),
        Constraint("caller", MUST_BE, "strelka", 1),
        Constraint("vartype", MUST_BE, "indel", 1),
        Constraint("somatic", MUST_BE, "yes", 1),
        Constraint("caller", MUST_BE, "jointsnvmix", 2),
        Constraint("vartype", MUST_BE, "indel", 2),
        Constraint("somatic", MUST_BE, "yes", 2),
        Consequence("vartype", SET_TO, "indel"),
        help="Call variants with all implemented somatic variant callers.",
        ),

    ModuleNode(
        "make_vcf_recalibration_report_snp",
        [VCFFolder, NGS.ReferenceGenome], VCFRecalibrationReport,
        OptionDef("vcf_recal_dbsnp"),
        OptionDef("vcf_recal_mills_indels"),
        OptionDef("vcf_recal_1kg_indels"),
        OptionDef("vcf_recal_omni"),
        Constraint("vcf_recalibrated", MUST_BE, "no", 0),
        Constraint("vartype", CAN_BE_ANY_OF, ["all", "snp"], 0),
        Consequence("vartype", SET_TO, "snp"),
        Constraint("somatic", CAN_BE_ANY_OF, YESNO, 0),
        Consequence("somatic", SAME_AS_CONSTRAINT),
        Constraint("caller", CAN_BE_ANY_OF, CALLERS, 0),
        Consequence("caller", SAME_AS_CONSTRAINT),
        help="VariantRecalibrator",
        ),
    ModuleNode(
        "recalibrate_variants_snp",
        [VCFFolder, NGS.ReferenceGenome, VCFRecalibrationReport], VCFFolder,
        Constraint("vcf_recalibrated", MUST_BE, "no", 0),
        Consequence("vcf_recalibrated", SET_TO, "yes"),
        Constraint("caller", CAN_BE_ANY_OF, CALLERS, 0),
        Constraint("caller", SAME_AS, 0, 2),
        Consequence("caller", SAME_AS_CONSTRAINT),
        Constraint("vartype", CAN_BE_ANY_OF, ["all", "snp"], 0),
        Constraint("vartype", MUST_BE, "snp", 2),
        Consequence("vartype", SAME_AS_CONSTRAINT, 2),
        Constraint("somatic", CAN_BE_ANY_OF, YESNO, 0),
        Constraint("somatic", SAME_AS, 0, 2),
        Consequence("somatic", SAME_AS_CONSTRAINT),
        help="GATK ApplyRecalibration",
        ),
    
    ModuleNode(
        "filter_snps_only_multivcf",
        ManySampleVCFFile, ManySampleVCFFile,
        Constraint("vartype", MUST_BE, "all"),
        Consequence("vartype", SET_TO, "snp"),
        ),

    ModuleNode(
        "annotate_with_annovar_vcffolder",
        VCFFolder, AnnotatedVCFFolder,
        OptionDef("buildver", help="E.g. hg19.  See annovar docs."),
        ),

    ModuleNode(
        "merge_vcf_folder",
        VCFFolder, ManySampleVCFFile,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("caller", CAN_BE_ANY_OF, CALLERS),
        Consequence("caller", SAME_AS_CONSTRAINT),
        Constraint("vartype", CAN_BE_ANY_OF, VARTYPES),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("backfilled", CAN_BE_ANY_OF, BACKFILLS),
        Consequence("backfilled", SAME_AS_CONSTRAINT),
        ),
    
    ModuleNode(
        "annotate_multivcf_annovar",
        ManySampleVCFFile, AnnotatedMultiVCFFile,
        OptionDef("buildver", help="E.g. hg19.  See annovar docs."),
        ),
    
    ModuleNode(
        "extract_positions_from_multivcf_file",
        ManySampleVCFFile, PositionsFile,
        #Constraint("backfilled", MUST_BE, "no"),
        Constraint("caller", CAN_BE_ANY_OF, CALLERS),
        Constraint("vartype", CAN_BE_ANY_OF, VARTYPES),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        ),

    ModuleNode(
        "backfill_vcf_folder",
        # The first one will be backfilled with the information from
        # the second.
        [VCFFolder, VCFFolder], VCFFolder,
        OptionDef(
            "backfill_common_only", default="no",
            help="Backfill the samples that are in common.  Ignore samples "
            'that only occur in one file.  Must be "yes" or "no".'),
        DefaultAttributesFrom(0),
        Constraint("backfilled", MUST_BE, "no", 0),
        Consequence("backfilled", SET_TO, "yes"),
        #Constraint("vartype", CAN_BE_ANY_OF, VARTYPE_NOT_CONSENSUS, 0),
        Constraint("vartype", CAN_BE_ANY_OF, VARTYPES, 0),
        Constraint("vartype", SAME_AS, 0, 1),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        Constraint("caller", CAN_BE_ANY_OF, CALLERS, 0),
        Consequence("caller", SAME_AS_CONSTRAINT),
        Constraint("backfilled", MUST_BE, "consensus", 1),
        Constraint("is_consensus", MUST_BE, "no", 0),
        Constraint("is_consensus", MUST_BE, "yes", 1),
        Consequence("is_consensus", SAME_AS_CONSTRAINT, 1),
        #Constraint("vartype", MUST_BE, "consensus", 1),
        Constraint("caller", CAN_BE_ANY_OF, CALLERS, 1),
        ),
    #ModuleNode(
    #    "backfill_multivcf_file",
    #    [AnnotatedMultiVCFFile, MultiVCFFile],
    #    AnnotatedMultiVCFFile,
    #    DefaultAttributesFrom(0),
    #    Constraint("backfilled", MUST_BE, "no", 0),
    #    Consequence("backfilled", SET_TO, "yes"),
    #    Constraint("vartype", MUST_BE, "consensus", 1),
    #    Constraint("caller", MUST_BE, "varscan", 1),
    #    ),

    ModuleNode(
        "summarize_coverage_at_positions",
        VCFFolder, PositionSpecificDepthOfCoverage,
        Constraint("filtered_calls", CAN_BE_ANY_OF, YESNO),
        Consequence("filtered_calls", SAME_AS_CONSTRAINT),
        Constraint("backfilled", MUST_BE, "consensus"),
        Constraint("caller", MUST_BE, "varscan"),
        Constraint("is_consensus", MUST_BE, "yes"),
        Constraint("vartype", CAN_BE_ANY_OF, VARTYPES),
        Consequence("vartype", SAME_AS_CONSTRAINT),
        #Constraint("vartype", MUST_BE, "consensus"),
        Constraint("coordinates_from", CAN_BE_ANY_OF, COORDINATES_FROM),
        Consequence("coordinates_from", SAME_AS_CONSTRAINT),
        Constraint("aligner", CAN_BE_ANY_OF, NGS.ALIGNERS),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        ),
    ]
