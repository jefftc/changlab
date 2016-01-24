# DataTypes:
# VCFFolder
# VCFRecalibrationReport
# AnnovarFolder
#
# Modules:
# call_variants_mpileup
# call_variants_GATK
# call_variants_platypus
#
# 
# Recalibrate variant scores with GATK.
# https://www.broadinstitute.org/gatk/guide/article?id=2805
    

from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS

CALLERS = ["none", "mpileup", "gatk", "platypus"]

VCFFolder = DataType(
    "VCFFolder",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "caller", CALLERS, "none", "mpileup",
        help="Which variant caller was used."),

    AttributeDef(
        "vartype", ["both", "snp", "indel"], "both", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "vcf_recalibrated", ["no", "yes"], "no", "yes",
        help="Whether quality scores are ready for filtering."),
    )

VCFRecalibrationReport = DataType(
    "VCFRecalibrationReport",
    AttributeDef(
        "vartype", ["snp", "indel"], "snp", "snp",
        help="What kind of variants are held in this file."),
    AttributeDef(
        "caller", CALLERS, "none", "mpileup",
        help="Which variant caller was used."),
    )

AnnovarFolder = DataType(
    "AnnovarFolder",
    )

all_data_types=[
    VCFFolder,
    VCFRecalibrationReport,
    AnnovarFolder,
    ]

all_modules = [
    ModuleNode(
        "call_variants_mpileup",
        NGS.BamFolder, VCFFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("duplicates_marked", MUST_BE, "yes"),
        Constraint("has_read_groups", MUST_BE, "yes"),
        
        Consequence("caller", SET_TO, "mpileup"),
        Consequence("vcf_recalibrated", SET_TO, "yes"),
        Consequence("vartype", SET_TO, "snp"),
        help="use mpileup to call variants"),
    ModuleNode(
        "call_variants_GATK",
        [NGS.BamFolder, NGS.ReferenceGenome], VCFFolder,
        
        # Pipeline: read groups -> sort -> mark dups -> realign ->
        # recalibrate -> call variants
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Constraint("base_recalibrated", CAN_BE_ANY_OF, ["no", "yes"], 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        
        Consequence("caller", SET_TO, "gatk"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "both"),
        help="Use GATK HaplotypeCaller to call variants."),
    ModuleNode(
        "call_variants_platypus",
        [NGS.BamFolder, NGS.ReferenceGenome], VCFFolder,
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("sorted", MUST_BE, "coordinate", 0),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("dict_added", MUST_BE, "yes", 1),
        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Consequence("caller", SET_TO, "platypus"),
        Consequence("vcf_recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "snp"),
        help="Use GATK HaplotypeCaller to call variants."),

    ModuleNode(
        "make_vcf_recalibration_report_snp",
        [VCFFolder, NGS.ReferenceGenome], VCFRecalibrationReport,
        OptionDef("vcf_recal_dbsnp"),
        OptionDef("vcf_recal_mills_indels"),
        OptionDef("vcf_recal_1kg_indels"),
        OptionDef("vcf_recal_omni"),
        Constraint("vcf_recalibrated", MUST_BE, "no", 0),
        Constraint("vartype", CAN_BE_ANY_OF, ["both", "snp"], 0),
        Consequence("vartype", SET_TO, "snp"),
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
        Constraint("vartype", CAN_BE_ANY_OF, ["both", "snp"], 0),
        Constraint("vartype", MUST_BE, "snp", 2),
        Consequence("vartype", SAME_AS_CONSTRAINT, 2),
        help="ApplyRecalibration",
        ),
    ModuleNode(
        "annotate_with_annovar",
        VCFFolder, AnnovarFolder,
        OptionDef("buildver", help="E.g. hg19.  See annovar docs."),
        ),

    
    ## ModuleNode(
    ##     "annotate_vcf_folder",
    ##     NGS.VcfFolder, NGS.VcfFolder,
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
    ##     Constraint("vcf_annotate", MUST_BE, "no"),
    ##     Constraint("vcf_filter", MUST_BE, "yes"),
    ##     Constraint("reheader", CAN_BE_ANY_OF, ["bcftool", "standard"]),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence("vcf_filter", SAME_AS_CONSTRAINT),
    ##     Consequence("reheader", SAME_AS_CONSTRAINT),
    ##     Consequence("vcf_annotate", SET_TO, "yes"),
    ##     Consequence("recalibrated", SAME_AS_CONSTRAINT),
    ##     help="annotate vcf file"),
    ]
