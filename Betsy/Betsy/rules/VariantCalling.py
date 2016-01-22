#Calling_variants_rule
from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS

VCFFolder = DataType(
    "VCFFolder",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "caller", ["none", "mpileup", "gatk", "platypus"], "none", "mpileup",
        help="Which variant caller was used."),

    AttributeDef(
        "recalibrated", ["no", "yes"], "no", "yes",
        help="Whether quality scores are ready for filtering."),
    AttributeDef(
        "vartype", ["both", "snp", "indel"], "both", "both",
        help="What kind of variants are held in this file."),

    
    #AttributeDef(
    #    "base_recalibrated", ["yes", "no"], "no", "no",
    #    help="recalibrated or not"),
    #AttributeDef(
    #    "indel_realigned", ["yes", "no"], "no", "no",
    #    help="realigned or not"),
    #AttributeDef(
    #    "read", ["single", "paired"], "single", "single",
    #    help="single or pair read"),
    #AttributeDef(
    #    "vcf_filter", ["yes", "no"], "no", "no", help="filter VcfFile or not"),
    #AttributeDef(
    #    "reheader", ["standard", "bcftool"], "standard", "standard",
    #    help="method to convert to VcfFile"),
    #AttributeDef(
    #    "vcf_annotate", ["yes", "no"], "no", "no",
    #    help="annotate VcfFile or not"),
    help="VCF file"
    )


all_data_types=[
    VCFFolder,
    ]

all_modules = [
    ModuleNode(
        "call_variants_mpileup",
        NGS.BamFolder, VCFFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("duplicates_marked", MUST_BE, "yes"),
        Constraint("has_read_groups", MUST_BE, "yes"),
        #Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("caller", SET_TO, "mpileup"),
        #Consequence("vcf_filter", SET_TO, "no"),
        #Consequence("reheader", SET_TO, "bcftool"),
        #Consequence("vcf_annotate", SET_TO, "no"),
        #Consequence("recalibrated", SAME_AS_CONSTRAINT),
        help="use mpileup to call variants"),
    ModuleNode(
        "call_variants_GATK",
        [NGS.BamFolder, NGS.ReferenceGenome], VCFFolder,
        
        # Pipeline: read groups -> sort -> mark dups -> realign ->
        # recalibrate -> call variants
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("indel_realigned", MUST_BE, "yes", 0),
        Constraint("base_recalibrated", MUST_BE, "yes", 0),
        
        Consequence("caller", SET_TO, "gatk"),
        Consequence("recalibrated", SET_TO, "no"),
        Consequence("vartype", SET_TO, "both"),
        help="use GATK to call variants"),
    ## ModuleNode(
    ##     "filter_vcf_folder",
    ##     NGS.VcfFolder, NGS.VcfFolder,
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("base_recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
    ##     Constraint("vcf_annotate", MUST_BE, "no"),
    ##     Constraint("vcf_filter", MUST_BE, "no"),
    ##     Constraint("reheader", MUST_BE, "bcftool"),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence("vcf_filter", SET_TO, "yes"),
    ##     Consequence("reheader", SAME_AS_CONSTRAINT),
    ##     Consequence("vcf_annotate", SAME_AS_CONSTRAINT),
    ##     Consequence("base_recalibrated", SAME_AS_CONSTRAINT),
    ##     help="filter vcf file"),
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
