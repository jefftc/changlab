#Calling_variants_rule
from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS

all_data_types=[]

all_modules = [
    ModuleNode(
        "call_variants_mpileup",
        NGS.BamFolder, NGS.VcfFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("duplicates_marked", MUST_BE, "yes"),
        Constraint("has_header", MUST_BE, "yes"),
        #Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        #Consequence("vcf_filter", SET_TO, "no"),
        # XXX What is this for?
        Consequence("reheader", SET_TO, "bcftool"),
        #Consequence("vcf_annotate", SET_TO, "no"),
        #Consequence("recalibrated", SAME_AS_CONSTRAINT),
        help="use mpileup to call variants"),
    ModuleNode(
        "call_variants_GATK",
        NGS.BamFolder, NGS.VcfFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("sorted", MUST_BE, "coordinate"),
        Constraint("duplicates_marked", MUST_BE, "yes"),
        Constraint("has_header", MUST_BE, "yes"),
        #Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("vcf_filter", SET_TO, "yes"),
        Consequence("reheader", SET_TO, "standard"),
        Consequence("vcf_annotate", SET_TO, "no"),
        #Consequence("recalibrated", SAME_AS_CONSTRAINT),
        help="use GATK to call variants"),
    ModuleNode(
        "filter_vcf_folder",
        NGS.VcfFolder, NGS.VcfFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("vcf_annotate", MUST_BE, "no"),
        Constraint("vcf_filter", MUST_BE, "no"),
        Constraint("reheader", MUST_BE, "bcftool"),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("vcf_filter", SET_TO, "yes"),
        Consequence("reheader", SAME_AS_CONSTRAINT),
        Consequence("vcf_annotate", SAME_AS_CONSTRAINT),
        Consequence("recalibrated", SAME_AS_CONSTRAINT),
        help="filter vcf file"),
    ModuleNode(
        "annotate_vcf_folder",
        NGS.VcfFolder, NGS.VcfFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
        Constraint("vcf_annotate", MUST_BE, "no"),
        Constraint("vcf_filter", MUST_BE, "yes"),
        Constraint("reheader", CAN_BE_ANY_OF, ["bcftool", "standard"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("vcf_filter", SAME_AS_CONSTRAINT),
        Consequence("reheader", SAME_AS_CONSTRAINT),
        Consequence("vcf_annotate", SET_TO, "yes"),
        Consequence("recalibrated", SAME_AS_CONSTRAINT),
        help="annotate vcf file"),
    
    ## ModuleNode(
    ##     "align_sequence",
    ##     NGS.FastqFile, NGS.SaiFile,
    ##     Constraint("read", CAN_BE_ANY_OF, ["single", "pair1", "pair2"]),
    ##     Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19", "mm9", "dm3"]),
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Consequence("read", SAME_AS_CONSTRAINT),
    ##     Consequence("ref", SAME_AS_CONSTRAINT),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     help="algin sequence in FastqFile, generate SaiFile"
    ##     ),
    #ModuleNode(
    #    "is_Bam_folder_sorted",
    #    NGS.BamFolder, NGS.BamFolder,
    #    Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    #    Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
    #    Constraint("duplicates_marked", MUST_BE, "no"),
    #    Constraint("indexed", MUST_BE, "no"),
    #    Constraint("sorted", MUST_BE, "unknown"),
    #    Constraint("sample_type", MUST_BE, "RNA"),
    #    Consequence("contents", SAME_AS_CONSTRAINT),
    #    Consequence("ref", SAME_AS_CONSTRAINT),
    #    Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
    #    Consequence("indexed", SAME_AS_CONSTRAINT),
    #    Consequence("sorted", BASED_ON_DATA, ["yes", "no"]),
    #    #Consequence("sample_type", SAME_AS_CONSTRAINT),   
    #    help="check bam folder sorted or not"
    #    ),
    ]
