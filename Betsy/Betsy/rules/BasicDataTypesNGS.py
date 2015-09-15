# SampleGroupFile
# FastqFolder
# SamFolder
# BamFolder
# VcfFolder
# SaiFolder     # from BWA

# FastqFile
# SamFile
# BamFile


from Betsy.bie3 import *
import BasicDataTypes as BDT

ALIGNERS = ["unknown", "bowtie1", "bowtie2", "bwa"]
#REFERENCE_GENOMES = ["hg18", "hg19", "mm9", "dm3"]

COMPRESSION = ["unknown", "no", "gz", "bz2", "xz"]
COMPRESSION_NOT_UNKNOWN = [x for x in COMPRESSION if x != "unknown"]

SORT_ORDERS = ["no", "coordinate", "name", "contig"]
SORTED = [x for x in SORT_ORDERS if x != "no"]
# coordinate  samtools sort
# name        by read name samtools sort -n
# contig      match contig ordering of reference genome (Picard)

ORIENTATION = ["unknown", "single", "paired_fr", "paired_rf", "paired_ff"]
ORIENTATION_NOT_UNKNOWN = [x for x in ORIENTATION if x != "unknown"]


ReferenceGenome = DataType(
    "ReferenceGenome",
    help="Should be FASTA file with reference genome.",
    )

Bowtie1IndexedGenome = DataType(
    "Bowtie1IndexedGenome",
    help="Indexed for bowtie1.",
    )

Bowtie2IndexedGenome = DataType(
    "Bowtie2IndexedGenome",
    help="Indexed for bowtie2.",
    )

Bowtie1AlignmentSummary = DataType(
    "Bowtie1AlignmentSummary",
    help="Summarizes the alignment from bowtie1.",
    )

Bowtie2AlignmentSummary = DataType(
    "Bowtie2AlignmentSummary",
    help="Summarizes the alignment from bowtie2.",
    )

BWAIndexedGenome = DataType(
    "BWAIndexedGenome",
    help="Indexed for BWA.",
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
    AttributeDef("sorted", SORT_ORDERS, "no", "no"),
    AttributeDef(
        "duplicates_marked", ["yes", "no"], "no", "no",
        help="mark duplicate or not"),
    AttributeDef(
        "recalibrated", ["yes", "no"], "no", "no",
        help="recalibrated or not"),
    AttributeDef(
        "has_header", ["yes", "no"], "no", "no",
        help="fix header or not"),
    AttributeDef(
        "has_read_groups", ["yes", "no"], "no", "no",
        help="Whether the file contains read groups.",
        ),
    ]


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
        "adapter_trimmed", ["yes", "no"], "no", "no",
        help="Whether the adapters are trimmed."),
    AttributeDef(
        "reads_merged", ["yes", "no"], "no", "no",
        help="Whether reads for a sample are merged into one file."),
    #AttributeDef(
    #    "orientation", ORIENTATION, "unknown", "unknown",
    #    help="Either single-end reads, paired-end reads with orientation.  "
    #    "See the bowtie manual for a description of the orientation.",
    #    ),
    help="A folder containing FASTQ files."
    )



SamFile = DataType("SamFile", *SAM_ATTRIBUTES)

SamFolder = DataType(
    "SamFolder", *SAM_ATTRIBUTES, help="A folder containing SAM files.")

BamFile = DataType("BamFile", *BAM_ATTRIBUTES)

BamFolder = DataType(
    "BamFolder", *BAM_ATTRIBUTES, help="A folder containing BAM files.")

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
    help=".sai file generated by BWA."
    )

VcfFolder = DataType(
    "VcfFolder",
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    AttributeDef(
        "recalibrated", ["yes", "no"], "no", "no",
        help="recalibrated or not"),
    AttributeDef(
        "read", ["single", "paired"], "single", "single",
        help="single or pair read"),
    #AttributeDef(
    #    "ref", REFERENCE_GENOMES, "hg19", "hg19",
    #    help="ref species"),
    AttributeDef(
        "vcf_filter", ["yes", "no"], "no", "no", help="filter VcfFile or not"),
    AttributeDef(
        "reheader", ["standard", "bcftool"], "standard", "standard",
        help="method to convert to VcfFile"),
    AttributeDef(
        "vcf_annotate", ["yes", "no"], "no", "no",
        help="annotate VcfFile or not"),
    help="Vcf file"
    )

all_data_types = [
    ReferenceGenome,
    Bowtie1IndexedGenome,
    Bowtie2IndexedGenome,
    BWAIndexedGenome,

    Bowtie1AlignmentSummary,
    Bowtie2AlignmentSummary,
    
    #FastqFile,
    #SamFile,
    #BamFile,
    SampleGroupFile,
    FastqFolder,
    SamFolder,
    BamFolder,
    SaiFolder,
    VcfFolder,
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
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION, 1),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT),
        ),
    ModuleNode(
        "check_single_or_paired_orientation",
        [SampleGroupFile, FastqFolder, Bowtie2IndexedGenome], SampleGroupFile,
        Constraint("orientation", MUST_BE, "unknown", 0),
        Consequence("orientation", BASED_ON_DATA, ORIENTATION_NOT_UNKNOWN),
        Constraint("reads_merged", MUST_BE, "yes", 1),
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
        "convert_sam_to_bam_folder",
        SamFolder, BamFolder,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Constraint("aligner", CAN_BE_ANY_OF, ALIGNERS),
        #Constraint("read", CAN_BE_ANY_OF, ["single", "paired"]),
        #Constraint("has_header", CAN_BE_ANY_OF, ["yes", "no"]),
        #Constraint("has_read_groups", CAN_BE_ANY_OF, ["yes", "no"]),
        #Constraint("sorted", CAN_BE_ANY_OF, ["yes", "no"]),
        #Constraint("duplicates_marked", CAN_BE_ANY_OF, ["yes", "no"]),
        #Constraint("recalibrated", CAN_BE_ANY_OF, ["yes", "no"]),
      
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SAME_AS_CONSTRAINT),
        #Consequence("read", SAME_AS_CONSTRAINT),
        Consequence("has_header", SET_TO, "no"),
        Consequence("has_read_groups", SET_TO, "no"),
        Consequence("sorted", SET_TO, "no"),
        Consequence("duplicates_marked", SET_TO, "no"),
        Consequence("recalibrated", SET_TO, "no"),

        #Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
        #help="Convert SAM to BAM files.",
        ),
    ModuleNode(
        "index_bwa_reference",
        ReferenceGenome, BWAIndexedGenome,
        OptionDef(
            "assembly", default="genome",
            help="Optional name for the genome assembly, e.g. hg19",
            ),
        ),
    ModuleNode(
        "align_with_bwa",
        #[FastqFolder, SaiFile], SamFolder,
        [FastqFolder, SampleGroupFile, BWAIndexedGenome],
        SaiFolder,
        
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        #Consequence("orientation", SAME_AS_CONSTRAINT),
        
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        #Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT),

        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),

        #Consequence("read", SAME_AS_CONSTRAINT),
        #Consequence("aligner", SET_TO, "bwa"),
        #Consequence("sorted", SET_TO, "no"),
        #Consequence("duplicates_marked", SET_TO, "no"),
        #Consequence("recalibration", SET_TO, "no"),
        #Consequence("has_header", SET_TO, "no"),
        #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19", "mm9", "dm3"]),
        #Consequence("ref", SAME_AS_CONSTRAINT),
        #help="generate algiment in SaiFile to SamFile"
        ),
    ModuleNode(
        "index_bowtie2_reference",
        ReferenceGenome, Bowtie2IndexedGenome,
        OptionDef(
            "assembly", default="genome",
            help="Optional name for the genome assembly, e.g. hg19",
            ),
        ),
    ModuleNode(
        "align_with_bowtie2",
        [FastqFolder, SampleGroupFile, Bowtie2IndexedGenome],
        SamFolder,
        #OptionDef(
        #    "orientation", default="fr",
        #    help="Which orientation.  See bowtie2 manual.  "
        #    "Values: fr, rf, ff.",
        #    ),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("aligner", SET_TO, "bowtie2"),
        #Consequence("ref", SET_TO_ONE_OF, ["human", "mouse"]),
        help="Align to a reference genome with bowtie 2.",
        ),
    ModuleNode(
        "index_bowtie1_reference",
        ReferenceGenome, Bowtie1IndexedGenome,
        OptionDef(
            "assembly", default="genome",
            help="Optional name for the genome assembly, e.g. hg19",
            ),
        ),
    ModuleNode(
        "align_with_bowtie1",
        [FastqFolder, SampleGroupFile, Bowtie1IndexedGenome],
        SamFolder,
        #OptionDef(
        #    "bowtie1_orientation", default="fr",
        #    help="Which orientation.  See bowtie manual.  "
        #    "Values: fr, rf, ff.",
        #    ),
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 1),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
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
        "convert_sai_to_sam_folder",
        [FastqFolder, SaiFolder, BWAIndexedGenome, SampleGroupFile], SamFolder,

        #OptionDef(
        #    "bwa_reverse_orientation", default="false",
        #    help='Set to "true" to flip the left and right Fastq files '
        #    'for each pair.',
        #    ),
        Constraint("compressed", MUST_BE, "no", 0),
        Constraint("reads_merged", MUST_BE, "yes", 0),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS, 0),
        Constraint("contents", SAME_AS, 0, 1),
        Constraint("contents", SAME_AS, 0, 3),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("aligner", SET_TO, "bwa"),
        
        Constraint("orientation", CAN_BE_ANY_OF, ORIENTATION_NOT_UNKNOWN, 3),
        #Constraint("orientation", SAME_AS, 0, 1),
        
        #Constraint("read", MUST_BE, "pair1", 0),
        #Constraint("read", MUST_BE, "pair2", 1),
        #Constraint("read", MUST_BE, "pair1", 2),
        #Constraint("read", MUST_BE, "pair2", 3),
        #Consequence("read", SET_TO, "single"),
        
        #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19", "mm9", "dm3"], 0),
        #Constraint("ref", SAME_AS,0,1),
        #Constraint("ref", SAME_AS,0,2),
        #Constraint("ref", SAME_AS,0,3),
        #Consequence("ref", SAME_AS_CONSTRAINT,0),
        
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("contents", SAME_AS, 0, 1),
        #Constraint("contents", SAME_AS,0,2),
        #Constraint("contents", SAME_AS,0,3),
        
        #Consequence("sorted", SET_TO, "no"),
        #Consequence("duplicates_marked", SET_TO, "no"),
        #Consequence("recalibration", SET_TO, "no"),
        #Consequence("has_header", SET_TO, "no"),
        help="Convert bwa's .sai alignments into .sam format.",
        ),
    #ModuleNode(
    #    "convert_sai_to_sam_folder_paired",
    #    SaiFolder, SamFolder,
    #    Consequence("read", SET_TO, "paired"),
    #    Consequence("aligner", SET_TO, "bwa"),
    #    ),

    ModuleNode(
        "sort_bam_folder_by_coordinate",
        BamFolder, BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("duplicates_marked", MUST_BE, "no"),
        #Constraint("recalibrated", MUST_BE, "no"),
        #Constraint("has_header", MUST_BE, "no"),
        Constraint("sorted", CAN_BE_ANY_OF, ["no", "name", "contig"]),
        Constraint("indexed", MUST_BE, "no"),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("sorted", SET_TO, "coordinate"),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        #Consequence("has_header", SAME_AS_CONSTRAINT),
        #Consequence("recalibrated", SAME_AS_CONSTRAINT),
        #Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        help="sort sam file and generate bam file "
        ),
    ModuleNode(
        "add_read_groups_to_bam_folder",
        BamFolder, BamFolder,
        Constraint("has_read_groups", MUST_BE, "no"),
        Consequence("has_read_groups", SET_TO, "yes"),
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
        Constraint("indexed", MUST_BE, "no"),
        Consequence("indexed", SAME_AS_CONSTRAINT),
        help="mark duplicates in SamFile"
        ),
    ModuleNode(
        "fix_header_GATK",
        #BamFile, BamFile,
        BamFolder, BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("duplicates_marked", MUST_BE, "yes"),
        #Constraint("recalibrated", MUST_BE, "no"),
        Constraint("has_header", MUST_BE, "no"),
        #Constraint("sorted", MUST_BE, "yes"),
        #Consequence("sorted", SAME_AS_CONSTRAINT),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("has_header", SET_TO, "yes"),
        #Consequence("recalibrated", SAME_AS_CONSTRAINT),
        #Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        help="use GATK to fix header"
        ),
    ## ModuleNode(
    ##     "recalibrate_base_quality_score",
    ##     #BamFile, BamFile,
    ##     BamFolder, BamFolder,
        
    ##     Constraint("recalibrated", MUST_BE, "no"),
    ##     Consequence("recalibrated", SET_TO, "yes"),

    ##     Constraint("sorted", MUST_BE, "yes"),
    ##     Consequence("sorted", SAME_AS_CONSTRAINT),

    ##     #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
    ##     #Constraint("duplicates_marked", MUST_BE, "yes"),
    ##     #Constraint("has_header", MUST_BE, "yes"),
    ##     #Consequence("contents", SAME_AS_CONSTRAINT),
    ##     #Consequence("has_header", SAME_AS_CONSTRAINT),
    ##     #Consequence("ref", SAME_AS_CONSTRAINT),
    ##     #Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
    ##     help="recalibrated sam file"
    ##     ),
    ## ModuleNode(
    ##     "sort_bam_folder",
    ##     BamFolder, BamFolder,
    ##     Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
    ##     Constraint("duplicates_marked", MUST_BE, "no"),
    ##     Constraint("indexed", MUST_BE, "no"),
    ##     Constraint("sorted", MUST_BE, "no"),
    ##     Constraint("sample_type", MUST_BE, "RNA"),
    ##     Consequence("contents", SAME_AS_CONSTRAINT),
    ##     Consequence("ref", SAME_AS_CONSTRAINT),
    ##     Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
    ##     Consequence("indexed", SAME_AS_CONSTRAINT),
    ##     Consequence("sorted", SET_TO, "yes"),
    ##     Consequence("sample_type", SAME_AS_CONSTRAINT),   
    ##     help="sort bam folder"  ),
    
    ## ModuleNode(
    ##     "flag_dups_in_bam_folder",
    ##     BamFolder, BamFolder,
    ##     #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    ##     #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
    ##     Constraint("duplicates_marked", MUST_BE, "no"),
    ##     #Constraint("indexed", MUST_BE, "no"),
    ##     #Constraint("sorted", MUST_BE, "yes"),
    ##     #Constraint("sample_type", MUST_BE, "RNA"),
    ##     #Consequence("contents", SAME_AS_CONSTRAINT),
    ##     #Consequence("ref", SAME_AS_CONSTRAINT),
    ##     Consequence("duplicates_marked", SET_TO, "yes"),
    ##     #Consequence("indexed", SAME_AS_CONSTRAINT),
    ##     #Consequence("sorted", SAME_AS_CONSTRAINT),
    ##     #Consequence("sample_type", SAME_AS_CONSTRAINT),   
    ##     help="mark duplicates in bam folder",
    ##     ),
    
    ModuleNode(
        "index_bam_folder",
        BamFolder, BamFolder,
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
        #Constraint("duplicates_marked", MUST_BE, "yes"),
        Constraint("indexed", MUST_BE, "no"),
        Constraint("sorted", CAN_BE_ANY_OF, SORTED),
        #Constraint("sample_type", MUST_BE, "RNA"),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Consequence("ref", SAME_AS_CONSTRAINT),
        #Consequence("duplicates_marked", SAME_AS_CONSTRAINT),
        Consequence("indexed", SET_TO, "yes"),
        Consequence("sorted", SAME_AS_CONSTRAINT),
        #Consequence("sample_type", SAME_AS_CONSTRAINT),   
        help="index bam folder",
        ),
    ]
