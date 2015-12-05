from Betsy.bie3 import *
import BasicDataTypes as BDT
import BasicDataTypesNGS as NGS

FastQCFolder = DataType(
    "FastQCFolder",
    help="Folder that holds FastQC results.",
    )

FastQCSummary = DataType(
    "FastQCSummary",
    help="An Excel file that merges and summarizes FastQC results."
    )

RNASeQCResults = DataType(
    "RNASeQCResults",
    )

RNASeQCSummary = DataType(
    "RNASeQCSummary",
    )

all_data_types = [
    FastQCFolder,
    FastQCSummary,
    RNASeQCResults,
    RNASeQCSummary,
    ]

all_modules = [
    ModuleNode(
        "run_fastqc",
        NGS.FastqFolder, FastQCFolder,
        # Actually, will work on gzip'd data.
        Constraint("compressed", MUST_BE, "no"),
        help="Run FastQC on a folder of FASTQ files.",
        ),
    ModuleNode(
        "summarize_fastqc_results",
        FastQCFolder, FastQCSummary,
        help="Merge and summarize the results from a FastQC folder.",
        ),
    
    ModuleNode(
        "run_RNA_SeQC",
        [NGS.BamFolder, NGS.ReferenceGenome], RNASeQCResults,
        #OptionDef(
        #    "RNA_ref", help="ref file for RNA_SeQC"),
        #OptionDef(
        #    "RNA_gtf", help="gtf file for RNA_SeQC"),
        OptionDef(
            "rna_seqc_gtf_file", 
            help="Gene annotations in GTF format.",
            ),
        #Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        #Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("ref", CAN_BE_ANY_OF, ["hg18", "hg19"]),
        Constraint("indexed", MUST_BE, "yes", 0),
        Constraint("sorted", MUST_BE, "contig", 0),
        Constraint("duplicates_marked", MUST_BE, "yes", 0),
        Constraint("has_read_groups", MUST_BE, "yes", 0),
        #Constraint("sample_type", MUST_BE, "RNA"),

        Constraint("samtools_indexed", MUST_BE, "yes", 1),
        Constraint("dict_added", MUST_BE, "yes", 1),
        
        help="run RNA-SeQC"),
    ModuleNode(
        "summarize_rnaseqc_results",
        RNASeQCResults, RNASeQCSummary,
        help="Summarize the results from a RNASeQC analysis.",
        ),
    ]
