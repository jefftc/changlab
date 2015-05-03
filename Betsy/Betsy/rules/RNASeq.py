#RNASeq
from Betsy.bie3 import *
import BasicDataTypesNGS
import Database
import GeneExpProcessing
list_files = []
all_modules = [
    Module(
        'is_sam_folder',
         BasicDataTypesNGS.RNASeqFile,BasicDataTypesNGS.RNASeqFile,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("format_type",MUST_BE,'not_fastqfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("format_type",BASED_ON_DATA,['not_samfolder','samfolder']),
         help=("extract rna files with different format")
        ),
    Module(
        'is_bam_folder',
         BasicDataTypesNGS.RNASeqFile,BasicDataTypesNGS.RNASeqFile,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("format_type",MUST_BE,'not_samfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("format_type",BASED_ON_DATA,['not_bamfolder','bamfolder']),
         help=("extract rna files with different format")
        ),
    Module(
        'is_fastq_folder',
         BasicDataTypesNGS.RNASeqFile,BasicDataTypesNGS.RNASeqFile,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("format_type",MUST_BE,'unknown'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("format_type",BASED_ON_DATA,['not_fastqfolder','fastqfolder']),
         help=("extract rna files with different format")
        ),
    Module(
        'extract_rna_files_sam',
         BasicDataTypesNGS.RNASeqFile,BasicDataTypesNGS.SamFolder,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("format_type",MUST_BE,'samfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("ref",SET_TO_ONE_OF,['human','mouse']),
         help=("extract rna files with different format")
        ),
    Module(
        'extract_rna_files_bam',
         BasicDataTypesNGS.RNASeqFile,BasicDataTypesNGS.BamFolder,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("format_type",MUST_BE,'bamfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("ref",SET_TO_ONE_OF,['human','mouse']),
         help=("extract rna files with bam format")
        ),
    Module(
        'extract_rna_files_fastq',
         BasicDataTypesNGS.RNASeqFile,BasicDataTypesNGS.FastqFolder,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("format_type",MUST_BE,'fastqfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help=("extract rna files with fa or fastq format")
        ),
    Module(
        'align_with_bowtie',
         [BasicDataTypesNGS.FastqFolder,BasicDataTypesNGS.SampleGroupFile],
         BasicDataTypesNGS.SamFolder,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Consequence("ref",SET_TO_ONE_OF,['human','mouse']),
         help=("convert fastq folder into sam folder with bowtie")
        ),
    Module(
        'convert_sam_to_bam',
         BasicDataTypesNGS.SamFolder,BasicDataTypesNGS.BamFolder,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         #Consequence("ref",SET_TO_ONE_OF,['human','mouse']),
         help=("convert sam folder into bam folder")
        ),
     Module(
        'normalize_with_rsem',
         [BasicDataTypesNGS.BamFolder,BasicDataTypesNGS.SampleGroupFile],
         GeneExpProcessing._SignalFile_Postprocess,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Consequence("preprocess",SET_TO,"rsem"),
         Consequence("logged",SET_TO,"unknown"),
         Consequence("predataset",SET_TO,"no"),
         Consequence("format",SET_TO,"tdf"),
         help=("process BamFolder , generate SignalFile_Postprocess with preprocess rsem")
        ),


    ]
