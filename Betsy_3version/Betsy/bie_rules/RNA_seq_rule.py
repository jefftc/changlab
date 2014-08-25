#RNA_seq_rule
from Betsy.bie3 import *
import SignalFile_rule

RNA_SeqFile=DataType(
    'RNA_SeqFile',
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    AttributeDef("format_type",['samfolder','bamfolder','fastqfolder'],
                 'samfolder','samfolder',help="format type"),
    help="RNA Seq File"
    )
FastqFolder=DataType(
    'FastqFolder',
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="RNA seq Fastq folder"
    )
SamFolder=DataType(
    'SamFolder',
    AttributeDef("ref",['human','mouse'], "human","human",
                 help='ref species'),
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    AttributeDef("sample_type",['RNA','DNA'],
                 'RNA','RNA',help="RNA or DNA type"),
    help="RNA seq Sam folder"
    )

BamFolder=DataType(
    'BamFolder',
    AttributeDef("ref",['human','mouse'], "human","human",
                 help='ref species'),
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    AttributeDef("sample_type",['RNA','DNA'],
                 'RNA','RNA',help="RNA or DNA type"),
    help="RNA seq Bam folder"
    )
SampleGroupFile=DataType(
    'SampleGroupFile',
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="File contains sample group infomation"
    )
list_files = [RNA_SeqFile, FastqFolder, SamFolder, BamFolder,SampleGroupFile]

all_modules = [
    Module(
        'extract_rna_files_sam',
         RNA_SeqFile,SamFolder,
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
         Constraint("format_type",MUST_BE,'samfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("ref",SET_TO_ONE_OF,['human','mouse']),
         help=("extract rna files with different format")
        ),
    Module(
        'extract_rna_files_bam',
         RNA_SeqFile,BamFolder,
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
         Constraint("format_type",MUST_BE,'bamfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("ref",SET_TO_ONE_OF,['human','mouse']),
         help=("extract rna files with bam format")
        ),
    Module(
        'extract_rna_files_fastq',
         RNA_SeqFile,FastqFolder,
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
         Constraint("format_type",MUST_BE,'fastqfolder'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help=("extract rna files with fa or fastq format")
        ),
    Module(
        'normalize_with_rsem_sam',
         [SamFolder,SampleGroupFile],
         SignalFile_rule.SignalFile_Postprocess,
         
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Consequence("preprocess",SET_TO,"rsem"),
         Consequence("logged",SET_TO,"unknown"),
         Consequence("predataset",SET_TO,"no"),
         Consequence("format",SET_TO,"tdf"),
         help=("process SamFolder , generate SignalFile_Postprocess with preprocess rsem")
        ),
     Module(
        'normalize_with_rsem_bam',
         [BamFolder,SampleGroupFile],
         SignalFile_rule.SignalFile_Postprocess,
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Consequence("preprocess",SET_TO,"rsem"),
         Consequence("logged",SET_TO,"unknown"),
         Consequence("predataset",SET_TO,"no"),
         Consequence("format",SET_TO,"tdf"),
         help=("process BamFolder , generate SignalFile_Postprocess with preprocess rsem")
        ),
    Module(
        'normalize_with_rsem_fastq',
         [FastqFolder,SampleGroupFile],
         SignalFile_rule.SignalFile_Postprocess,
         OptionDef("fastq_ref",'human',help='ref for fastq file to align,human or mouse'),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Consequence("preprocess",SET_TO,"rsem"),
         Consequence("logged",SET_TO,"unknown"),
         Consequence("predataset",SET_TO,"no"),
         Consequence("format",SET_TO,"tdf"),
         help=("process FastqFolder , generate SignalFile_Postprocess with preprocess rsem")
        ),
    ]
