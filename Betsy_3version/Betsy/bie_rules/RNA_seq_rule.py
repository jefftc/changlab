#RNA_seq_rule
from Betsy.bie3 import *
import SignalFile_rule

RNA_SeqFile=DataType(
    'RNA_SeqFile',
    AttributeDef("rna_format",['unknown','fa_or_sam'],
                 "unknown","unknown",help='rna file format'),
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                 'unspecified','unspecified',help="contents"),
    help="RNA Seq File"
    )



list_files = [RNA_SeqFile]

all_modules = [
    Module(
        'extract_rna_files',
         RNA_SeqFile,RNA_SeqFile,
         Constraint("rna_format",MUST_BE,'unknown'),
         Consequence("rna_format",SET_TO,'fa_or_sam'),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help=("extract rna files with fa or sam format")
        ),
    Module(
        'normalize_with_rsem',
         RNA_SeqFile,SignalFile_rule.SignalFile_Postprocess,
         Constraint("rna_format",MUST_BE,'fa_or_sam'),
         Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("preprocess",SET_TO,"rsem"),
         Consequence("logged",SET_TO,"unknown"),
         Consequence("predataset",SET_TO,"no"),
         Consequence("format",SET_TO,"tdf"),
         help=("process RNA_SeqFile , generate SignalFile_Postprocess with preprocess rsem")
        ),
    
    ]
