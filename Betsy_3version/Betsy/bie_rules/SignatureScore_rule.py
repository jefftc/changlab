#SignatureScore
from Betsy.bie3 import *
import SignalFile_rule
SignatureScore = DataType(
    'SignatureScore',AttributeDef("contents",SignalFile_rule.CONTENTS,
                                  'unspecified','unspecified',help="contents"),
    help="Signature score file")
list_files = [SignatureScore]
all_modules = [
    Module(
        'score_pathway_with_scoresig',
        [SignalFile_rule.SignalFile,SignalFile_rule.SignalFile],SignatureScore,
        OptionDef('platform_value','HG_U133A',help="platform to add"),
        Constraint("format",MUST_BE,'tdf',0),
        Constraint("preprocess",MUST_BE,'rma',0),
        Constraint("quantile_norm",MUST_BE,'yes',0),
        Constraint("logged",MUST_BE,'yes',0),
        Constraint("platform",MUST_BE,"yes",0),
        Constraint("duplicate_probe",MUST_BE,'high_var_probe',0),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("preprocess",MUST_BE,'mas5',1),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("platform",MUST_BE,'yes',1),
        Constraint("duplicate_probe",MUST_BE,'high_var_probe',1),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Consequence('contents',SAME_AS_CONSTRAINT,0),
        help="score pathway iwth scoresig method"
        )]
