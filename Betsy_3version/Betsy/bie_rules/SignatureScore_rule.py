#SignatureScore
from Betsy.bie3 import *
import SignalFile_rule
SignatureScore = DataType(
    'SignatureScore')
list_files = [SignatureScore]
all_modules = [
    Module(
        'score_pathway_with_scoresig',
        [SignalFile_rule.PrettySignalFile,SignalFile_rule.PrettySignalFile],SignatureScore,
        UserInputDef('platform_value','HG_U133A'),
        Constraint("format",MUST_BE,'tdf',0),
        Constraint("preprocess",MUST_BE,'rma',0),
        Constraint("logged",MUST_BE,'yes',0),
        Constraint("missing_values",MUST_BE,'no',0),
        Constraint("platform",MUST_BE,"yes",0),
        Constraint("duplicate_probe",MUST_BE,'high_var_probe',0),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("preprocess",MUST_BE,'mas5',1),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("missing_values",MUST_BE,'no',1),
        Constraint("platform",MUST_BE,'yes',1),
        Constraint("duplicate_probe",MUST_BE,'high_var_probe',1)
        )]
