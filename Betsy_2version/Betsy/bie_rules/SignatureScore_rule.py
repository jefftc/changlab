#SignatureScore
import bie
import SignalFile2_rule
SignatureScore = bie.DataType(
    'SignatureScore',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True))
list_files = [SignatureScore]
all_modules = [
    bie.Module(
        'score_pathway_with_scoresig',
        [SignalFile2_rule.SignalFile2(format='tdf',preprocess='rma',logged='yes',missing_values='no',
                     platform='HG_U133A',duplicate_probe='high_var_probe'),
         SignalFile2_rule.SignalFile2(format='tdf',preprocess='mas5',logged='yes',missing_values='no',
                     platform='HG_U133A',duplicate_probe='high_var_probe')],
        SignatureScore)]
