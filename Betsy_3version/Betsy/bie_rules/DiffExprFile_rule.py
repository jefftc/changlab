#DiffExprFile_rule
from Betsy.bie3 import *
import SignalFile_rule

DiffExprFile=DataType(
    'DiffExprFile',
    AttributeDef("diff_expr",['t_test','sam'],"t_test",'t_test'),
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified'))

list_files = [DiffExprFile]

all_modules = [
    Module(
        'calc_diffexp_with_ttest',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],DiffExprFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("diff_expr",SET_TO,'t_test'),
        Consequence('contents',SAME_AS_CONSTRAINT,0)),
    
    Module(
        'calc_diffexp_with_sam',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],DiffExprFile,
        UserInputDef("sam_delta_value",0),
        UserInputDef("sam_foldchange_value",0),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("diff_expr",SET_TO,'sam'),
        Consequence('contents',SAME_AS_CONSTRAINT,0)),
    ]
