#DiffExprFile_rule
from Betsy.bie3 import *
import SignalFile_rule
import SignalFile2_rule

DiffExprFile=DataType(
    'DiffExprFile',
    Attribute("diff_expr",['t_test','sam'],"t_test",'t_test'),
    Attribute("sam_delta",["yes","no"],"no","no"),
    Attribute("sam_delta",["yes","no"],"no","no"))

list_files = [DiffExprFile]

all_modules = [
    Module(
        'calc_diffexp_with_ttest',
        [SignalFile_rule.ClassLabelFile,SignalFile2_rule.SignalFile2],DiffExprFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Consequence("diff_expr",SET_TO,'t_test')),
    
    Module(
        'calc_diffexp_with_sam',
        [SignalFile_rule.ClassLabelFile,SignalFile2_rule.SignalFile2],DiffExprFile,
        UserInput("sam_delta_value",0),
        UserInput("sam_foldchange_value",0),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Consequence("diff_expr",SET_TO,'sam')),
    ]
