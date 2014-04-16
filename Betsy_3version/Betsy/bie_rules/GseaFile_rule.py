from Betsy.bie3 import *
import SignalFile_rule

GseaFile = DataType('GseaFile',
                    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified'))
list_files = [GseaFile]
all_modules=[
    Module(
        'annotate_genes_with_gsea',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],GseaFile,
         Constraint("cls_format",MUST_BE,'cls',0),
         Constraint("contents",CAN_BE_ANY_OF,["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],0),
         Constraint("format",MUST_BE,'gct',1),
         Constraint("logged",MUST_BE,'yes',1),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0))]
    


    


