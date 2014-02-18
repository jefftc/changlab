from Betsy.bie3 import *
import SignalFile_rule

GseaFile = DataType('GseaFile')
list_files = [GseaFile]
all_modules=[
    Module(
        'annotate_genes_with_gsea',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.PrettySignalFile],GseaFile,
         Constraint("cls_format",MUST_BE,'cls',0),
         Constraint("format",MUST_BE,'gct',1),
         Constraint("logged",MUST_BE,'yes',1))]
    


    


