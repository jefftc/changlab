from Betsy import bie
import SignalFile_rule
import SignalFile2_rule

GseaFile = bie.DataType(
    'GseaFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    )
list_files = [GseaFile]
all_modules=[
    bie.Module(
        'annotate_genes_with_gsea',
        [SignalFile_rule.ClassLabelFile(cls_format='cls'),
         SignalFile2_rule.SignalFile2(format='gct',logged='yes')],
        GseaFile)]
    


    


