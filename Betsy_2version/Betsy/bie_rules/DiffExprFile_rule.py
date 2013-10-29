#DiffExprFile_rule
from Betsy import bie
import SignalFile_rule
import SignalFile2_rule

DiffExprFile=bie.DataType(
    'DiffExprFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(diff_expr=['t_test','sam'],DEFAULT='t_test'),
    bie.Attribute(sam_delta=bie.ANYATOM,DEFAULT='0'),
    bie.Attribute(sam_foldchange=bie.ANYATOM,DEFAULT='0'))

list_files = [DiffExprFile]

all_modules = [
    bie.Module(
        'calc_diffexp_with_ttest',
        [SignalFile_rule.ClassLabelFile(cls_format='cls'),
         SignalFile2_rule.SignalFile2(logged='yes',format='tdf',gene_order='no')],
        DiffExprFile(diff_expr='t_test')),
    
    bie.Module(
        'calc_diffexp_with_sam',
        [SignalFile_rule.ClassLabelFile(cls_format='cls'),
         SignalFile2_rule.SignalFile2(logged='yes',format='tdf',gene_order='no')],
        DiffExprFile(diff_expr='sam',
                    sam_delta=bie.ANYATOM,sam_foldchange=bie.ANYATOM))
    ]
