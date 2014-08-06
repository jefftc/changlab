#DiffExprFile_rule
from Betsy.bie3 import *
import SignalFile_rule

DiffExprFile=DataType(
    'DiffExprFile',
    AttributeDef("gene_order",['diff_ttest','diff_sam','diff_ebayes','diff_fold_change'],
                 "diff_ttest",'diff_ttest',help="differential method"),
    AttributeDef("contents",SignalFile_rule.CONTENTS,'diff_unspecified',
                 'diff_unspecified',help='contents'),
    help="Differential expression result file")

list_files = [DiffExprFile]

all_modules = [
    Module(
        'calc_diffexp_with_ttest',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],DiffExprFile,
        UserInputDef("diffexp_foldchange_value",0,
                     help="fold change value for differential expression analysis"),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("gene_order",SET_TO,'diff_ttest'),
        Consequence('contents',SAME_AS_CONSTRAINT,0),
        help="calculate the differential expression with ttest method"),
    
    Module(
        'calc_diffexp_with_sam',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],DiffExprFile,
        UserInputDef("sam_delta_value",1.0,
                     help="delta value for sam differential expression method"),
        UserInputDef("diffexp_foldchange_value",0,
                     help="fold change value for differential expression analysis"),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("gene_order",SET_TO,'diff_sam'),
        Consequence('contents',SAME_AS_CONSTRAINT,0),
        help="calculate the differential expression with sam method"),
    
    Module(
        'calc_diffexp_with_ebayes',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],DiffExprFile,
        UserInputDef("diffexp_foldchange_value",0,
                     help="fold change value for differential expression analysis"),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("gene_order",SET_TO,'diff_ebayes'),
        Consequence('contents',SAME_AS_CONSTRAINT,0),
        help="calculate the differential expression with ebayes method"),
    Module(
        'calc_diffexp_with_fold_change',
        [SignalFile_rule.ClassLabelFile,SignalFile_rule.SignalFile],DiffExprFile,
        UserInputDef("diffexp_foldchange_value",0,
                     help="fold change value for differential expression analysis"),
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS,0),
        Constraint("logged",MUST_BE,'yes',1),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("gene_order",MUST_BE,'no',1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("gene_order",SET_TO,'diff_fold_change'),
        Consequence('contents',SAME_AS_CONSTRAINT,0),
        help="calculate the differential expression with fold change method"),
    Module(
        'generate_genelist_from_diffexprfile',
        DiffExprFile,SignalFile_rule.GeneListFile,
        UserInputDef("select_gene_by_foldchange","",
                     help="select gene by foldchange value,etc.5"),
        UserInputDef("select_gene_by_p_value","",
                     help="select gene by p-value value, etc. 0.05"),
        UserInputDef("select_gene_by_fdr","",
                     help="select gene by fdr value, etc.0.05"),
        Constraint("gene_order",CAN_BE_ANY_OF,['diff_ttest',
                                               'diff_sam','diff_ebayes','diff_fold_change']),
        Constraint(
            "contents", CAN_BE_ANY_OF, SignalFile_rule.CONTENTS),
        Consequence("gene_order",SAME_AS_CONSTRAINT),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="generate genelist from the result file of differential expression analysis"
       ),
    ]
