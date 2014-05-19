# ReportFile
from Betsy.bie3 import *
from Betsy.bie_rules import SignalFile_rule,ClusterFile_rule,ClassifyFile_rule,PcaAnalysis_rule,plot_rule
from Betsy.bie_rules import GenesetAnalysis_rule, DiffExprFile_rule,GatherFile_rule,GseaFile_rule

ReportFile = DataType(
    'ReportFile',
    AttributeDef("report_type",['normalize_file','batch_effect_remove',
                            'classify','cluster','diffgenes',
                            'heatmap','geneset','all'],'normalize_file','normalize_file')
    )

list_files = [ReportFile]

                                  
all_modules = [
    Module(
        'make_normalize_report',
         [
         SignalFile_rule.SignalFile,
         plot_rule.IntensityPlot,
         plot_rule.ControlPlot,
         PcaAnalysis_rule.PcaPlot,
         plot_rule.ActbPlot,
         PcaAnalysis_rule.PcaPlot],
          
         ReportFile,
         Constraint('preprocess',CAN_BE_ANY_OF,['mas5','agilent','loess','unknown'],0),
         Constraint("annotate",MUST_BE,"yes",0),
         Constraint("contents",CAN_BE_ANY_OF,["train0","train1", "test",
                             "class0,class1,test","class0",
                             "class1", "class0,class1","unspecified"],0),
         Constraint('quantile_norm',MUST_BE,'no',5),
         Constraint('combat_norm',MUST_BE,'no',5),
         Constraint('shiftscale_norm',MUST_BE,'no',5),
         Constraint('bfrm_norm',MUST_BE,'no',5),
         Constraint('dwd_norm',MUST_BE,'no',5),
         Constraint('gene_center',MUST_BE,'no',5),
         Constraint('gene_normalize',MUST_BE,'no',5),
         Constraint('unique_genes',MUST_BE,'no',5),
         Constraint('platform',MUST_BE,'no',5),
         Constraint('group_fc',MUST_BE,'no',5),
         Constraint('num_features',MUST_BE,'no',5),
        
##         Constraint('quantile_norm',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('combat_norm',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('shiftscale_norm',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('bfrm_norm',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('dwd_norm',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('gene_center',CAN_BE_ANY_OF,['median','mean','no'],0),#
##         Constraint('gene_normalize',CAN_BE_ANY_OF,['variance','sum_of_squares','no'],0),#
##         Constraint('unique_genes',CAN_BE_ANY_OF,['average_genes', 'high_var', 'first_gene'],0),#
##         Constraint('platform',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('group_fc',CAN_BE_ANY_OF,['yes','no'],0),#
##         Constraint('num_features',CAN_BE_ANY_OF,['yes','no'],0),#
##        
##         Constraint('quantile_norm',SAME_AS,0,3),#
##         Constraint('combat_norm',SAME_AS,0,3),#
##         Constraint('shiftscale_norm',SAME_AS,0,3),#
##         Constraint('bfrm_norm',SAME_AS,0,3),#
##         Constraint('dwd_norm',SAME_AS,0,3),#
##         Constraint('gene_center',SAME_AS,0,3),#
##         Constraint('gene_normalize',SAME_AS,0,3),#
##         Constraint('unique_genes',SAME_AS,0,3),#
##         Constraint('platform',SAME_AS,0,3),#
##         Constraint('group_fc',SAME_AS,0,3),#
##         Constraint('num_features',SAME_AS,0,3),#
         
         Constraint('contents',SAME_AS,0,1),
         Constraint('contents',SAME_AS,0,2),
         Constraint('contents',SAME_AS,0,3),
         Constraint('contents',SAME_AS,0,4),
         Constraint('contents',SAME_AS,0,5),
         Constraint("preprocess",SAME_AS,0,3),
         Constraint("preprocess",SAME_AS,0,4),
         Constraint("preprocess",SAME_AS,0,5),
         Consequence('report_type',SET_TO,'normalize_file'),
        ),
    
    Module(
        'make_normalize_report',
         [SignalFile_rule.SignalFile,
         plot_rule.IntensityPlot,
         plot_rule.ControlPlot,
         PcaAnalysis_rule.PcaPlot,
         plot_rule.ActbPlot,
         PcaAnalysis_rule.PcaPlot],
          
         ReportFile,
         Constraint('preprocess',MUST_BE,'rma',0),
         Constraint("annotate",MUST_BE,"yes",0),
         Constraint("contents",CAN_BE_ANY_OF,["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],0),
         Constraint('quantile_norm',MUST_BE,'yes',5),
         Constraint('quantile_norm',MUST_BE,'yes',0),
        Constraint('quantile_norm',MUST_BE,'yes',3),
         Constraint('combat_norm',MUST_BE,'no',5),
         Constraint('shiftscale_norm',MUST_BE,'no',5),
         Constraint('bfrm_norm',MUST_BE,'no',5),
         Constraint('dwd_norm',MUST_BE,'no',5),
         Constraint('gene_center',MUST_BE,'no',5),
         Constraint('gene_normalize',MUST_BE,'no',5),
         Constraint('unique_genes',MUST_BE,'no',5),
         Constraint('platform',MUST_BE,'no',5),
         Constraint('group_fc',MUST_BE,'no',5),
         Constraint('num_features',MUST_BE,'no',5),
         Constraint('contents',SAME_AS,0,1),
         Constraint('contents',SAME_AS,0,2),
         Constraint('contents',SAME_AS,0,3),
         Constraint('contents',SAME_AS,0,4),
         Constraint('contents',SAME_AS,0,5),
         Constraint("preprocess",SAME_AS,0,3),
         Constraint("preprocess",SAME_AS,0,4),
         Constraint("preprocess",SAME_AS,0,5),
         Consequence('report_type',SET_TO,'normalize_file'),
        ),
                                      
    Module(
        'make_normalize_report_illumina',
       [ SignalFile_rule.SignalFile,
         PcaAnalysis_rule.PcaPlot,
         PcaAnalysis_rule.PcaPlot,
         plot_rule.IntensityPlot,
         plot_rule.ActbPlot,
         plot_rule.BiotinPlot,
         plot_rule.HousekeepingPlot,
         plot_rule.Hyb_barPlot,
         SignalFile_rule.ControlFile],ReportFile,
         Constraint("contents",CAN_BE_ANY_OF,["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],0),
         Constraint('preprocess',MUST_BE,'illumina',0),
         Constraint('annotate',MUST_BE,'yes',0),
         Constraint('preprocess',MUST_BE,'illumina',1),
         Constraint('quantile_norm',MUST_BE,'no',1),
         Constraint('combat_norm',MUST_BE,'no',1),
         Constraint('shiftscale_norm',MUST_BE,'no',1),
         Constraint('bfrm_norm',MUST_BE,'no',1),
         Constraint('dwd_norm',MUST_BE,'no',1),
         Constraint('gene_center',MUST_BE,'no',1),
         Constraint('gene_normalize',MUST_BE,'no',1),
         Constraint('unique_genes',MUST_BE,'no',1),
         Constraint('platform',MUST_BE,'no',1),
         Constraint('group_fc',MUST_BE,'no',1),
         Constraint('num_features',MUST_BE,'no',1),
         Constraint('preprocess',MUST_BE,'illumina',2),
         Constraint('preprocess',MUST_BE,'illumina',8),
         Constraint('format',MUST_BE,'gct',8),
         Constraint("logged",MUST_BE,'no',8),
         Constraint('contents',SAME_AS,0,1),
         Constraint('contents',SAME_AS,0,2),
         Constraint('contents',SAME_AS,0,3),
         Constraint('contents',SAME_AS,0,4),
         Constraint('contents',SAME_AS,0,5),
         Constraint('contents',SAME_AS,0,6),
         Constraint('contents',SAME_AS,0,7),
         Constraint('contents',SAME_AS,0,8),
         Constraint("preprocess",SAME_AS,0,1),
         Constraint("preprocess",SAME_AS,0,2),
         Constraint("preprocess",SAME_AS,0,4),
         Consequence('report_type',SET_TO,'normalize_file')),
    
    Module(
        'make_cluster_report',
        [ClusterFile_rule.ClusterFile,
         ClusterFile_rule.Heatmap],ReportFile,
         Constraint("cluster_alg",CAN_BE_ANY_OF,['som','pca','kmeans','hierarchical'],0),
         Constraint("distance",CAN_BE_ANY_OF,['correlation','euclidean'],0),
         Constraint("cluster_alg",SAME_AS,0,1),
         Constraint("distance",SAME_AS,0,1),
         Consequence('report_type',SET_TO,'cluster'),),
    Module(
        'make_classify_report',
        [SignalFile_rule.SignalFile,
         ClassifyFile_rule.ClassifyFile,
         ClassifyFile_rule.ClassifyFile,
         ClassifyFile_rule.PredictionPlot,
         ClassifyFile_rule.PredictionPlot,
         ClassifyFile_rule.PredictionPCAPlot,
         ClassifyFile_rule.ClassifyFile,
         ClassifyFile_rule.ClassifyFile,
         ClassifyFile_rule.PredictionPlot,
         ClassifyFile_rule.PredictionPlot,
         ClassifyFile_rule.PredictionPCAPlot
         ],ReportFile,
         Constraint('contents',MUST_BE,'class0,class1,test',0),
         Constraint('logged',MUST_BE,'yes',0),
         Constraint("format",MUST_BE,'gct',0),
         Constraint("classify_alg",MUST_BE,'svm',1),
         Constraint("actual_label",MUST_BE,'yes',1),
         Constraint("classify_alg",MUST_BE,'svm',2),
         Constraint("actual_label",MUST_BE,'no',2),
         Constraint("loocv",MUST_BE,'yes',2),
         Constraint("classify_alg",MUST_BE,'svm',3),
         Constraint("actual_label",MUST_BE,'yes',3),
         Constraint("loocv",MUST_BE,'no',3),
         Constraint("classify_alg",MUST_BE,'svm',4),
         Constraint("actual_label",MUST_BE,'no',4),
         Constraint("loocv",MUST_BE,'yes',4),
         Constraint("classify_alg",MUST_BE,'svm',5),
         Constraint("actual_label",MUST_BE,'yes',5),
         Constraint("loocv",MUST_BE,'no',5),
         Constraint("classify_alg",MUST_BE,'weighted_voting',6),
         Constraint("actual_label",MUST_BE,'yes',6),
         Constraint("loocv",MUST_BE,'no',6),
         Constraint("classify_alg",MUST_BE,'weighted_voting',7),
         Constraint("actual_label",MUST_BE,'no',7),
         Constraint("loocv",MUST_BE,'yes',7),
         Constraint("classify_alg",MUST_BE,'weighted_voting',8),
         Constraint("actual_label",MUST_BE,'yes',8),
         Constraint("loocv",MUST_BE,'no',8),
         Constraint("classify_alg",MUST_BE,'weighted_voting',9),
         Constraint("actual_label",MUST_BE,'no',9),
         Constraint("loocv",MUST_BE,'yes',9),
         Constraint("classify_alg",MUST_BE,'weighted_voting',10),
         Constraint("actual_label",MUST_BE,'yes',10),
         Constraint("loocv",MUST_BE,'no',10),
         Consequence('report_type',SET_TO,'classify')),
    Module(
        'make_heatmap_report',
        ClusterFile_rule.Heatmap,ReportFile,
        Constraint("cluster_alg",MUST_BE,'no_cluster_alg'),
        Consequence("report_type",SET_TO,'heatmap')),
    Module(
        'make_geneset_report',
        [GenesetAnalysis_rule.GenesetAnalysis,
         GenesetAnalysis_rule.GenesetPlot],ReportFile,
         Consequence("report_type",SET_TO,'geneset')),
        
    Module(
        'make_diffgenes_report',
        [DiffExprFile_rule.DiffExprFile,DiffExprFile_rule.DiffExprFile,
         ClusterFile_rule.Heatmap,GatherFile_rule.GatherFile,
         GseaFile_rule.GseaFile],ReportFile,
         UserInputDef("hm_width",50),
         UserInputDef("hm_height",1),
         Constraint("diff_expr",MUST_BE,'t_test',0),
         Constraint("diff_expr",MUST_BE,'sam',1),
         Constraint("cluster_alg",MUST_BE,'no_cluster_alg',2),
         #Constraint("hm_width",MUST_BE,"yes",2),
         #Constraint("hm_height",MUST_BE,"yes",2),
         Consequence("report_type",SET_TO,'diffgenes')),
    Module(
        'make_batch_effect_report',
        [SignalFile_rule.SignalFile,SignalFile_rule.SignalFile_Merge,
         SignalFile_rule.SignalFile,SignalFile_rule.SignalFile_Merge,
         SignalFile_rule.SignalFile_Merge],ReportFile,
         Constraint("quantile_norm",MUST_BE,'yes',0),
         Constraint("dwd_norm",MUST_BE,'no',0),
         Constraint("bfrm_norm",MUST_BE,'no',0),
         Constraint("combat_norm",MUST_BE,'no',0),
         Constraint("shiftscale_norm",MUST_BE,'no',0),
         Constraint("quantile_norm",MUST_BE,'yes',1),
         Constraint("dwd_norm",MUST_BE,'yes',1),
         Constraint("bfrm_norm",MUST_BE,'no',1),
         Constraint("combat_norm",MUST_BE,'no',1),
         Constraint("shiftscale_norm",MUST_BE,'no',1),
         Constraint("quantile_norm",MUST_BE,'yes',2),
         Constraint("dwd_norm",MUST_BE,'no',2),
         Constraint("bfrm_norm",MUST_BE,'no',2),
         Constraint("combat_norm",MUST_BE,'no',2),
         Constraint("shiftscale_norm",MUST_BE,'yes',2),
         Constraint("quantile_norm",MUST_BE,'yes',3),
         Constraint("dwd_norm",MUST_BE,'no',3),
         Constraint("bfrm_norm",MUST_BE,'yes',3),
         Constraint("combat_norm",MUST_BE,'no',3),
         Constraint("shiftscale_norm",MUST_BE,'no',3),
         Constraint("quantile_norm",MUST_BE,'yes',4),
         Constraint("dwd_norm",MUST_BE,'no',4),
         Constraint("bfrm_norm",MUST_BE,'no',4),
         Constraint("combat_norm",MUST_BE,'yes',4),
         Constraint("shiftscale_norm",MUST_BE,'no',4),
         Consequence("report_type",SET_TO,'batch_effect_remove'))

        ]
  
        
