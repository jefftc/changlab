# ReportFile
from Betsy.bie3 import *
from Betsy.bie_rules import SignalFile_rule,ClusterFile_rule,ClassifyFile_rule,PcaAnalysis_rule,plot_rule,DavidFile_rule
from Betsy.bie_rules import GenesetAnalysis_rule, DiffExprFile_rule,GatherFile_rule,GseaFile_rule,ReportFile_rule,EMTAnalysis_rule
from Betsy.bie_rules import ClinicalAnalysis_rule,SignatureScore_rule
NetworkFile = DataType(
    'NetworkFile',
    help="network file"
    )
NetworkFile_Test = DataType(
    'NetworkFile_Test',
    help="network file"
    )
list_files = [NetworkFile,NetworkFile_Test]

                                  
all_modules = [
##    Module('group_all',
##           [ReportFile_rule.ReportFile,
##            ReportFile_rule.ReportFile,
##            ReportFile_rule.ReportFile,
##            ReportFile_rule.ReportFile,
##            ReportFile_rule.ReportFile,
##            ReportFile_rule.ReportFile,
##            ReportFile_rule.ReportFile,
##            EMTAnalysis_rule.EMTAnalysis],
##            NetworkFile,
##            Constraint('report_type',MUST_BE,'classify',0),
##            Constraint('report_type',MUST_BE,'geneset',1),
##            Constraint('report_type',MUST_BE,'heatmap',2),
##            Constraint('report_type',MUST_BE,'diffgenes',3),
##            Constraint('report_type',MUST_BE,'cluster',4),
##            Constraint('report_type',MUST_BE,'batch_effect_remove',5),
##            Constraint('report_type',MUST_BE,'normalize_file',6),),
    Module('group_all_module_2mintues',
           [SignalFile_rule.SignalFile,#0
            SignalFile_rule.SignalFile,#1
            SignalFile_rule.SignalFile,#2
            SignalFile_rule.SignalFile,#3
            SignalFile_rule.SignalFile,#4
            SignalFile_rule.SignalFile,#5
            SignalFile_rule.ControlFile,#6
            ],
            NetworkFile_Test,
            Constraint('preprocess', MUST_BE, 'mas5',0),
            Constraint('preprocess', MUST_BE, 'rma',1),
            Constraint('preprocess', MUST_BE, 'unknown',2),
            Constraint('preprocess', MUST_BE, 'illumina',3),
            Constraint('preprocess', MUST_BE, 'loess',4),
            Constraint('preprocess', MUST_BE, 'agilent',5),
            Constraint("contents", MUST_BE, "class0,class1", 0),
            Constraint("contents", MUST_BE, "class0,class1", 1),
            Constraint("contents", MUST_BE, "class0,class1", 2),
            Constraint("contents", MUST_BE, "class0,class1", 3),
            Constraint("contents", MUST_BE, "class0,class1", 4),
            Constraint("contents", MUST_BE, "class0,class1", 5),
            Constraint("contents", MUST_BE, "class0,class1", 6),
           Constraint('quantile_norm', MUST_BE, 'yes', 0),
            Constraint('quantile_norm', MUST_BE, 'yes', 1),
            Constraint('quantile_norm', MUST_BE, 'yes', 2),
            Constraint('quantile_norm', MUST_BE, 'yes', 3),
            Constraint('quantile_norm', MUST_BE, 'yes', 4),
            Constraint('quantile_norm', MUST_BE, 'yes', 5),
            Constraint('combat_norm', MUST_BE, 'yes', 0),
            Constraint('shiftscale_norm', MUST_BE, 'yes', 1),
            Constraint('dwd_norm', MUST_BE, 'yes', 2),
            Constraint('bfrm_norm', MUST_BE, 'yes', 3),
            Constraint('missing_algorithm', MUST_BE, 'median_fill', 2),
            Constraint('gene_order', MUST_BE, 't_test_p', 2),
            Constraint('gene_order', MUST_BE, 'diff_ttest', 3),
            Constraint('duplicate_probe',MUST_BE,'closest_probe',2),),
    Module('group_all_module',
           [SignalFile_rule.SignalFile,#0
            SignalFile_rule.SignalFile,#1
            SignalFile_rule.SignalFile,#2
            SignalFile_rule.SignalFile,#3
            SignalFile_rule.SignalFile,#4
            SignalFile_rule.SignalFile,#5
            SignalFile_rule.ControlFile,#6
            #---------------------------
            plot_rule.IntensityPlot,#7
            plot_rule.ActbPlot,#8
            plot_rule.Hyb_barPlot,#9
            plot_rule.ControlPlot,#10
            plot_rule.BiotinPlot,#11
            plot_rule.HousekeepingPlot,#12
            #--------------------------
            GatherFile_rule.GatherFile,#13
            GseaFile_rule.GseaFile,#14
            GenesetAnalysis_rule.GenesetPlot,#15
            GenesetAnalysis_rule.GenesetPlot,#16
            DavidFile_rule.DavidFile,#17
            EMTAnalysis_rule.EMTAnalysis,#18
            DiffExprFile_rule.DiffExprFile,#19
            DiffExprFile_rule.DiffExprFile,#20
            DiffExprFile_rule.DiffExprFile,#21
            ClinicalAnalysis_rule.ClinicalAnalysis,#22
            PcaAnalysis_rule.PcaPlot,#23
            SignatureScore_rule.SignatureScore,#24
            ClusterFile_rule.Heatmap,#25
            ClusterFile_rule.Heatmap,#26
            ClusterFile_rule.Heatmap,#27
            ClusterFile_rule.Heatmap,#28
            ClusterFile_rule.Heatmap,#29
            ClassifyFile_rule.PredictionPCAPlot,#30
            ClassifyFile_rule.PredictionPlot,#31
            ClassifyFile_rule.PredictionPCAPlot,#32
            ClassifyFile_rule.PredictionPlot,#33
            ClassifyFile_rule.PredictionPCAPlot,#34
            ClassifyFile_rule.PredictionPlot,#35
            
            ],
            NetworkFile,
            Constraint('preprocess', MUST_BE, 'mas5',0),
            Constraint('preprocess', MUST_BE, 'rma',1),
            Constraint('preprocess', MUST_BE, 'unknown',2),
            Constraint('preprocess', MUST_BE, 'illumina',3),
            Constraint('preprocess', MUST_BE, 'loess',4),
            Constraint('preprocess', MUST_BE, 'agilent',5),
            Constraint("contents", MUST_BE, "class0,class1", 0),
            Constraint("contents", MUST_BE, "class0,class1", 1),
            Constraint("contents", MUST_BE, "class0,class1", 2),
            Constraint("contents", MUST_BE, "class0,class1", 3),
            Constraint("contents", MUST_BE, "class0,class1", 4),
            Constraint("contents", MUST_BE, "class0,class1", 5),
            Constraint("contents", MUST_BE, "class0,class1", 6),
            Constraint("contents", MUST_BE, "class0,class1", 7),
            Constraint("contents", MUST_BE, "class0,class1", 8),
            Constraint("contents", MUST_BE, "class0,class1", 9),
            Constraint("contents", MUST_BE, "class0,class1", 10),
            Constraint("contents", MUST_BE, "class0,class1", 11),
            Constraint("contents", MUST_BE, "class0,class1", 12),
            Constraint("contents", MUST_BE, "class0,class1", 13),
            Constraint("contents", MUST_BE, "class0,class1", 14),
            Constraint("contents", MUST_BE, "class0,class1", 15),
            Constraint("contents", MUST_BE, "class0,class1", 16),
            Constraint("contents", MUST_BE, "class0,class1", 17),
            Constraint("contents", MUST_BE, "class0,class1", 18),
            Constraint("contents", MUST_BE, "class0,class1", 19),
            Constraint("contents", MUST_BE, "class0,class1", 20),
            Constraint("contents", MUST_BE, "class0,class1", 21),
            Constraint("contents", MUST_BE, "class0,class1", 22),
            Constraint("contents", MUST_BE, "class0,class1", 23),
            Constraint("contents", MUST_BE, "class0,class1", 24),
            Constraint("contents", MUST_BE, "class0,class1", 25),
            Constraint("contents", MUST_BE, "class0,class1", 26),
            Constraint("contents", MUST_BE, "class0,class1", 27),
            Constraint("contents", MUST_BE, "class0,class1", 28),
            Constraint("contents", MUST_BE, "class0,class1", 29),
           
            Constraint('quantile_norm', MUST_BE, 'yes', 0),
            Constraint('quantile_norm', MUST_BE, 'yes', 1),
            Constraint('quantile_norm', MUST_BE, 'yes', 2),
            Constraint('quantile_norm', MUST_BE, 'yes', 3),
            Constraint('quantile_norm', MUST_BE, 'yes', 4),
            Constraint('quantile_norm', MUST_BE, 'yes', 5),
            Constraint('combat_norm', MUST_BE, 'yes', 0),
            Constraint('shiftscale_norm', MUST_BE, 'yes', 1),
            Constraint('dwd_norm', MUST_BE, 'yes', 2),
            Constraint('bfrm_norm', MUST_BE, 'yes', 3),
            Constraint('missing_algorithm', MUST_BE, 'median_fill', 2),
            Constraint('gene_order', MUST_BE, 't_test_p', 2),
            Constraint('gene_order', MUST_BE, 'diff_ttest', 3),
            Constraint('duplicate_probe',MUST_BE,'closest_probe',2),
            Constraint("allgenes",MUST_BE,'yes',15),
            Constraint("allgenes",MUST_BE,'no',16),
            Constraint("gene_order",MUST_BE,'diff_sam',19),
            Constraint("gene_order",MUST_BE,'diff_ebayes',20),
            Constraint("gene_order",MUST_BE,'diff_fold_change',21),
            Constraint("cluster_alg",MUST_BE,'som',26),
            Constraint("cluster_alg",MUST_BE,'pca',27),
            Constraint("cluster_alg",MUST_BE,'kmeans',28),
            Constraint("cluster_alg",MUST_BE,'hierarchical',29),
            Constraint("classify_alg",MUST_BE,'svm',30),
            Constraint("actual_label",MUST_BE,'yes',30),
            Constraint("classify_alg",MUST_BE,'svm',31),
            Constraint("loocv",MUST_BE,'yes',31),
            Constraint("classify_alg",MUST_BE,'weighted_voting',32),
            Constraint("actual_label",MUST_BE,'yes',32),
            Constraint("classify_alg",MUST_BE,'weighted_voting',33),
            Constraint("loocv",MUST_BE,'yes',33),
            Constraint("classify_alg",MUST_BE,'random_forest',34),
            Constraint("actual_label",MUST_BE,'yes',34),
            Constraint("classify_alg",MUST_BE,'random_forest',35),
            Constraint("loocv",MUST_BE,'yes',35),
            
        
           ),
    
        ]
  
        
