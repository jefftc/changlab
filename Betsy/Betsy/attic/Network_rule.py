# ReportFile
#from Betsy.bie3 import *
#from Betsy.bie_rules import GeneExpProcessing,ArrayPlatforms,Clustering,Classification,GOAnalysis,PcaAnalysis,BasicDataTypes
#from Betsy.bie_rules import GeneExpSignature, DiffExp,GSEAAnalysis,ClinicalOutcomes,Heatmap

from Betsy.bie3 import DataType, ModuleNode, Constraint
from Betsy.bie3 import MUST_BE

from Betsy.rules import GeneExpSignature
from Betsy.rules import GeneExpProcessing
from Betsy.rules import ArrayPlatforms
from Betsy.rules import ClinicalOutcomes
from Betsy.rules import Clustering

from Betsy.rules import DiffExp
from Betsy.rules import Classification

from Betsy.rules import GOAnalysis
from Betsy.rules import GSEAAnalysis
from Betsy.rules import PcaAnalysis

NetworkFile = DataType(
    'NetworkFile',
    help="network file"
    )
NetworkFile_Test = DataType(
    'NetworkFile_Test',
    help="network file"
    )
all_data_types = [
    NetworkFile,
    NetworkFile_Test,
    ]

                                  
all_modules = [
    ModuleNode(
        'group_all_module_2minutes',
        [
            GeneExpProcessing.SignalFile,#0
            GeneExpProcessing.SignalFile,#1
            GeneExpProcessing.SignalFile,#2
            GeneExpProcessing.SignalFile,#3
            GeneExpProcessing.SignalFile,#4
            GeneExpProcessing.SignalFile,#5
            ArrayPlatforms.ControlFile,#6
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
        Constraint('gene_order', MUST_BE, 'ttest_p', 2),
        Constraint('gene_order', MUST_BE, 'ttest_fdr', 3),  # WAS diff_ttest(?)
        Constraint('duplicate_probe',MUST_BE,'closest_probe',2),),
    ModuleNode(
        'group_all_module',
        [
            GeneExpProcessing.SignalFile,#0
            GeneExpProcessing.SignalFile,#1
            GeneExpProcessing.SignalFile,#2
            GeneExpProcessing.SignalFile,#3
            GeneExpProcessing.SignalFile,#4
            GeneExpProcessing.SignalFile,#5
            ArrayPlatforms.ControlFile,#6
            #---------------------------
            GeneExpProcessing.IntensityPlot,#7
            ArrayPlatforms.ActbPlot,#8
            ArrayPlatforms.Hyb_barPlot,#9
            ArrayPlatforms.ControlPlot,#10
            ArrayPlatforms.BiotinPlot,#11
            ArrayPlatforms.HousekeepingPlot,#12
            #--------------------------
            GOAnalysis.GatherFile,#13
            GSEAAnalysis.GseaFile,#14
            GeneExpSignature.GenesetPlot,#15
            GeneExpSignature.GenesetPlot,#16
            GOAnalysis.DavidFile,#17
            ClinicalOutcomes.EMTAnalysis,#18
            DiffExp.DiffExprFile,#19
            DiffExp.DiffExprFile,#20
            DiffExp.DiffExprFile,#21
            ClinicalOutcomes.ClinicalAnalysis,#22
            PcaAnalysis.PcaPlot,#23
            GeneExpSignature.SignatureScore,#24
            Clustering.Heatmap,#25
            Clustering.Heatmap,#26
            #Clustering.Heatmap,#27
            Clustering.Heatmap,#27
            Clustering.Heatmap,#28
            Classification.PredictionPCAPlot,#29
            Classification.PredictionPlot,#30
            Classification.PredictionPCAPlot,#31
            Classification.PredictionPlot,#32
            Classification.PredictionPCAPlot,#33
            Classification.PredictionPlot,#34
            GeneExpProcessing.SignalFile,#35 tcga
            ],
        NetworkFile,
        Constraint('preprocess', MUST_BE, 'mas5',0),
        Constraint('preprocess', MUST_BE, 'rma',1),
        Constraint('preprocess', MUST_BE, 'unknown',2),
        Constraint('preprocess', MUST_BE, 'illumina',3),
        Constraint('preprocess', MUST_BE, 'loess',4),
        Constraint('preprocess', MUST_BE, 'agilent',5),
        #Constraint('preprocess', MUST_BE, 'RSEM_exons',36),
        Constraint('preprocess', MUST_BE, 'tpm',35),
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
        #Constraint("contents", MUST_BE, "class0,class1", 13),
        Constraint("contents", MUST_BE, "class0,class1", 14),
        Constraint("contents", MUST_BE, "class0,class1", 15),
        Constraint("contents", MUST_BE, "class0,class1", 16),
        #Constraint("contents", MUST_BE, "class0,class1", 17),
        Constraint("contents", MUST_BE, "class0,class1", 18),
        #Constraint("contents", MUST_BE, "class0,class1", 19),
        #Constraint("contents", MUST_BE, "class0,class1", 20),
        #Constraint("contents", MUST_BE, "class0,class1", 21),
        Constraint("contents", MUST_BE, "class0,class1", 22),
        Constraint("contents", MUST_BE, "class0,class1", 23),
        Constraint("contents", MUST_BE, "class0,class1", 24),
        Constraint("contents", MUST_BE, "class0,class1", 25),
        Constraint("contents", MUST_BE, "class0,class1", 26),
        Constraint("contents", MUST_BE, "class0,class1", 27),
        Constraint("contents", MUST_BE, "class0,class1", 28),
        Constraint("contents", MUST_BE, "class0,class1", 35),
	    Constraint('quantile_norm', MUST_BE, 'yes', 0),
        Constraint('gene_center', MUST_BE, 'mean', 0),
        Constraint('gene_normalize', MUST_BE, 'variance', 0),
        Constraint('gene_order', MUST_BE, 'class_neighbors', 0),
        Constraint('annotate', MUST_BE, 'yes', 0),
        Constraint('rename_sample', MUST_BE, 'yes', 0),
        Constraint('platform', MUST_BE, 'yes', 0),
        Constraint('rename_sample', MUST_BE, 'yes', 0),
        Constraint('num_features', MUST_BE, 'yes', 0),
        Constraint('unique_genes', MUST_BE, 'high_var', 0),
        Constraint('duplicate_probe', MUST_BE, 'high_var_probe', 0),
        Constraint('group_fc', MUST_BE, 'yes', 0),
        Constraint('filter_missing_values', MUST_BE, 'yes', 0),
        #bie3.Attribute(rulebase.SignalFile, "gene_center", "mean"),
        #bie3.Attribute(rulebase.SignalFile, "gene_normalize", "variance"),
        #bie3.Attribute(rulebase.SignalFile,"gene_order",'class_neighbors'),
        #bie3.Attribute(rulebase.SignalFile,"annotate",'yes'),
        #bie3.Attribute(rulebase.SignalFile,"rename_sample",'yes'),
        #bie3.Attribute(rulebase.SignalFile,"platform",'yes'),
        #bie3.Attribute(rulebase.SignalFile,"rename_sample",'yes'),
        #bie3.Attribute(rulebase.SignalFile,"num_features",'yes'),
        #bie3.Attribute(rulebase.SignalFile,"unique_genes",'high_var'),
        #bie3.Attribute(rulebase.SignalFile,"duplicate_probe",'high_var_probe'),
        #bie3.Attribute(rulebase.SignalFile,"group_fc",'yes'),
        #bie3.Attribute(rulebase.SignalFile,"filter",'yes'),
        
        Constraint('quantile_norm', MUST_BE, 'yes', 1),
        Constraint('quantile_norm', MUST_BE, 'yes', 2),
        Constraint('quantile_norm', MUST_BE, 'yes', 3),
        Constraint('quantile_norm', MUST_BE, 'yes', 4),
        Constraint('quantile_norm', MUST_BE, 'yes', 5),
        Constraint('quantile_norm', MUST_BE, 'yes', 35),
        Constraint('combat_norm', MUST_BE, 'yes', 0),
        Constraint('filter_and_threshold', MUST_BE, 'yes', 0),
        Constraint('logged', MUST_BE, 'no', 0),
        Constraint('shiftscale_norm', MUST_BE, 'yes', 1),
        Constraint('dwd_norm', MUST_BE, 'yes', 2),
        Constraint('bfrm_norm', MUST_BE, 'yes', 3),
        Constraint('missing_algorithm', MUST_BE, 'median_fill', 35),
        Constraint('gene_order', MUST_BE, 'ttest_p', 2),
        Constraint('gene_order', MUST_BE, 'ttest_fdr', 3), # was diff_ttest
        Constraint('platform', MUST_BE, 'yes', 2),
        Constraint('duplicate_probe',MUST_BE,'closest_probe',2),
        Constraint("allgenes",MUST_BE,'yes',15),
        Constraint("allgenes",MUST_BE,'no',16),
        #Constraint("gene_order",MUST_BE,'diff_sam',19),
        #Constraint("gene_order",MUST_BE,'diff_ebayes',20),
        #Constraint("gene_order",MUST_BE,'diff_fold_change',21),
        Constraint("cluster_alg",MUST_BE,'som',26),
        #Constraint("cluster_alg",MUST_BE,'pca',27),
        Constraint("cluster_alg",MUST_BE,'kmeans',27),
        Constraint("cluster_alg",MUST_BE,'hierarchical',28),
        Constraint("classify_alg",MUST_BE,'svm',30),#30
        Constraint("actual_label",MUST_BE,'yes',30),#30
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
        #Constraint('annotate', MUST_BE, 'yes', 1),
        #Constraint('annotate', MUST_BE, 'yes', 2),
        #Constraint('annotate', MUST_BE, 'yes', 3),
        #Constraint('annotate', MUST_BE, 'yes', 4),
        #Constraint('annotate', MUST_BE, 'yes', 5),
        ),
    ]
