# ReportFile
from Betsy import bie
from Betsy.bie_rules import SignalFile_rule,SignalFile2_rule,ClusterFile_rule,ClassifyFile_rule,PcaAnalysis_rule
from Betsy.bie_rules import GenesetAnalysis_rule, DiffExprFile_rule,GatherFile_rule,GseaFile_rule

ReportFile = bie.DataType(
    'ReportFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(report_type=['normalize','batch_effect_remove',
                            'classify','cluster','diffgenes',
                            'heatmap','geneset','all'],DEFAULT='normalize'),
    bie.Attribute(preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"],
                  DEFAULT='unknown')
    )

list_files = [ReportFile]

#contents=["train0", "train1", "test", "class0,class1,test",
#                  "class0", "class1", "class0,class1",
#                  "no"],preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
#                  logged=['yes','no'],
#                   quantile_norm=['yes','no'],bfrm_norm=['yes','no'],combat_norm=['yes','no'],
#                   shiftscale_norm=['yes','no'],dwd_norm=['yes','no'],gene_center=['mean','median','no','unknown'],
#                   gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
#                   missing_algorithm=["none", "median_fill", "zero_fill"],
#                   unique_genes=["no", "average_genes", "high_var", "first_gene"],
#                   duplicate_probe=["no", "yes", "closest_probe", "high_var_probe"],
#                   predataset=["no", "yes"],num_features=bie.ANYATOM),platform=bie.ANYATOM,
#                   filter=bie.ANYATOM,group_fc=bie.ANYATOM)
                   
                   
all_modules = [
    #bie.Module(
    #    'make_normalize_report',
    #    [
    #     SignalFile2_rule.IntensityPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     SignalFile2_rule.ControlPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     SignalFile2_rule.SignalFile2(annotate='yes',preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     PcaAnalysis_rule.PcaAnalysis(annotate='no',preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     SignalFile2_rule.ActbPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     PcaAnalysis_rule.PcaAnalysis(annotate='no',preprocess=["unknown", "agilent", "mas5", "rma", "loess"],
    #            quantile_norm='no',combat_norm='no',shiftscale_norm='no',bfrm_norm='no',dwd_norm='no',gene_center='no',
    #            gene_normalize='no',),
    #    ],
    #     ReportFile(report_type='normalize',preprocess=["unknown", "agilent", "mas5", "rma", "loess"])),
    #bie.Module(
    #    'make_normalize_report',
    #    [
    #     #SignalFile2_rule.ControlPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     SignalFile2_rule.SignalFile2(annotate='yes',preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     PcaAnalysis_rule.PcaAnalysis(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
         #SignalFile2_rule.ActbPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #     PcaAnalysis_rule.PcaAnalysis(preprocess=["unknown", "agilent", "mas5", "rma", "loess"],
    #           quantile_norm='no',combat_norm='no',shiftscale_norm='no',bfrm_norm='no',dwd_norm='no',gene_center='no',
    #            gene_normalize='no',unique_genes="no",
    #            platform='no', group_fc='no'),
    #     SignalFile2_rule.IntensityPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
    #    
    #    ],
     #    ReportFile(report_type='normalize',preprocess=["unknown", "agilent", "mas5", "rma", "loess"])),
   bie.Module(
        'make_normalize_report',
        [
        SignalFile2_rule.SignalFile2(annotate='yes',preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
        
         
        PcaAnalysis_rule.PcaAnalysis(preprocess=["unknown", "agilent", "mas5", "rma", "loess"],
               quantile_norm='no',combat_norm='no',shiftscale_norm='no',bfrm_norm='no',dwd_norm='no',gene_center='no',
                gene_normalize='no',unique_genes="no",
                platform='no', group_fc='no'),
        PcaAnalysis_rule.PcaAnalysis(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
        SignalFile2_rule.IntensityPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
        ],
         ReportFile(report_type='normalize',preprocess=["unknown", "agilent", "mas5", "rma", "loess"])),
                                      
    bie.Module(
        'make_normalize_report_illumina',
       [
        SignalFile2_rule.SignalFile2(preprocess='illumina',annotate='yes'),
         PcaAnalysis_rule.PcaPlot(preprocess='illumina',quantile_norm='no',
              combat_norm='no',shiftscale_norm='no',bfrm_norm='no',dwd_norm='no'),
         PcaAnalysis_rule.PcaPlot,
         SignalFile2_rule.IntensityPlot,
         SignalFile2_rule.ActbPlot,
         SignalFile2_rule.BiotinPlot,
         SignalFile2_rule.HousekeepingPlot,
         SignalFile2_rule.Hyb_barPlot,
         SignalFile_rule.ControlFile(preprocess='illumina',format="gct", logged="no")
         ],
         ReportFile(report_type='normalize',preprocess='illumina')),
    
    bie.Module(
        'make_cluster_report',
        [ClusterFile_rule.ClusterFile,
         ClusterFile_rule.Heatmap],
        ReportFile(report_type='cluster',preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"])),
    bie.Module(
        'make_classify_report',
        [SignalFile2_rule.SignalFile2(
            contents='class0,class1,test',logged='yes',format='gct'),
         ClassifyFile_rule.ClassifyFile(classify_alg='svm',actual_label='yes'),
         ClassifyFile_rule.ClassifyFile(classify_alg='svm',loocv='yes',actual_label='no'),
         ClassifyFile_rule.PredictionPlot(classify_alg='svm',actual_label='yes',loocv='no'),      
         ClassifyFile_rule.PredictionPlot(classify_alg='svm',loocv='yes',actual_label='no'),
         ClassifyFile_rule.PredictionPCAPlot(classify_alg='svm',loocv='no',actual_label='yes'),
         ClassifyFile_rule.ClassifyFile(classify_alg='weighted_voting',loocv='no',actual_label='yes'),   
         ClassifyFile_rule.ClassifyFile(classify_alg='weighted_voting',loocv='yes',actual_label='no'),
         ClassifyFile_rule.PredictionPlot(classify_alg='weighted_voting',loocv='no',actual_label='yes'),
         ClassifyFile_rule.PredictionPlot(classify_alg='weighted_voting',loocv='yes',actual_label='no'),
         ClassifyFile_rule.PredictionPCAPlot(classify_alg='weighted_voting',loocv='no',actual_label='yes'),
         ],
        ReportFile(report_type='classify')),
    bie.Module(
        'make_heatmap_report',
        ClusterFile_rule.Heatmap(cluster_alg='no_cluster_alg'),
        ReportFile(report_type='heatmap',preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"])),
    bie.Module(
        'make_geneset_report',
        [GenesetAnalysis_rule.GenesetAnalysis,
         GenesetAnalysis_rule.GenesetPlot],
        ReportFile(report_type='geneset',preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"])),
    bie.Module(
        'make_diffgenes_report',
        [DiffExprFile_rule.DiffExprFile(diff_expr='t_test'),
         DiffExprFile_rule.DiffExprFile(diff_expr='sam'),
         ClusterFile_rule.Heatmap(cluster_alg='no_cluster_alg',hm_width='50',
                                  hm_height='1'),
         GatherFile_rule.GatherFile,
         GseaFile_rule.GseaFile
         ],
         ReportFile(report_type='diffgenes',preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"])),
    bie.Module(
        'make_batch_effect_report',
        [SignalFile_rule.SignalFile(format='tdf',logged='yes',missing_values='no',
                                    quantile_norm='yes',dwd_norm='no',bfrm_norm='no',combat_norm='no',shiftscale_norm='no'),
         SignalFile_rule.SignalFile(format='tdf',logged='yes',missing_values='no',
                                    quantile_norm='yes',dwd_norm='yes',bfrm_norm='no',combat_norm='no',shiftscale_norm='no'),
         SignalFile_rule.SignalFile(format='tdf',logged='yes',missing_values='no',
                                    quantile_norm='yes',shiftscale_norm='yes',dwd_norm='no',combat_norm='no',bfrm_norm='no'),
         SignalFile_rule.SignalFile(format='tdf',logged='yes',missing_values='no',
                                    quantile_norm='yes',bfrm_norm='yes',dwd_norm='no',combat_norm='no',shiftscale_norm='no'),
         SignalFile_rule.SignalFile(format='tdf',logged='yes',missing_values='no',
                                    quantile_norm='yes',combat_norm='yes',dwd_norm='no',bfrm_norm='no',shiftscale_norm='no')],
        ReportFile(report_type='batch_effect_remove',preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"])),
#  bie.Module(
  #      'combine all',
 #      [ReportFile(report_type='classify'),
 #       ReportFile(report_type='normalize'),
  #       ReportFile(report_type='geneset'),
 #       ReportFile(report_type='diffgenes'),
 #       ReportFile(report_type='batch_effect_remove'),
  #       ReportFile(report_type='heatmap'),
   #      ReportFile(report_type='cluster')],
 #     ReportFile(report_type='all',preprocess=["unknown", 'illumina', "agilent", "mas5", "rma", "loess"])),
         
        ]
  
        
