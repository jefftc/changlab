# ReportFile
import bie
from bie_rules import SignalFile_rule,SignalFile2_rule,ClusterFile_rule,ClassifyFile_rule,PcaAnalysis_rule

ReportFile = bie.DataType(
    'ReportFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(report_type=['normalize','batch_effect_remove',
                            'classify','cluster','diffgenes',
                            'heatmap','geneset'],DEFAULT='normalize'),
    )

list_files = [ReportFile]

all_modules = [
    bie.Module(
        'make_normalize_report',
        [SignalFile2_rule.SignalFile2(annotate='yes'),
         PcaAnalysis_rule.PcaPlot(process='before'),
         PcaAnalysis_rule.PcaPlot(process='after')],
        ReportFile(report_type='normalize')),
    
    bie.Module(
        'make_cluster_report',
        [ClusterFile_rule.ClusterFile,
         ClusterFile_rule.Heatmap],
        ReportFile(report_type='cluster')),
    bie.Module(
        'make_classify_report',
        [SignalFile2_rule.SignalFile2(
            contents='class0,class1,test',logged='yes',format='gct'),
         ClassifyFile_rule.ClassifyFile(classify_alg='svm',actual_label='yes'),
         ClassifyFile_rule.ClassifyFile(classify_alg='svm',loocv='yes',actual_label='no'),
         ClassifyFile_rule.ClassifyFile(classify_alg='weighted_voting',loocv='no',actual_label='yes'),
         ClassifyFile_rule.ClassifyFile(classify_alg='weighted_voting',loocv='yes',actual_label='no'),         ],
        ReportFile(report_type='classify')),
         
        ]

         
