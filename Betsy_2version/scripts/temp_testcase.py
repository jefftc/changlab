#case1
##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##         filename='/home/xchen/chencode/betsy_test/all_aml_train.res')]
##goal_datatype = rulebase.SignalFile
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="rma",missing_values='no')



#case2
##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res')]
##goal_datatype = rulebase.SignalFile
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="rma",
##    missing_values='no',bfrm_norm='yes',quantile_norm='yes'
##    )
        
#case3

##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls')]
##goal_datatype = rulebase.SignalFile2
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="rma",missing_algorithm='zero_fill',
##    missing_values='no',)#group_fc='2',gene_center='mean',gene_normalize='variance')
    #platform='HG_U133A',annotate='yes',duplicate_probe='closest_probe')
#,unique_genes='high_var')
#,),
   # ,,)

#case4
##in_data = [SignalFile_rule.SignalFile(contents='class0',preprocess="rma", format="jeffs",
##                      filename='/home/xchen/chencode/betsy_test/er.l2.mas5.train0'),
##           SignalFile_rule.SignalFile(contents='class1',preprocess="rma", format="jeffs",
##                      filename='/home/xchen/chencode/betsy_test/er.l2.mas5.train1')]
##goal_datatype = SignalFile_rule.SignalFile
##goal_attributes = dict(contents='class0,class1',
##    format='tdf', logged='yes', preprocess="rma",
##    missing_values='no')

#case5
##in_data=[rulebase.GEOSeries,
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls')]
##goal_datatype = rulebase.SignalFile
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="rma",missing_values='no',GSEID='GSE8286',
##        combat_norm='yes')

#case6
##in_data=[rulebase.ExpressionFiles(filename='/home/xchen/chencode/betsy_test/6991010018')]
###in_data=[rulebase.SignalFile(format='tdf',logged='yes',missing_values='no')]
##goal_datatype = rulebase.SignalFile2
####goal_datatype=rulebase.ControlFile
####goal_datatype=Heatmap
####goal_attributes = dict(format='tdf',logged='yes',
####                       cluster_alg='som',missing_values='no',preprocess='illumina')
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="illumina",missing_values='no',platform='HG_U133A',
##    duplicate_probe='closest_probe',num_features='500'
##)

#case7
##in_data=[rulebase.ExpressionFiles(filename='/home/xchen/chencode/betsy_test/agilent_expression')]
##goal_datatype = rulebase.SignalFile
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="agilent",missing_values='no',
##)

#case 8

##in_data=[rulebase.ExpressionFiles(filename='/home/xchen/chencode/betsy_test/GSE4189')]
##goal_datatype = rulebase.SignalFile
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="loess",missing_values='no',missing_algorithm='zero_fill',
##)
#case 9

##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res')]
##goal_datatype = rulebase.ClusterFile
##goal_datatype = rulebase.Heatmap
##goal_attributes = dict(cluster_alg='hierarchical',
##    format='tdf', logged='yes', preprocess="rma",
##    missing_values='no',predataset='yes')

#case 10
##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/Chang_AR00410.pcl'),
##         rulebase.RenameFile(filename='/home/xchen/chencode/betsy_test/rename_list_file.txt')]
##goal_datatype = rulebase.SignalFile2
##goal_attributes = dict(
##    format='tdf', logged='yes', preprocess="rma",
##    missing_values='no',rename_sample='yes',gene_center='mean')

##case11

##in_data=[rulebase.ClassLabelFile(contents='class0,class1',filename='/home/xchen/chencode/betsy_test/all_aml_train.cls'),
##         rulebase.ClassLabelFile(contents='test',filename='/home/xchen/chencode/betsy_test/all_aml_test.cls'),
##         rulebase.SignalFile2(contents='class0,class1',format='tdf',filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.SignalFile2(contents='test',format='tdf',filename='/home/xchen/chencode/betsy_test/all_aml_test.res')
##         ]
##goal_datatype=rulebase.ClassifyFile
###goal_datatype=rulebase.SignalFile2
###goal_attributes=dict(contents='class0,class1',align='yes_for_train')
##goal_attributes=dict(classify_alg='weighted_voting',actual_label='yes')

#case12

##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls',cls_format='cls')]
##goal_datatype = rulebase.DiffExprFile
##goal_attributes = dict(
##     logged='yes', preprocess="rma",
##    missing_values='no',diff_expr='t_test')
###case13
##in_data=[rulebase.GenesetFile(filename='/home/xchen/chencode/betsy_test/genesets.gmt'),
##         rulebase.SignalFile(format='res',filename='/home/xchen/chencode/betsy_test/se2fplate6_48.illu.gz')]
##goal_datatype=rulebase.GenesetPlot
##goal_attributes=dict(logged='yes',quantile_norm='yes',gene_center='mean',
##                     gene_normalize='variance',unique_genes='high_var',format='tdf',
##                     missing_values='no',annotate='yes',geneset='E2F1n_affy_150_UP')

#case14

##in_data=[rulebase.SignalFile(preprocess="rma", format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         SignalFile_rule.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls')]
##goal_datatype = rulebase.GseaFile
##goal_attributes = dict(
##     logged='yes', preprocess="rma",
##    missing_values='no')
#case15


##in_data=[rulebase.CELFiles(filename='/home/xchen/chencode/betsy_test/GSE8286_folder')]
###in_data=[SignalFile2(format='tdf',preprocess='rma',logged='yes',
###                     platform='HG_U133A',duplicate_probe='high_var_probe'),
###         SignalFile2(format='tdf',preprocess='mas5',logged='yes',
###                     platform='HG_U133A',duplicate_probe='high_var_probe')]
##goal_datatype = rulebase.SignatureScore 
##goal_attributes = dict()

#case16

##in_data=[rulebase.SignalFile(format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls')]
##in_data=[rulebase.GeneListFile(filename='/home/xchen/chencode/betsy_test/gene_list.txt')]
##goal_datatype=DavidFile_rule.DavidFile
##goal_attributes=dict()

#case17

##in_data=[rulebase.SignalFile(format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls')]
##goal_datatype=rulebase.PcaPlot
##goal_attributes=dict(process='before')
   
#case18

##in_data=[rulebase.SignalFile(format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(cls_format='cls',filename='/home/xchen/chencode/betsy_test/all_aml_train.cls')
##         ]
##in_data=[rulebase.ExpressionFiles(filename='/home/xchen/chencode/betsy_test/6991010018')]
##goal_datatype=rulebase.ReportFile
##goal_attributes=dict(report_type='normalize',format='tdf',logged='yes',
##                     missing_values='no',preprocess='illumina')

#case19
##in_data=[rulebase.ClassLabelFile(contents='class0,class1',cls_format='cls',filename='/home/xchen/chencode/betsy_test/all_aml_train.cls'),
##         rulebase.ClassLabelFile(contents='test',cls_format='cls',filename='/home/xchen/chencode/betsy_test/all_aml_test.cls'),
##         rulebase.SignalFile(contents='class0,class1',format='tdf',filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.SignalFile(contents='test',format='tdf',filename='/home/xchen/chencode/betsy_test/all_aml_test.res')
##         ]
##goal_datatype=rulebase.ReportFile
###goal_datatype=rulebase.PredictionPCAPlot
###goal_datatype=rulebase.PcaAnalysis
###goal_attributes = dict(process='after',contents='test')
###goal_attributes = dict(contents='test',logged='yes',format='gct')
###goal_attributes = dict(classify_alg='svm',actual_label='yes',loocv='no')
##goal_attributes=dict(report_type='classify')

#case20
##in_data=[rulebase.ClassLabelFile(cls_format='label',filename='/home/xchen/chencode/betsy_test/chang_AR.label.txt'),
##         rulebase.ExpressionFiles(filename='/home/xchen/chencode/betsy_test/Chang_AR00410')]
##goal_datatype=rulebase.SignalFile
##goal_attributes=dict(preprocess='illumina',format='tdf',dwd_norm='yes',logged='yes',missing_values='no')

#case21
##in_data=[rulebase.SignalFile(format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res')]
##goal_datatype=rulebase.ReportFile
##goal_attributes=dict(report_type='heatmap')
##   
#case22
##in_data=[rulebase.GenesetFile(filename='/home/xchen/chencode/betsy_test/genesets.gmt'),
##         rulebase.SignalFile(format='res',filename='/home/xchen/chencode/betsy_test/se2fplate6_48.illu.gz')]
##goal_datatype=rulebase.ReportFile
##goal_attributes=dict(logged='yes',quantile_norm='yes',gene_center='mean',
##                     gene_normalize='variance',unique_genes='high_var',format='tdf',
##                     missing_values='no',annotate='yes',geneset='E2F1n_affy_150_UP',report_type='geneset')

#23
##in_data=[rulebase.SignalFile(format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls',cls_format='cls')]
##goal_datatype = rulebase.ReportFile
##goal_attributes=dict(
##    report_type='diffgenes',gene_order='t_test_p')

#goal_attributes=dict(logged='yes', preprocess="rma",
#    missing_values='no',diff_expr='t_test')
##in_data=[rulebase.ExpressionFiles(filename='/home/xchen/chencode/betsy_test/6991010018'),
##         rulebase.ClassLabelFile(cls_format='cls',filename='/home/xchen/chencode/betsy_test/6991010018.cls')]
##goal_datatype = rulebase.ReportFile
##goal_attributes = dict(
##      preprocess="illumina",#missing_values='no',logged='yes',
##      report_type='normalize')#,gene_order='t_test_p')#,report_type='normalize')
##goal_datatype=rulebase.PcaPlot
##goal_attributes=dict(logged='yes', preprocess="rma",
##    missing_values='no',process='before',format='tdf')
##
##in_data=[rulebase.SignalFile(format="jeffs",
##                    filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.ClassLabelFile(filename='/home/xchen/chencode/betsy_test/all_aml_train.cls',cls_format='cls')]
##goal_datatype=rulebase.ReportFile
##goal_attributes=dict(report_type='batch_effect_remove')
##network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
##network = bie.optimize_network(network)
##network = bie.prune_network_by_start(network, in_data)
####network = bie.prune_network_by_internal(
####        network, rulebase.SignalFile(quantile_norm="yes", bfrm_norm="no"))
####network = bie.prune_network_by_internal(
####       network, rulebase.SignalFile(quantile_norm="yes", dwd_norm="no"))
####network = bie.prune_network_by_internal(
####       network, rulebase.SignalFile(quantile_norm="yes", shiftscale_norm="no"))
####network = bie.prune_network_by_internal(
####       network, rulebase.SignalFile(quantile_norm="yes", combat_norm="no"))
##bie._print_network(network)
##bie._plot_network_gv("out.png", network)
######import pickle
####f=file('network_classify_weightedvoting','rb')
########pickle.dump(network,f)
##network = pickle.load(f)
##f.close()
####print 'done'
#run_pipeline(network,in_data)
##def test_bie():
##    in_data=[rulebase.ClassLabelFile(contents='class0,class1',cls_format='cls',filename='/home/xchen/chencode/betsy_test/all_aml_train.cls'),
##         rulebase.ClassLabelFile(contents='test',cls_format='cls',filename='/home/xchen/chencode/betsy_test/all_aml_test.cls'),
##         rulebase.SignalFile(contents='class0,class1',format='tdf',filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
##         rulebase.SignalFile(contents='test',format='tdf',filename='/home/xchen/chencode/betsy_test/all_aml_test.res')
##         ]
##    goal_datatype=rulebase.ReportFile
##    goal_attributes=dict(report_type='classify')
##    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
##    network = bie.optimize_network(network)
##    network = bie.prune_network_by_start(network, in_data)
##    bie._print_network(network)
##    bie._plot_network_gv("out.png", network)
##    import pickle
##    f=file('network_classify_report','wb')
##    pickle.dump(network,f)
##    f.close()
##import cProfile; cProfile.run("test_bie()")

##def test_bie():
##    in_data=[rulebase.GenesetFile(filename='/home/xchen/chencode/betsy_test/genesets.gmt'),
##         rulebase.SignalFile(format='res',filename='/home/xchen/chencode/betsy_test/se2fplate6_48.illu.gz')]
##    #goal_datatype=rulebase.GenesetPlot
####    goal_attributes=dict(logged='yes',quantile_norm='yes',gene_center='mean',
####                     gene_normalize='variance',unique_genes='high_var',format='tdf',
####                     missing_values='no',annotate='yes',geneset='E2F1n_affy_150_UP')
##    goal_datatype=rulebase.ReportFile
##    goal_attributes=dict(logged='yes',quantile_norm='yes',gene_center='mean',
##                     gene_normalize='variance',unique_genes='high_var',format='tdf',
##                     missing_values='no',annotate='yes',geneset='E2F1n_affy_150_UP',report_type='geneset')
##    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
##    network = bie.optimize_network(network)
##    network = bie.prune_network_by_start(network, in_data)
##    bie._print_network(network)
##    bie._plot_network_gv("out.png", network)
##    import pickle
##    #f=file('network_genesetplot','wb')
##    f=file('network_genesetnetwork','wb')
##    pickle.dump(network,f)
##    f.close()
from Betsy import rulebase
from Betsy import bie
def test_bie():
     in_data=[
        rulebase.SignalFile(
            preprocess="rma", format="jeffs",
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
        rulebase.ClassLabelFile(
            filename='/home/xchen/chencode/betsy_test/all_aml_train.cls',
            cls_format='cls')
        ]
     goal_datatype = rulebase.SignalFile
     goal_attributes = dict(
        format='tdf', logged='yes', preprocess="rma",
        missing_values='no',combat_norm='yes',quantile_norm='yes'
        )

     network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
     network = bie.optimize_network(network)
     network = bie.prune_network_by_start(network, in_data)
 
     #order by quantile and then combat
     network = bie.prune_network_by_internal(
        network, rulebase.SignalFile(quantile_norm="yes", combat_norm="no"))

     bie._print_network(network)
     bie._plot_network_gv("out.png", network)

import cProfile; cProfile.run("test_bie()")
#test_bie()
