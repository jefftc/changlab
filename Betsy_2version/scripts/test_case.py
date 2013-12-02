from Betsy import rulebase
from Betsy import bie

def run_case1():
    # case1 (to generate the classify report network, which takes 757s)

    in_data = [
        rulebase.ClassLabelFile(
            contents='class0,class1', cls_format='cls',
            filename='/home/xchen/chencode/betsy_test/all_aml_train.cls'),
        rulebase.ClassLabelFile(
            contents='test', cls_format='cls',
            filename='/home/xchen/chencode/betsy_test/all_aml_test.cls'),
        rulebase.SignalFile(
            contents='class0,class1', format='tdf',
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
        rulebase.SignalFile(
            contents='test', format='tdf',
            filename='/home/xchen/chencode/betsy_test/all_aml_test.res')
        ]
    
    goal_datatype = rulebase.ReportFile
    goal_attributes = dict(report_type='classify')

    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)
    bie.print_network(network)
    bie.plot_network_gv("out.png", network)


def run_case2():
    #case2 (to generate the GenesetPlot network, which takes 1006s)
    
    in_data=[
        rulebase.GenesetFile(
            filename='/home/xchen/chencode/betsy_test/genesets.gmt'),
        rulebase.SignalFile(
            format='res',
            filename='/home/xchen/chencode/betsy_test/se2fplate6_48.illu.gz')
        ]
    
    goal_datatype = rulebase.GenesetPlot
    goal_attributes = dict(
        logged='yes',quantile_norm='yes',gene_center='mean',
        gene_normalize='variance',unique_genes='high_var',format='tdf',
        missing_values='no',annotate='yes',geneset='E2F1n_affy_150_UP')


    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)
    bie.print_network(network)
    bie.plot_network_gv("out.png", network)


def run_case3():
    #case 3 (to generate the SignalFile with quantile and combat method)
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
    network = bie.select_start_node(network, in_data)
    #order by quantile and then combat
    network = bie.remove_data_node(
        network, rulebase.SignalFile(quantile_norm="no", combat_norm="yes"))
    network = bie.optimize_network(network)

    bie.print_network(network)
    bie.plot_network_gv("out.png", network)


def run_case4():
    #case 4 (to generate the batch effect remove report)

    in_data=[
        rulebase.SignalFile(
            format="jeffs",
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
        rulebase.ClassLabelFile(
            filename='/home/xchen/chencode/betsy_test/all_aml_train.cls',
            cls_format='cls')]
    goal_datatype=rulebase.ReportFile
    goal_attributes=dict(report_type='batch_effect_remove')

    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)

    network = bie.remove_data_node(
        network, rulebase.SignalFile(quantile_norm="no", bfrm_norm="yes"))
    network = bie.remove_data_node(
        network, rulebase.SignalFile(quantile_norm="no", dwd_norm="yes"))
    network = bie.remove_data_node(
        network,
        rulebase.SignalFile(quantile_norm="no", shiftscale_norm="yes"))
    network = bie.remove_data_node(
        network, rulebase.SignalFile(quantile_norm="no", combat_norm="yes"))
    
    network = bie.optimize_network(network)

    bie.print_network(network)
    bie.plot_network_gv("out.png", network)


def run_case5():
    #case5 (to generate normalize report with quantile)

    # Why are there two analyze_samples_pca modules?
    in_data=[
        rulebase.SignalFile(
            format="jeffs",
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res')
        ]
    goal_datatype=rulebase.ReportFile
    goal_attributes=dict(
        report_type='normalize', format='tdf', logged='yes',
        missing_values='no', quantile_norm='yes',gene_center='mean',gene_normalize='variance',
        #missing_algorithm='zero_fill',
        #filter='20',num_features='500'
        )
    
    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    bie.print_network(network)
    bie.plot_network_gv("out1.png", network)

def run_case6():
    #case6 (to generate normalize report with illumina)
    in_data=[
        rulebase.ExpressionFiles(
            filename='/home/xchen/chencode/betsy_test/6991010018')
        ]
    goal_datatype=rulebase.ReportFile
    goal_attributes=dict(
        report_type='normalize',format='tdf',logged='yes',
        missing_values='no',preprocess='illumina')

    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)

    bie.print_network(network)
    bie.plot_network_gv("out.png", network)

def run_case7():
    # For optimization.
    in_data = [
        rulebase.SignalFile(
            contents='class0', format='tdf',
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
        ]
    
    goal_datatype = rulebase.SignalFile2
    goal_attributes = dict(contents='class0',quantile_norm='yes',
                           bfrm_norm='yes',gene_center='mean',
                           gene_normalize='variance')

    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    bie.print_network(network)
    bie.plot_network_gv("out.png", network)

def run_case7():
    # To optimize the optimize_network function.
    in_data = [
        rulebase.SignalFile(
            contents='class0', format='tdf',
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res'),
        ]

    goal_datatype = rulebase.SignalFile2
    goal_attributes = dict(
        contents='class0', quantile_norm='yes',
        bfrm_norm='yes', gene_center='mean',
        gene_normalize='variance')

    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    bie.print_network(network)
    bie.plot_network_gv("out.png", network)

def run_case8():
    #preprocess value cannot pass to the Heatmap
    in_data = [
        rulebase.ExpressionFiles(
            filename='/home/xchen/chencode/betsy_test/6991010018'),
        ]
    #goal_datatype = rulebase.Heatmap
    goal_datatype = rulebase.ReportFile
    goal_attributes = dict(
        report_type='cluster',format='tdf',logged='yes',
        missing_values='no',preprocess='illumina',quantile_norm='yes',bfrm_norm='yes')
    
    # No ReportFile with report_type='heatmap' and
    # preprocess='illumina'.  To avoid ambiguity, should make sure no
    # DataTypes have the same names for attributes.  e.g. don't reuse
    # "preprocess" for multiple DataTypes.


    #goal_attributes = dict(
    #    report_type='normalize', format='tdf', logged='yes',
    #    missing_values='no', preprocess='illumina')
    #    #missing_values='no')

  #  goal_attributes = dict(
   #     report_type='heatmap', format='tdf', logged='yes',
  #      missing_values='no', preprocess='unknown', quantile_norm="yes")
  #      #missing_values='no')


    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    bie.print_network(network)
    bie.plot_network_gv("out.png", network)

def run_case9():
    #case9 (to test the different order of PcaAnalysis in report input list)
    in_data=[
        rulebase.SignalFile(
            format="jeffs",
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res')
        ]
    goal_datatype=rulebase.ReportFile
    goal_attributes=dict(
        report_type='normalize', format='tdf', logged='yes',
        missing_values='no', quantile_norm='yes',gene_center='mean',
        gene_normalize='variance',
        )
    new_module = bie.Module(
        'make_normalize_report',
        [
        rulebase.SignalFile2(annotate='yes',preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
        rulebase.PcaAnalysis(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
        rulebase.IntensityPlot(preprocess=["unknown", "agilent", "mas5", "rma", "loess"]),
        rulebase.PcaAnalysis(preprocess=["unknown", "agilent", "mas5", "rma", "loess"],
               quantile_norm='no',combat_norm='no',shiftscale_norm='no',bfrm_norm='no',dwd_norm='no',gene_center='no',
                gene_normalize='no',unique_genes="no",
                platform='no', group_fc='no'),
        ],
        rulebase.ReportFile(report_type='normalize',preprocess=["unknown", "agilent", "mas5", "rma", "loess"]))
    all_modules = rulebase.all_modules
    all_modules.append(new_module)    
    network = bie.backchain(
        all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    bie.print_network(network)
    bie.plot_network_gv("out.png", network)
    
    
def main():
    import cProfile 
    #run_case1()
    #run_case2()
    #run_case3()
    #run_case4()
    #run_case5()
    #run_case6()
    #run_case7()
    #run_case8()
    run_case9()
    #cProfile.run("run_case7()")

if __name__ == '__main__':
    main()
