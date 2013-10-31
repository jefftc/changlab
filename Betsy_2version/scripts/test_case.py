from Betsy import rulebase
from Betsy import bie

def run_case1():
    #case1 (to generate the classify report network,which take 757s)

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
    bie._print_network(network)
    bie._plot_network_gv("out.png", network)


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
    bie._print_network(network)
    bie._plot_network_gv("out.png", network)


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
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)

    #order by quantile and then combat
    network = bie.prune_network_by_internal(
        network, rulebase.SignalFile(quantile_norm="yes", combat_norm="no"))

    bie._print_network(network)
    bie._plot_network_gv("out.png", network)


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
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)

    network = bie.prune_network_by_internal(
        network, rulebase.SignalFile(quantile_norm="yes", bfrm_norm="no"))
    network = bie.prune_network_by_internal(
        network, rulebase.SignalFile(quantile_norm="yes", dwd_norm="no"))
    network = bie.prune_network_by_internal(
        network,
        rulebase.SignalFile(quantile_norm="yes", shiftscale_norm="no"))
    network = bie.prune_network_by_internal(
        network, rulebase.SignalFile(quantile_norm="yes", combat_norm="no"))

    bie._print_network(network)
    bie._plot_network_gv("out.png", network)


def run_case5():
    #case5 (to generate normalize report with quantile)
    in_data=[
        rulebase.SignalFile(
            format="jeffs",
            filename='/home/xchen/chencode/betsy_test/all_aml_train.res')
        ]
    goal_datatype=rulebase.ReportFile
    goal_attributes=dict(
        report_type='normalize', format='tdf', logged='yes',
        missing_values='no', quantile_norm='yes')
    
    network = bie.backchain(
        rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)
    bie._print_network(network)
    bie._plot_network_gv("out.png", network)

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
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)

    bie._print_network(network)
    bie._plot_network_gv("out.png", network)


def main():
    run_case1()
    #run_case2()
    #run_case3()
    #run_case4()
    #run_case5()
    #run_case6()


if __name__ == '__main__':
    main()
