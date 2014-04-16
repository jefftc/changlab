from Betsy import rulebase
from Betsy import bie3

def run_case1():
    in_data = rulebase.GEOSeries
    out_data = rulebase.SignalFile.output(preprocess="rma",
        format="tdf", logged="yes",gene_center='mean',#annotate='yes',
        missing_values="no",quantile_norm='yes',#contents="class0,class1"
                                          )
    
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    

    print "INPUT:"
    print in_data
    print
    
    print "OUTPUT:"
    print out_data
    print
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)


def run_case2():
    in_data = rulebase.GEOSeries

    # Will generate network back to illumina preprocessing if
    # SignalFile2 is given.  Problem is that SignalFile cannot be
    # shiftscale normalized.
    out_data = rulebase.SignalFile.output(
        preprocess="illumina",
        format="tdf", logged="yes",
        missing_values="no", shiftscale_norm='yes'
        )
    
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)

    print "INPUT:"
    print in_data
    print
    
    print "OUTPUT:"
    print out_data
    print
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case3():
    in_data = rulebase.GEOSeries
    out_data = rulebase.SignalFile.output(preprocess="illumina",
        format="tdf",  logged="yes",
        missing_values="no")
    
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)

    print "INPUT:"
    print in_data
    print
    
    print "OUTPUT:"
    print out_data
    print
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case4():
    in_data = rulebase.GEOSeries

    # The SignalFile2 should be created by the reorder_genes module.
    # However, it can not create it because gene_center="no", by
    # default.  reorder_genes produces a SignalFile2 with
    # gene_center="unknown", which conflicts.  SignalFile2 has no way
    # to check_gene_center.
    #
    # Work around is to make gene_center="unknown" and
    # gene_normalize="unknown".  Better solution is to rethink how the
    # SignalFiles work.
    out_data = rulebase.PrettySignalFile.output(
        preprocess="illumina",
        format="tdf", logged="yes",
        missing_values="no", gene_order='t_test_p',
        #gene_center="unknown", gene_normalize="unknown",
        )
    
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)

    print "INPUT:"
    print in_data
    print
    
    print "OUTPUT:"
    print out_data
    print
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case5():
    """ for each module,the attributes not mentioned will
    be set to its default input value."""
    in_data = rulebase.GEOSeries
    out_data = rulebase.SignalFile.output(
        preprocess="agilent", format="tdf",  quantile_norm='yes')

    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    
    print "INPUT:"
    print in_data
    print
    
    print "OUTPUT:"
    print out_data
    print

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
                                                                    

def run_case6():
    network = bie3.backchain(
        rulebase.all_modules, rulebase.ActbPlot,
        bie3.Attribute(rulebase.SignalFile, "preprocess", "rma"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case7():
    network = bie3.backchain(
        rulebase.all_modules, rulebase.SignalFile,
        bie3.Attribute(rulebase.SignalFile,"contents","class0,class1"),
         bie3.Attribute(rulebase.SignalFile,"preprocess","rma"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    
def run_case8():
    #test ClusterFile

    # Heatmap requires SignalFile to be logged.  Explicitly
    # specificying logged=yes changes the network, even though they
    # should in principle be the same.
    in_data = rulebase.GEOSeries
    network = bie3.backchain(
        rulebase.all_modules, rulebase.Heatmap,
        ###specify this attribute or not make the network different
        bie3.Attribute(rulebase.SignalFile, "logged", "yes"),
        )
    network = bie3.optimize_network(network)

    print "INPUT:"
    print in_data
    print
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    

def run_case9():
    # command1 (command 1 and command 2 suppose to have the same
    # result, but they are not)

    # command 1
    out_data = rulebase.SignalFile.output(
        preprocess="rma",quantile_norm='yes',
        gene_center='mean',gene_normalize='variance')
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    bie3.print_network(network, open("out1.log", 'w'))
    bie3.plot_network_gv("out1.png", network)

    # command 2 
    network = bie3.backchain(  
        rulebase.all_modules, rulebase.SignalFile,
        bie3.Attribute(rulebase.SignalFile,"preprocess","rma"),
        bie3.Attribute(rulebase.SignalFile,"quantile_norm","yes"),
        bie3.Attribute(rulebase.SignalFile,'gene_center',"mean"),
        bie3.Attribute(rulebase.SignalFile,'gene_normalize',"variance"))
    network = bie3.optimize_network(network)
    bie3.print_network(network, open("out2.log", 'w'))
    bie3.plot_network_gv("out2.png", network)


def run_case10():
    # the SignalFile has several preprocess not only 'mas5'
    out_data = rulebase.SignalFile.output(
        preprocess='mas5', contents="class0,class1")

    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)


def run_case11():
    # New version of bie3 (2/20/14) runs too closly and generates
    # "network too large" error.  Older version finishes quickly.

    if 0:
        # No problems.
        out_data = rulebase.SignalFile.output()
        network = bie3.backchain(rulebase.all_modules, out_data)
        network = bie3.optimize_network(network)
    else:
        # network too large.
        out_data = rulebase.PrettySignalFile.output()
        network = bie3.backchain(rulebase.all_modules, out_data)
        network = bie3.optimize_network(network)
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)


def run_case12():
    # the branches to merge module has only one GeoSeries, it supposed
    # to have two, one is contents=class0, one is contents=class1
    
    #out_data = rulebase.SignalFile.output(
    #    contents='class0,class1',preprocess='mas5')
    out_data = rulebase.SignalFile.output(
        bfrm_norm='no', combat_norm='no', contents='class1',
        dwd_norm='no', filter='no', format='tdf', gene_center='no',
        gene_normalize='no', logged='yes', missing_algorithm='zero_fill',
        missing_values='no', predataset='no', preprocess='mas5',
        processing_step='merge', quantile_norm='no', shiftscale_norm='no')
    
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    
def run_case13():
    '''test cluster report,
       ClusterFile cluster_alg should be pca, but it shows
       four different cluster algorithms'''
    network = bie3.backchain(  
        rulebase.all_modules, rulebase.ReportFile,
         bie3.Attribute(rulebase.PrettySignalFile,"preprocess","mas5"),
         bie3.Attribute(rulebase.ReportFile,"report_type","cluster"),
         bie3.Attribute(rulebase.PrettySignalFile,"quantile_norm","yes"),
         bie3.Attribute(rulebase.ClusterFile,"cluster_alg","pca"),
         bie3.Attribute(rulebase.Heatmap,"cluster_alg","pca"),
        )
    network = bie3.optimize_network(network)
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case14():
    '''test normalize report,
       requires PSF preprocess=unknown and contents=test,
       but preprocess can be any of ['rma', 'mas5', 'agilent',
       'loess', 'unknown'] and contents can any of
       ['train0', 'train1', 'test', 'class0,class1,test',
       'class0', 'class1', 'class0,class1','unspecified']'''
    network = bie3.backchain(  
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.ReportFile,"report_type","normalize_file"),
        bie3.Attribute(rulebase.PrettySignalFile,"preprocess","unknown"),
        bie3.Attribute(rulebase.PrettySignalFile,"contents","test"),
        
        )
    network = bie3.optimize_network(network)
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case15():
    """want PSF preprocess=illumina, but PSF that goes
       into rank_genes_by_class_neighbors has preprocess
       unknown"""
    out_data = rulebase.PrettySignalFile.output(
        gene_order='class_neighbors',preprocess='illumina')                                 
    network = bie3.backchain(  
        rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    

def run_case16():
    """the difference between node 59 and node 191 is the sf_processing_step,
       if we have an input SignalFile as node 66, the pipeline will go to node 59
       but no way to go module 6."""
    out_data = rulebase.SignalFile.output(
        dwd_norm='yes')                                 
    network = bie3.backchain(  
        rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case17():
    '''test the 'network too large error',
       I have changed the MAX_NETWORK_SIZE to 10024, the
       out network is about 768 nodes and does not pop
       'network too large' error'''
    network = bie3.backchain(  
        rulebase.all_modules, rulebase.ClassifyFile,
        bie3.Attribute(rulebase.ClassifyFile,"classify_alg","weighted_voting"),
        bie3.Attribute(rulebase.PrettySignalFile,"quantile_norm","yes")
        )
    network = bie3.optimize_network(network)
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    
def main(): 
    #run_case1()
    #run_case2()
    #run_case3()
    #run_case4()
    #run_case5()
    #run_case6()
    #run_case7()
    #run_case8()
    #run_case9()
    #run_case10()
    #run_case11()
    #run_case12()
    #run_case13()
    #run_case14()
    #run_case15()
    #run_case16()
    run_case17()
if __name__ == '__main__':
    main()
