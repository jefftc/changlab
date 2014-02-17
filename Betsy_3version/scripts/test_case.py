from Betsy import rulebase
from Betsy import bie3

def run_case1():
    in_data = rulebase.GEOSeries
    out_data = rulebase.SignalFile.output(preprocess="rma",
        format="tdf", logged="yes",#gene_center='mean',annotate='yes',
        missing_values="no",quantile_norm='yes',contents="class0,class1")
    
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
    out_data = rulebase.SignalFile2.output(
        preprocess="illumina",
        format="tdf", logged="yes",
        missing_values="no", shiftscale_norm='yes')
    
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
    out_data = rulebase.SignalFile1.output(preprocess="illumina",
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
    out_data = rulebase.SignalFile2.output(
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
    out_data = rulebase.SignalFile2.output(
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
        bie3.Attribute(rulebase.SignalFile1,"contents","class0,class1"),
         bie3.Attribute(rulebase.SignalFile1,"preprocess","rma"),
         bie3.Attribute(rulebase.SignalFile, "preprocess", "rma"),
         bie3.Attribute(rulebase.SignalFile, "preprocess", "rma")
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    
def run_case8():
    #test ClusterFile
    in_data = rulebase.GEOSeries
    network = bie3.backchain(
        rulebase.all_modules, rulebase.Heatmap,
        ###specify this attribtue or not make the network different
        bie3.Attribute(rulebase.SignalFile, "logged", "yes"),
        )
    network = bie3.optimize_network(network)

    print "INPUT:"
    print in_data
    print
    
    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    
def main(): 
    run_case1()
    #run_case2()
    #run_case3()
    #run_case4()
    #run_case5()
    #run_case6()
    #run_case7()
    #run_case8()
    

if __name__ == '__main__':
    main()
