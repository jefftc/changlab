from Betsy import rulebase
from Betsy import bie3

def run_case1():
    in_data = rulebase.GEOSeries
    out_data = rulebase.SignalFile.output(preprocess="illumina",
        format="tdf",  logged="yes",
        missing_values="no", contents="class0,class1")
    
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
    out_data = rulebase.SignalFile.output(preprocess="illumina",
        format="tdf",  logged="yes",
        missing_values="no",shiftscale_norm='yes')
    
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
    out_data = rulebase.SignalFile2.output(preprocess="illumina",
        format="tdf",  logged="yes",
        missing_values="no",gene_order='t_test_p')
    
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

def main(): 
    #run_case1()
    #run_case2()
    #run_case3()
    run_case4()
if __name__ == '__main__':
    main()
