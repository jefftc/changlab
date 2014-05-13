from Betsy import rulebase
from Betsy import bie3

def run_case01():
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


def run_case02():
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

def run_case03():
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

def run_case04():
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

def run_case05():
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


def run_case06():
    network = bie3.backchain(
        rulebase.all_modules, rulebase.ActbPlot,
        bie3.Attribute(rulebase.SignalFile, "preprocess", "rma"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case07():
    network = bie3.backchain(
        rulebase.all_modules, rulebase.SignalFile,
        bie3.Attribute(rulebase.SignalFile,"contents","class0,class1"),
         bie3.Attribute(rulebase.SignalFile,"preprocess","rma"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case08():
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


def run_case09():
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
        bie3.Attribute(rulebase.PrettySignalFile, "preprocess", "mas5"),
        bie3.Attribute(rulebase.ReportFile, "report_type", "cluster"),
        bie3.Attribute(rulebase.PrettySignalFile, "quantile_norm", "yes"),
        bie3.Attribute(rulebase.ClusterFile, "cluster_alg", "pca"),
        bie3.Attribute(rulebase.Heatmap, "cluster_alg", "pca"),
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
        #bie3.Attribute(rulebase.PrettySignalFile,"quantile_norm","yes")
        
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case15():
    # want PSF preprocess=illumina, but PSF that goes into
    # rank_genes_by_class_neighbors has preprocess unknown.
    #
    # Problem: Why does no PSF with preprocess=illumina point to
    # rank_genes_by_class_neighbors?
    # Answer: rank_genes_by_class_neighbors takes PrettySignalFile.
    # In default PrettySignalFile, output preprocess=unknown.

    #out_data = rulebase.PrettySignalFile.output(
    #    gene_order='class_neighbors', preprocess='illumina')
    #network = bie3.backchain(rulebase.all_modules, out_data)
    #network = bie3.optimize_network(network)
    #bie3.write_network("test.network", network)
    #network = bie3.read_network("test.network")
    #bie3.complete_network(network)
    #network = bie3.optimize_network(network)

    out_data = rulebase.PrettySignalFile.output(
        gene_order='class_neighbors', preprocess='illumina')
    network = bie3.backchain(rulebase.all_modules, out_data)

    network = bie3.complete_network(network)
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)


def run_case16():
    """the difference between node 59 and node 191 is the
       sf_processing_step, if we have an input SignalFile as node 66,
       the pipeline will go to node 59 but no way to go module 6."""

    #  59.  SignalFile  gene_normalize="no"
    #                   sf_processing_step="processed"
    #  66.  SignalFile  gene_normalize="unknown"
    #                   sf_processing_step="normalize"
    # 191.  SignalFile  gene_normalize="no"
    #                   sf_processing_step="merge"

    # Problem: Input file with gene_normalize="unknown" cannot be used
    # to normalize_samples_with_dwd.
    # SignalFile [59] should be acceptable as input for
    # normalize_samples_with_dwd.
    #  66 -> check_gene_normalize -> 59 -> convert_label_to_cls ->
    #    21 -> normalize_samples_with_dwd [6]
    # 191 -> normalize_samples_with_dwd [6]
    #
    # normalize_samples_with_dwd requires sf_processing_step to be
    # "merge".
    #
    # Is this a problem?  Node 66 should not be an input.  Inputs
    # should have an earlier processing step (e.g. "postprocess").  In
    # the network, nodes higher up do go into
    # normalize_samples_with_dwd [6].

    # Processing steps:
    # postprocess -> impute -> merge -> normalize -> processed

    out_data = rulebase.SignalFile.output(dwd_norm='yes')
    network = bie3.backchain(rulebase.all_modules, out_data)
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



def run_case18():
    """Result that generates:
    A network with 232 nodes. Node 2 and Node 3 are following:
    2.  Data(ClusterFile, cluster_alg='pca', contents='unspecified',
    distance=['correlation', 'euclidean'])
    3.  Data(Heatmap, cluster_alg='pca', color='red_green',
    contents='unspecified', distance=['correlation', 'euclidean'],
    hm_height='yes', hm_width='yes')

    Result I expected:
    distance in Node 2 and Node 3 should be set to default because we
    did not specify it.

    That is: distance='correlation'.

    JC: The distance is specified in the make_cluster_report Module:
    Constraint("distance", CAN_BE_ANY_OF, ['correlation','euclidean'], 0),

    Defaults are used only if no other information is available.

    """

    network = bie3.backchain(
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.PrettySignalFile, "preprocess", "mas5"),
        bie3.Attribute(rulebase.ReportFile, "report_type", "cluster"),
        bie3.Attribute(rulebase.PrettySignalFile, "quantile_norm", "yes"),
        bie3.Attribute(rulebase.ClusterFile, "cluster_alg", "pca"),
        bie3.Attribute(rulebase.Heatmap, "cluster_alg", "pca"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)


def run_case19():
    """Result that generates:
    A network with 127 nodes. Node 2 and Node 3 are following:
    2.  Data(ClusterFile, cluster_alg=['som', 'pca', 'kmeans', 'hierarchica
        l'], contents='unspecified', distance=['correlation', 'euclidean'])
    3.  Data(Heatmap, cluster_alg=['som', 'pca', 'kmeans', 'hierarchical'],
        color='red_green', contents='unspecified', distance=['correlatio
        n', 'euclidean'], hm_height='yes', hm_width='yes')

    Result I expected:
    distance and cluster_alg in Node 2 and Node 3 should be set
    to default because we did not specify it.
    That is: distance='correlation', cluster_alg = 'kmeans'.

    JC: The distance and cluster_alg is specified in the
    make_cluster_report Module:
    Constraint("cluster_alg",CAN_BE_ANY_OF,['som','pca','kmeans',
      'hierarchical'],0),
    Constraint("distance", CAN_BE_ANY_OF, ['correlation','euclidean'], 0),

    Defaults are used only if no other information is available.

    """

    network = bie3.backchain(
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.PrettySignalFile, "preprocess", "mas5"),
        bie3.Attribute(rulebase.ReportFile, "report_type", "cluster"),
        )
    network = bie3.optimize_network(network)
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case20():
    """Result that generates:
    A network with 368 nodes.

    Problem: The network is different from the result using old
    bie3.py.  The output of Module 49 goes to multiple different
    SignalFile.  It should only go to one SignalFile. In that
    SignalFile, the attributes that are not specified in the
    get_illumina_signal module are set to default.

    JC: I believe this is the correct behavior because each one of the
    output files can be generated by the get_illumina_signal Module,
    as the Module has been described.

    The output Data objects of get_illumina_signal have varying values
    for attributes predataset, quantile_normalize, gene_normalize,
    etc.  The get_illumina_signal needs more consequences to describe
    the values of these parameters.  E.g. There should be a
    Consequence that sets gene_normalize to "no".

    Made some changes to address case22, and now get_illumina_signal
    is not generated.  Not sure what is the issue.  Will look again
    after implementation of new Signal files (and removal of
    processing_step attribute).

    """
    #out_data = rulebase.PrettySignalFile.output(
    #    preprocess='illumina', missing_algorithm="zero_fill",
    #    missing_values='no', logged='yes', quantile_norm="yes",
    #    predataset='yes')
    #out_data = rulebase.SignalFile.output(
    #    preprocess='illumina',
    #    format="gct",
    #    missing_values="unknown",
    #    logged="yes",
    #    )
    #    #missing_values="no",
    #    #missing_algorithm="zero_fill",
    #    #quantile_norm="yes")
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case21():
    """AssertionError: Module make_normalize_report requires a
    PrettySignalFile with preprocess=unknown, but user requests it to
    be mas5.

    Problem: I have added the constraint of preprocess for PcaPlot to
    be SAME_AS PrettySignalFile. But for other attributes, it still
    get the error.  Do we need to constraint all the attributes in
    PrettySignalFile and PcaPlot?

    JC: Fixed.  Will accept user constraints now.
    
    """
    network = bie3.backchain(
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.ReportFile,"report_type","normalize_file"),
        bie3.Attribute(rulebase.PrettySignalFile,"preprocess","mas5"),
        bie3.Attribute(rulebase.PrettySignalFile,"contents","test"),
        bie3.Attribute(rulebase.PrettySignalFile,'gene_center',"median"),
        )
    network = bie3.optimize_network(network)
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case22():
    """Result to get: only the node 0 in the network.
    Need to change the priority of the attributes value:
    1. constraint for priority
    2. get from output
    3. user input
    4. default

    JC: I'm not sure this problem will be fixed with a change in the
    priority.  I thought there was another case where PrettySignalFile
    was an internal node?

    This is only generating 1 node in the network because if
    PrettySignalFile gene_order=t_test_p, then transfer will no longer
    be able to generate it.  It requires gene_order=no.

    """
    network = bie3.backchain(
        rulebase.all_modules, rulebase.PrettySignalFile,
        #bie3.Attribute(rulebase.PrettySignalFile,"gene_order","t_test_p"),
       )
    network = bie3.complete_network(network)
    network = bie3.optimize_network(network)
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def main():

    #run_case01()
    #run_case02()
    #run_case03()
    #run_case04()
    #run_case05()
    #run_case06()
    #run_case07()
    #run_case08()
    #run_case09()
    #run_case10()
    #run_case11()
    #run_case12()
    #run_case13()
    #run_case14()
    run_case15()
  #  run_case16()
    #run_case17()

    #run_case18()
    #run_case19()
    run_case20()
    #run_case21()
    #run_case22()


if __name__ == '__main__':
    main()
    #import cProfile; cProfile.run("main()")
