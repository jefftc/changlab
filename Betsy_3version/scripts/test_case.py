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

def run_case23():
    """cannot trace back to GeoSeries to generate the
    ExpressionFile and preprocess with illumina.

    Expected a network with the nodes:
    DATA                        MODULE
    GEOSeries               ->  download_geo                 ->
    ExpressionFiles         ->  extract_illumina_idat_files  ->
    IDATFiles               ->  preprocess_illumina          ->
    ILLUFolder              ->  get_illumina_signal          ->
    SignalFile_Postprocess  ->  convert_signal_to_tdf        ->
    SignalFile_Postprocess

    However, we currently only get a network:
    DATA                        MODULE
    SignalFile_Postprocess  ->  check_for_log                ->
    SignalFile_Postprocess

    """
    out_data = rulebase.SignalFile.output(
        preprocess="illumina",
        format="tdf",
        #logged="no",
        logged="yes",
        )

    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case24():
    """Generate a network start from SignalFile_order, cannot trace back
    to SignalFile_Postprocess.
    
    Expected a network with the nodes:
    DATA                               MODULE
    SignalFile_Postprocess       ->   convert_signal_to_tdf          ->
    SignalFile_Postprocess       ->   check_for_log                  ->
    SignalFile_Postprocess       ->   log_signal                     ->
    SignalFile_Postprocess       ->   convert_postprocess_impute     ->
    SignalFile_Impute            ->   fill_missing_with_zeros        ->
    SignalFile_Impute            ->   convert_impute_merge           ->
    SignalFile_Merge             ->   convert_merge_normalize        ->
    SignalFile_Normalize         ->   check_gene_center              ->
    SignalFile_Normalize         ->   check_gene_normalize           ->
    SignalFile_Normalize         ->   convert_normalize_order        ->
    SignalFile_Order,ClassLableFile-> rank_genes_by_sample_ttest     ->
    GeneListFile,SignalFile_Order->   reorder_genes                  ->
    SignalFile_Order             ->   convert_order_annotate         ->
    SignalFile_Annotate          ->   convert_annotate_filter        ->
    SignalFile

    However, we currently get a network:
    DATA                               MODULE
    SignalFile_Order             ->   convert_order_annotate         ->
    SignalFile_Annotate          ->   convert_annotate_filter        ->
    SignalFile
    
    """
    out_data = rulebase.SignalFile.output(
        format="tdf", logged="yes",
        gene_order='t_test_p',
        )

    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)
    
def run_case25():
    """
    cannot trace back to GEOSeries and SignalFile_Postprocess to
    generate SignalFile_Merge with preprocess=aiglent.
    
    Expected a network generate from GeoSeries or SignalFile_Postprocess
    The Data path in the network is like: 
    GEOSeries -> SignalFile_Postprocess ->
    SignalFile_Impute -> SignalFile_Merge -> (plot_actb_line) -> 
    ActPlot
    
    However, we currently get a network with only one node
    Data(ActbPlot, contents='unspecified')
    
    """
    network = bie3.backchain(
        rulebase.all_modules, rulebase.ActbPlot,
        bie3.Attribute(rulebase.SignalFile_Merge, "preprocess", "agilent"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case26():
    '''Test normalize report,
    
    Expected a network generated form SignalFile_Postprocess with
    contents='test' and preprocess="unknown".
    - The make_normalize_report module has 5 input Data nodes, which
      are SignalFile, IntensityPlot,ControlPlot,PcaPlot and ActbPlot.
    - The SignalFile is generated from:
      SignalFile_Postprocess->SignalFile_Impute->SignalFile_Merge->
      SignalFile_Normalize->SignalFile_Order->SignalFile_Annotate->SignalFile
    - The IntensityPlot,ControlPlot, PcaPlot are generated from SignalFile.
    - The ActbPlot is generated from SignalFile_Merge.

    However, we got a network which has three SignalFile_Postprocess
    with different values for "contents".
     
    Also the network has ExpressionFiles, AgilentFiles,GPRFiles,
    which lead the SignalFile has different "contents" values and
    "preprocess" values.
     
    The IntensityPlot, ControlPlot, PcaPlot and ActvPlot are not
    generated from any other Data.
     
    '''
    network = bie3.backchain(
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.ReportFile,"report_type","normalize_file"),
        bie3.Attribute(rulebase.SignalFile,"preprocess","unknown"),
        bie3.Attribute(rulebase.SignalFile_Merge,"preprocess","unknown"),
        bie3.Attribute(rulebase.SignalFile,"contents","test"),
        #bie3.Attribute(rulebase.SignalFile,"quantile_norm","yes")
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case27():
    """Three problems:
      1. Cannot trace back to GeoSeries to preprocess mas5.
      2. HeapMap node is not generated from ClusterFile.
      3. We require the cluster_alg is pca but different cluster_alg is shown.

      We expected a network generated from GEOSeries and go through to SignalFile,
      GEOSeries -> download_geo -> ExpressionFile ->......-> SignalFile ->
      Cluster_genes_by_pca -> ClusterFile->plot_heatmap -> HeatMap
      ClusterFile, HeatMap -> make_cluster_report->ReportFile

      However, we got a network which is from ExpressionFile, but not GEOSeries.
      The SignalFile can go to different cluster_alg but not the only one we specify.
      HeatMap is isolated from ClusterFile.
    """
    network = bie3.backchain(
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.SignalFile, "preprocess", "mas5"),
        bie3.Attribute(rulebase.ReportFile, "report_type", "cluster"),
        bie3.Attribute(rulebase.SignalFile, "quantile_norm", "yes"),
        bie3.Attribute(rulebase.ClusterFile, "cluster_alg", "pca"),
        bie3.Attribute(rulebase.Heatmap, "cluster_alg", "pca"),
        )
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case28():
    """get error when running this command

    File "test_case.py", line 704, in run_case28
    network = bie3.backchain(rulebase.all_modules, out_data)
    File "/home/xchen/chencode/Betsy_3version/Betsy/bie3.py", line 928, in backchain
    modules = _backchain_to_modules(moduledb, node, user_attributes)
    File "/home/xchen/chencode/Betsy_3version/Betsy/bie3.py", line 1872, in _backchain_to_modules
    if _can_module_produce_data(module, data, user_attributes):
    File "/home/xchen/chencode/Betsy_3version/Betsy/bie3.py", line 2533, in _can_module_produce_data
    if x.name == conseq2.name and x.input_index == const2.arg1]
    NameError: global name 'conseq2' is not defined

    """
    out_data = rulebase.SignalFile.output(group_fc='yes')
    network = bie3.backchain(rulebase.all_modules, out_data)
    network = bie3.optimize_network(network)    
    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)

def run_case29():
    """Expected a network generated from GEOSeries
    
    the path of Data is:
    GEOSeries -> ... -> SignalFile_Merge -> ...
    ->SignalFile_Order -> SignalFile
    SignalFile_Merge, ClassLabelFile -> (convert_label_cls)->ClassLabelFile
    SignalFile_Order,ClassLabelFile-> (rank_genes_sample_ttest)->GeneListFile
    SignalFile_Order,GeneListFile -> (reorder_genes)->SignalFile_Order

    However, we got a network which the input to convert_label_cls is
    not generated from GEOSeries, it is generated from
    SignalFile_Postprocess with preprocess=unknown, That is, we expect
    the node 17 and node 54 to be the same node

    JC: Node 17 and 54 have different preprocess.  In principle, we
    could generate a cls from SignalFiles with different
    preprocessing.  I think the issue is that node 17 should point to
    node 73.

    17  SignalFile_Merge  preprocess="illumina"
    54  SignalFile_Merge  preprocess="unknown"
    45  ClassLabelFile
    73  convert_label_to_cls

    Before optimization, ClassLabelFile (118) + SignalFile_Merge (25)
    should go into convert_label_to_cls (116, 117).

    """
    out_data = rulebase.SignalFile.output(
        preprocess="illumina",
        format="tdf", logged="yes",
         gene_order='t_test_p',
        )

    network = bie3.backchain(rulebase.all_modules, out_data)
    #bie3.write_network("out.network", network)
    #network = bie3.read_network("out.network")
    network = bie3.complete_network(network)
    network = bie3.optimize_network(network)

    bie3.print_network(network)
    bie3.plot_network_gv("out.png", network)


def run_case30():
    '''test normalize report,
    
    Expect a network generated from GEOSeries, the
    make_normalize_report has 6 input data:
    SignalFile,IntensityPlot,ControlPlot,PcaPlot,ActbPlot,PcaPlot.

    The first PcaPlot is generated from SignalFile and we require
    the attribute of quantile_norm, combat_norm,shiftscale_norm, bfrm_norm,
    dwd_norm, gene_center,gene_normalize,unique_genes,platform,group_fc, num_features,
    and duplicate_probes for both SignalFile and first PcaPlot are the same. If it is not
    specified by the user in the output, the value of these attributes will be set
    to output default.
    
    The second PcaPlot is generated from SignalFile and we require
    the attributes of quantile_norm, combat_norm,shiftscale_norm, bfrm_norm,
    dwd_norm, gene_center,gene_normalize,unique_genes,platform,group_fc, num_features,
    and duplicate_probes all set to 'no'.

    The reason of two PcaPlot is that we want to compare the SignalFile
    before any normalization and after normalization.
     
    However, the network we currently got is:
    the attributes of SignalFile, which we are not specified in the output,
    can set to different values,like:
    bfrm_norm=['yes', 'no']
    combat_norm=['yes', 'no']
    dwd_norm=['yes', 'no'] 
    gene_normalize=['variance', 'sum_of_squares', 'no'],
    group_fc=['yes', 'no'],  
    num_features=['yes', 'no'],
    platform=['yes', 'no'],
    shiftscale_norm=['yes', 'no'],
    unique_genes=['average_genes', 'high_var', 'first_gene'])
    duplicate_probes=["no", "closest_probe", "high_var_probe"]
    
    
    '''
    network = bie3.backchain(  
        rulebase.all_modules, rulebase.ReportFile,
        bie3.Attribute(rulebase.ReportFile,"report_type","normalize_file"),
        bie3.Attribute(rulebase.SignalFile,"preprocess","mas5"),
        bie3.Attribute(rulebase.SignalFile,"contents","test"),
        bie3.Attribute(rulebase.SignalFile,"quantile_norm","yes"),
        bie3.Attribute(rulebase.SignalFile,'gene_center',"median"),
        )
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
    #run_case15()
    #run_case16()
    #run_case17()

    #run_case18()
    #run_case19()
    #run_case20()
    #run_case21()
    #run_case22()
    #run_case23()
    #run_case24()
    #run_case25()
    #run_case26()
    #run_case27()
    #run_case28()
    #run_case29()
    run_case30()
    
    
if __name__ == '__main__':
    main()
    #import cProfile; cProfile.run("main()")
