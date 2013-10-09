#SignalFile2
import bie
import SignalFile_rule
import SignalFile1_rule

SignalFile2 = bie.DataType(
    "SignalFile2",
    # Properties of the format.
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(
        format=[ "tdf", "gct"],
        DEFAULT="tdf"),
    # Properties of the data.
    bie.Attribute(
        preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        DEFAULT="unknown"),
    bie.Attribute(
        missing_values=[ "no"],
        DEFAULT="no"),
    bie.Attribute(
        missing_algorithm=["none", "median_fill", "zero_fill"],
        DEFAULT="none", OPTIONAL=True),
    bie.Attribute(
        logged=[ "no", "yes"],
        DEFAULT="yes"),
    # Normalizing the genes.
    bie.Attribute(
        gene_center=["unknown", "no", "mean", "median"],
        DEFAULT="unknown"),
    bie.Attribute(
        gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
        DEFAULT="unknown"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    bie.Attribute(
        dwd_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        bfrm_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        quantile_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        shiftscale_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        combat_norm=["no", "yes"], DEFAULT="no"),

    # Annotations.
    bie.Attribute(annotate=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        unique_genes=["no", "average_genes", "high_var", "first_gene"],
        DEFAULT="no"),
    bie.Attribute(
        duplicate_probe=["no", "yes", "closest_probe", "high_var_probe"],
        DEFAULT="no"),
    bie.Attribute(rename_sample=["no", "yes"], DEFAULT="no"),


    # Unclassified.
    bie.Attribute(num_features=bie.ANYATOM, DEFAULT="all"),
    bie.Attribute(
        gene_order=[
            "no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
        DEFAULT="no"),
    bie.Attribute(predataset=["no", "yes"], DEFAULT="no"),
    bie.Attribute(platform=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(filter=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(group_fc=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(contents=["train0", "train1", "test", "class0,class1,test",
                        "class0", "class1", "class0,class1", "no"], DEFAULT="no")
    )
list_files=[SignalFile2]
all_modules = [
    bie.QueryModule(
        "transfer2",
        SignalFile1_rule.SignalFile1(
            format="tdf", logged="yes",
            missing_values="no",preprocess=["unknown", "illumina", "agilent",
                                "mas5", "rma", "loess"],
            contents=["train0", "train1", "test",
                        "class0", "class1", "class0,class1","class0,class1,test", "no"]),
        SignalFile2(
            format="tdf", logged="yes",
            missing_values="no", missing_algorithm=["none", "median_fill", "zero_fill"],          
                    preprocess=["unknown", "illumina", "agilent",
                             "mas5", "rma", "loess"],
                    predataset=["no", "yes"], rename_sample=["no", "yes"],
                    filter=bie.ANYATOM,
                    dwd_norm=["no", "yes"], bfrm_norm=["no", "yes"],
                    quantile_norm=["no", "yes"],
                    shiftscale_norm=["no", "yes"], combat_norm=["no", "yes"],
                    group_fc='no', gene_order='no', annotate="no",
                    num_features="all", platform="no",
                    duplicate_probe='no', unique_genes='no',
                    gene_center=["unknown", "no", "mean", "median"],
                    gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
                    contents=["train0", "train1", "test",
                        "class0", "class1", "class0,class1","class0,class1,test", "no"],
            )),
   bie.Module(
        "merge_files_for_classification",
        [SignalFile2(contents="class0,class1",format='gct',logged="yes",
                       missing_values='no'),
        SignalFile2(contents="test",format='gct',logged="yes",
                       missing_values='no')],
        SignalFile2(contents="class0,class1,test",format='gct',logged="yes",
                       missing_values='no')),
    bie.Module(  
        "filter_genes_by_fold_change_across_classes",
        [SignalFile_rule.ClassLabelFile,
         SignalFile2(
             format="tdf", logged="yes",
             group_fc="no", gene_order='no', annotate="no",
             num_features="all", platform="no",
            duplicate_probe='no', unique_genes='no')
         ],
        SignalFile2(
            format="tdf", logged="yes",
            group_fc=bie.ANYATOM, gene_order='no', annotate="no",
            num_features="all", platform="no",
            duplicate_probe='no', unique_genes='no')),

    bie.Module(  
        "rank_genes_by_class_neighbors",
        [SignalFile2(
            format="tdf", logged="yes",
            gene_order="no",annotate="no",
            num_features="all", platform="no",
            duplicate_probe='no', unique_genes='no'),
         SignalFile_rule.ClassLabelFile],
        SignalFile_rule.GeneListFile(gene_order="class_neighbors",
                    cn_num_neighbors=bie.ANYATOM,
                    cn_num_perm=bie.ANYATOM,
                    cn_user_pval=bie.ANYATOM,  
                    cn_mean_or_median=['mean', 'median'],
                    cn_ttest_or_snr=['t_test','snr'],
                    cn_filter_data=['yes','no'],
                    cn_min_threshold=bie.ANYATOM, 
                    cn_max_threshold=bie.ANYATOM, 
                    cn_min_folddiff=bie.ANYATOM, 
                    cn_abs_diff=bie.ANYATOM)),
    bie.Module(
         "rank_genes_by_sample_ttest",   
         [SignalFile2(format="tdf", logged="yes",
                   gene_order="no", annotate="no",
                    num_features="all", platform="no",
                    duplicate_probe='no', unique_genes='no'),
          SignalFile_rule.ClassLabelFile],
         SignalFile_rule.GeneListFile(gene_order=["t_test_p", "t_test_fdr"],
                      cn_num_neighbors=bie.ANYATOM,
                    cn_num_perm=bie.ANYATOM,
                    cn_user_pval=bie.ANYATOM,  
                    cn_mean_or_median=['mean', 'median'],
                    cn_ttest_or_snr=['t_test','snr'],
                    cn_filter_data=['yes','no'],
                    cn_min_threshold=bie.ANYATOM, 
                    cn_max_threshold=bie.ANYATOM, 
                    cn_min_folddiff=bie.ANYATOM, 
                    cn_abs_diff=bie.ANYATOM)),
    bie.Module(
         "reorder_genes",  
         [SignalFile2(format="tdf", logged="yes",
                     gene_order="no", annotate="no",
                     num_features="all", platform="no",
                     duplicate_probe='no', unique_genes='no'),
          SignalFile_rule.GeneListFile(gene_order=['t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list"])],
         SignalFile2(
             format="tdf", logged="yes", annotate="no",
             num_features="all", platform="no",
            duplicate_probe='no', unique_genes='no',
             gene_order=['t_test_p', "t_test_fdr",
                         'class_neighbors', "gene_list"])),
    bie.Module(
         'annotate_probes',
         SignalFile2(format="tdf", logged="yes", annotate="no",
                    num_features="all", platform="no",
                    duplicate_probe='no', unique_genes='no'),
         SignalFile2(format="tdf", logged="yes", annotate="yes",
                    num_features="all", platform="no",
                    duplicate_probe='no', unique_genes='no')),
    bie.Module(
        'remove_duplicate_genes',
        SignalFile2(format="tdf",  logged="yes", annotate="yes", unique_genes="no",
                   platform="no", duplicate_probe='no', num_features="all",),
        SignalFile2(format="tdf",  logged="yes", annotate="yes",
                   unique_genes=['average_genes', 'high_var', 'first_gene'],
                   platform="no", duplicate_probe='no', num_features="all")),
    bie.Module(
         'select_first_n_genes',
         SignalFile2(format="tdf", logged="yes", num_features="all", platform="no",
                    duplicate_probe='no'),
         SignalFile2(format="tdf", logged="yes", num_features=bie.ANYATOM, platform="no",
                    duplicate_probe='no')), 
     bie.Module(
         'add_crossplatform_probeid',
         SignalFile2(format="tdf", logged="yes", platform="no",
                    duplicate_probe='no'),
         SignalFile2(format="tdf", logged="yes", platform=bie.ANYATOM,
                    duplicate_probe='no')),
     bie.Module(
        'remove_duplicate_probes',
        SignalFile2(format="tdf", logged="yes", duplicate_probe='no'),
        SignalFile2(format="tdf", logged="yes", duplicate_probe='high_var_probe')),
    bie.Module(
         'select_probe_by_best_match',
         SignalFile2(format="tdf", logged="yes", duplicate_probe='no'),
         SignalFile2(format="tdf", logged="yes", duplicate_probe='closest_probe')),
    bie.Module( 
        'convert_signal_to_gct',
        SignalFile2(format="tdf"),
        SignalFile2(format="gct")),
    bie.Module( 
        'unlog_signal',
        SignalFile2(format="tdf", logged="yes"),
        SignalFile2(format="tdf", logged="no")),
    
   
    ]
