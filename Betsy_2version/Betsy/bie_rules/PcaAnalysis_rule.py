#PcaAnalysis
from Betsy import bie
import SignalFile_rule,SignalFile2_rule,SignalFile1_rule
PcaAnalysis = bie.DataType(
    'PcaAnalysis',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],DEFAULT='no'),
    #bie.Attribute(
    #    format=[ "tdf", "gct"],
    #    DEFAULT="tdf"),
    # Properties of the data.
    bie.Attribute(
        preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        DEFAULT="unknown"),
    bie.Attribute(
        missing_values=["no"],
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
    #bie.Attribute(annotate=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        unique_genes=["no", "average_genes", "high_var", "first_gene"],
        DEFAULT="no"),
    bie.Attribute(
        duplicate_probe=["no", "yes", "closest_probe", "high_var_probe"],
        DEFAULT="no"),
    #bie.Attribute(rename_sample=["no", "yes"], DEFAULT="no"),

    # Unclassified.
    bie.Attribute(num_features=bie.ANYATOM, DEFAULT="all"),
    #bie.Attribute(
    #    gene_order=[
    #        "no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
    #    DEFAULT="no"),
    bie.Attribute(predataset=["no", "yes"], DEFAULT="no"),
    bie.Attribute(platform=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(filter=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(group_fc=bie.ANYATOM, DEFAULT="no"),
    )

PcaPlot = bie.DataType(
    'PcaPlot',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(pca_gene_num=bie.ANYATOM,DEFAULT='500'),
    bie.Attribute(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],DEFAULT='no'),
    bie.Attribute(
        format=[ "tdf", "gct"],
        DEFAULT="tdf"),
    # Properties of the data.
    bie.Attribute(
        preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        DEFAULT="unknown"),
    bie.Attribute(
        missing_values=["no"],
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
    

list_files = [PcaAnalysis,PcaPlot]
all_modules = [
    bie.Module(
        'analyze_samples_pca',
        SignalFile2_rule.SignalFile2(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],format='tdf',logged='yes',
              preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
              quantile_norm=['yes','no'],bfrm_norm=['yes','no'],combat_norm=['yes','no'],
                   shiftscale_norm=['yes','no'],dwd_norm=['yes','no'],gene_center=['mean','median','no','unknown'],
                   gene_normalize=["unknown", "no", "variance", "sum_of_squares"], 
                   missing_algorithm=["none", "median_fill", "zero_fill"],
                   unique_genes=["no", "average_genes", "high_var", "first_gene"],
                   #num_features=bie.ANYATOM,
                   predataset=["no", "yes"],
                   #platform=bie.ANYATOM,
                   #group_fc=bie.ANYATOM,
                   annotate='no'
                   ),
        PcaAnalysis(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
                   logged='yes',
                   quantile_norm=['yes','no'],bfrm_norm=['yes','no'],combat_norm=['yes','no'],
                   shiftscale_norm=['yes','no'],dwd_norm=['yes','no'],
                   gene_center=['mean','median','no','unknown'],
                   gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
                    missing_algorithm=["none", "median_fill", "zero_fill"],
                    unique_genes=["no", "average_genes", "high_var", "first_gene"],
                    num_features=bie.ANYATOM,
                    filter=bie.ANYATOM,
                    predataset=["no", "yes"],
                    platform=bie.ANYATOM,
                    group_fc=bie.ANYATOM,
                    )),

    bie.Module(
        'plot_sample_pca_wo_label',
        PcaAnalysis(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
                   logged='yes',
                   quantile_norm=['yes','no'],bfrm_norm=['yes','no'],combat_norm=['yes','no'],
                   shiftscale_norm=['yes','no'],dwd_norm=['yes','no'],
                   gene_center=['mean','median','no','unknown'],
                   gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
                   missing_algorithm=["none", "median_fill", "zero_fill"],
                    unique_genes=["no", "average_genes", "high_var", "first_gene"],
                    predataset=["no", "yes"],
                    ),
        PcaPlot(contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                "no"],preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
                   logged='yes',
                   quantile_norm=['yes','no'],bfrm_norm=['yes','no'],combat_norm=['yes','no'],
                   shiftscale_norm=['yes','no'],dwd_norm=['yes','no'],
                   gene_center=['mean','median','no','unknown'],
                   gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
                   missing_algorithm=["none", "median_fill", "zero_fill"],
                    unique_genes=["no", "average_genes", "high_var", "first_gene"],
                    num_features=bie.ANYATOM,
                    filter=bie.ANYATOM,
                    predataset=["no", "yes"],
                    platform=bie.ANYATOM,
                    group_fc=bie.ANYATOM)), 
 ]
