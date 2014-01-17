#SignalFile1
from Betsy import bie
import SignalFile_rule
SignalFile1 = bie.DataType(
    "SignalFile1",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    # Properties of the format.
    bie.Attribute(
        format=[ "tdf", "pcl"],
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
        logged=["yes"],
        DEFAULT="yes"),
    bie.Attribute(predataset=["no", "yes"], DEFAULT="no"),
    bie.Attribute(filter=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(rename_sample=["no", "yes"], DEFAULT="no"),
    bie.Attribute(dwd_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(bfrm_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(quantile_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(shiftscale_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(combat_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        gene_center=["unknown", "no", "mean", "median"],
        DEFAULT="unknown"),
    bie.Attribute(
        gene_normalize=["unknown", "no", "variance", "sum_of_squares"],
        DEFAULT="unknown"),
    bie.Attribute(contents=["train0", "train1", "test", 'class0,class1,test',
                        "class0", "class1", "class0,class1", "no"], DEFAULT="no")
    )
list_files = [SignalFile1]
all_modules = [
    bie.QueryModule(
        "transfer1",
        SignalFile_rule.SignalFile(
            format="tdf", logged="yes",
            missing_values="no",
            preprocess=["unknown", "illumina", "agilent",
                                "mas5", "rma", "loess"],
            rename_sample=["no", "yes"],
            missing_algorithm=["none", "median_fill", "zero_fill"],
            predataset=["no", "yes"], 
            dwd_norm=["no", "yes"], bfrm_norm=["no", "yes"],
            quantile_norm=["no", "yes"],
            shiftscale_norm=["no", "yes"], combat_norm=["no", "yes"],
            contents=["train0", "train1", "test", 'class0,class1,test',
                        "class0", "class1", "class0,class1", "no"]),
        SignalFile1(
            format="tdf", logged="yes",
            missing_values="no",
            missing_algorithm=["none", "median_fill", "zero_fill"],           
            preprocess=["unknown", "illumina", "agilent",
                                "mas5", "rma", "loess"],
            predataset=["no", "yes"], rename_sample=["no", "yes"],
            filter=bie.ANYATOM,
            dwd_norm=["no", "yes"], bfrm_norm=["no", "yes"],
            quantile_norm=["no", "yes"],
            shiftscale_norm=["no", "yes"], combat_norm=["no", "yes"],
            gene_center='unknown',gene_normalize='unknown',
            contents=["train0", "train1", "test", 'class0,class1,test',
                        "class0", "class1", "class0,class1", "no"]
            )),
    bie.QueryModule(
        "check_gene_center",
        SignalFile1(
            format="tdf", logged="yes", missing_values="no",
            gene_center="unknown", gene_normalize='unknown'),
        SignalFile1(
            format="tdf", logged="yes", missing_values="no",
            gene_center=["no", "mean", "median"],gene_normalize='unknown')),
    bie.QueryModule(
        "check_gene_normalize",
        SignalFile1(
            format="tdf", logged="yes", missing_values="no",
            gene_normalize="unknown"),
        SignalFile1(
            format="tdf", logged="yes", missing_values="no",
            gene_normalize=["no", "variance", "sum_of_squares"])),
    bie.Module(   
        "convert_signal_to_pcl",
        SignalFile1(format='tdf'),
        SignalFile1(format='pcl')),
    bie.Module(
        "center_genes",
        SignalFile1(
            format="pcl", logged="yes", missing_values="no",
            gene_center="no",gene_normalize='unknown'),
        SignalFile1(
            format="tdf", logged="yes", missing_values="no",
            gene_center=["mean", "median"],gene_normalize='unknown')),
    bie.Module(
        "normalize_genes",
        SignalFile1(
            format="pcl", logged="yes", missing_values="no",
            gene_normalize="no"),
        SignalFile1(
            format="tdf", logged="yes", missing_values="no",
            gene_normalize=["variance", "sum_of_squares"]))
    ]
