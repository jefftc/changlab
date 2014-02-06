#SignalFile2
from Betsy.bie3 import *
import SignalFile_rule
import SignalFile1_rule

SignalFile2 = DataType(
    "SignalFile2",
    # Properties of the format.
    AttributeDef("format", ["tdf", "gct"], "tdf", "tdf"),
    
    # Properties of the data.
    AttributeDef(
        "preprocess",
        ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        "unknown", "unknown"),
    AttributeDef(
        "missing_values", ["no"], "no", "no"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill"),
    AttributeDef("filter", ["yes", "no"], "no", "no"),
    AttributeDef(
        "logged", ["no", "yes"], "yes", "yes"),
    
    # Normalizing the genes.
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "unknown", "no"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "unknown", "no"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no"),

    # Annotations.
    AttributeDef("annotate", ["no", "yes"], "no", "no"),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no"),
    AttributeDef(
        "duplicate_probe", ["no", "yes", "closest_probe", "high_var_probe"],
        "no", "no"),
    AttributeDef("rename_sample", ["no", "yes"], "no","no"),

    # Unclassified.
    AttributeDef("num_features", ["yes","no"], "no", "no"),
    AttributeDef(
        "gene_order",
        ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
       "no", "no"),
    AttributeDef("predataset", ["no", "yes"], "no", "no"),
    AttributeDef("platform", ["yes","no"], "no", "no"),
    AttributeDef("group_fc", ["yes","no"], "no","no"),
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified")
    )

IntensityPlot = DataType(
    'IntensityPlot',
    AttributeDef("pca_gene_num", ["no", "yes"], "no","no"),
    AttributeDef(
        "format", [ "tdf", "gct"],"tdf","tdf"),
    # Properties of the data.
    AttributeDef(
        "preprocess", ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        "unknown","unknown"),
    AttributeDef(
        "missing_values", ["no"],"no","no"),
    AttributeDef(
        "missing_algorithm", ["none", "median_fill", "zero_fill"],
        "none","none"),
    AttributeDef("filter", ["yes","no"], "no","no"),
    AttributeDef(
        "logged", [ "no", "yes"],"yes","yes"),
    # Normalizing the genes.
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "no","no"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "no","no"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    AttributeDef(
        "dwd_norm", ["no", "yes"], "no","no"),
    AttributeDef(
        "bfrm_norm", ["no", "yes"], "no","no"),
    AttributeDef(
        "quantile_norm", ["no", "yes"], "no","no"),
    AttributeDef(
        "shiftscale_norm", ["no", "yes"], "no","no"),
    AttributeDef(
        "combat_norm", ["no", "yes"], "no","no"),

    # Annotations.
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no"),
    AttributeDef(
        "duplicate_probe", ["no", "yes", "closest_probe", "high_var_probe"],
        "no","no"),
    # Unclassified.
    AttributeDef(
        "gene_order", [
            "no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
        "no","no"),
    AttributeDef("predataset", ["no", "yes"], "no","no"),
    AttributeDef("platform", ["yes","no"], "no","no"),
    AttributeDef("num_features", ["yes","no"], "no","no"),
    AttributeDef("group_fc", ["yes","no"], "no", "no"),
    AttributeDef(
        "contents",
        ["train0", "train1", "test", "class0,class1,test",
         "class0", "class1", "class0,class1", "unspecified"],
        "unspecified", "unspecified")
    )
    
ActbPlot = DataType(
    'ActbPlot',
    # Why is this necessary?
    AttributeDef(
        "preprocess",
        ["unknown", "agilent", "mas5", "rma", "loess", 'illumina'],
        'unknown', 'unknown'),
    )
Hyb_barPlot = DataType('Hyb_barPlot')
ControlPlot = DataType(
    'ControlPlot',
    AttributeDef(
        "preprocess", ["unknown", "agilent", "mas5", "rma", "loess"],
        'unknown', 'unknown'))
BiotinPlot = DataType('BiotinPlot')
HousekeepingPlot = DataType('HousekeepingPlot')

list_files=[
    SignalFile2,
    IntensityPlot,
    ActbPlot,
    Hyb_barPlot,
    BiotinPlot,
    HousekeepingPlot,
    ControlPlot,
    ]

all_modules = [
    Module(
        "transfer2",
        SignalFile1_rule.SignalFile1, SignalFile2,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("missing_values", MUST_BE,"no"),
        Constraint(
            "preprocess",
            CAN_BE_ANY_OF,
            ["unknown", "illumina", "agilent", "mas5", "rma", "loess"]),
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0", "class1",
             "class0,class1","class0,class1,test", "unspecified"]),
        Constraint(
            "gene_center", CAN_BE_ANY_OF, ["unknown", "no", "mean", "median"]),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["unknown", "no", "variance", "sum_of_squares"]),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"]),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("rename_sample", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("combat_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence("contents", SAME_AS_CONSTRAINT, 0),
        Consequence("format", SAME_AS_CONSTRAINT, 0),
        Consequence("logged", SAME_AS_CONSTRAINT, 0),
        Consequence("missing_values", SAME_AS_CONSTRAINT, 0),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT, 0),
        Consequence("preprocess", SAME_AS_CONSTRAINT, 0),
        Consequence("predataset", SAME_AS_CONSTRAINT, 0),
        Consequence("rename_sample", SAME_AS_CONSTRAINT, 0),
        Consequence("filter", SAME_AS_CONSTRAINT, 0),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("gene_center", SAME_AS_CONSTRAINT, 0),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT, 0),
        #Consequence("num_features", SET_TO, "no"),
        #Consequence("platform", SET_TO, "no"),
        #Consequence("duplicate_probe", SET_TO, "no"),
        #Consequence("unique_genes", SET_TO, "no"),
        #Consequence("group_fc", SET_TO, "no"),
        #Consequence("gene_order", SET_TO, "no"),
        #Consequence("annotate", SET_TO, "no"),
        ),
   Module(
        "merge_files_for_classification",
        [SignalFile2,SignalFile2],SignalFile2,
        Constraint("contents", MUST_BE, "class0,class1", 0),
        Constraint("format", MUST_BE, 'gct', 0),
        Constraint("logged", MUST_BE, "yes", 0),
        Constraint("missing_values", MUST_BE,'no', 0),
        Constraint("contents", MUST_BE, "test", 1),
        Constraint("format", MUST_BE, 'gct', 1),
        Constraint("logged", MUST_BE, "yes", 1),
        Constraint("missing_values", MUST_BE, 'no', 1),
        Consequence("contents", SET_TO, "class0,class1,test"),
        Consequence("format", SAME_AS_CONSTRAINT, 0),
        Consequence("logged", SAME_AS_CONSTRAINT, 0),
        Consequence("missing_values", SAME_AS_CONSTRAINT, 0)),
    Module(  
        "filter_genes_by_fold_change_across_classes",
        [SignalFile_rule.ClassLabelFile,SignalFile2],SignalFile2,
        UserInputDef("group_fc_num"),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("format", MUST_BE,"tdf",1),
        Constraint("logged", MUST_BE,"yes",1),
        Constraint("group_fc", MUST_BE,"no",1),
        Constraint("gene_order", MUST_BE,"no",1),
        Constraint("annotate", MUST_BE,"no",1),
        Constraint("num_features", MUST_BE,"no",1),
        Constraint("platform", MUST_BE,"no",1),
        Constraint("duplicate_probe", MUST_BE,"no",1),
        Constraint("unique_genes", MUST_BE,"no",1),
        Consequence("format",SAME_AS_CONSTRAINT,1),
        Consequence("logged",SAME_AS_CONSTRAINT,1),
        Consequence("gene_order",SAME_AS_CONSTRAINT,1),
        Consequence("annotate",SAME_AS_CONSTRAINT,1),
        Consequence("num_features",SAME_AS_CONSTRAINT,1),
        Consequence("platform",SAME_AS_CONSTRAINT,1),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT,1),
        Consequence("unique_genes",SAME_AS_CONSTRAINT,1),
        Consequence("group_fc",SET_TO,"yes")),

    Module(  
        "rank_genes_by_class_neighbors",
        [SignalFile_rule.ClassLabelFile,SignalFile2],SignalFile_rule.GeneListFile,
        UserInputDef("cn_num_neighbors",50),
        UserInputDef("cn_num_perm",100),
        UserInputDef("cn_user_pval",0.5),
        UserInputDef("cn_min_threshold",10),
        UserInputDef("cn_max_threshold",16000),
        UserInputDef("cn_min_folddiff",5),
        UserInputDef("cn_abs_diff",50),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("format", MUST_BE,"tdf",1),
        Constraint("logged", MUST_BE,"yes",1),
        Constraint("gene_order", MUST_BE,"no",1),
        Constraint("annotate", MUST_BE,"no",1),
        Constraint("num_features", MUST_BE,"no",1),
        Constraint("platform", MUST_BE,"no",1),
        Constraint("duplicate_probe", MUST_BE,"no",1),
        Constraint("unique_genes", MUST_BE,"no",1),
        Consequence("gene_order",SET_TO,"class_neighbors"),
        Consequence("cn_mean_or_median",SET_TO_ONE_OF, ['mean', 'median']),
        Consequence("cn_ttest_or_snr",SET_TO_ONE_OF, ['t_test','snr']),
        Consequence("cn_filter_data",SET_TO_ONE_OF, ['yes','no'])),
    Module(
         "rank_genes_by_sample_ttest",
         [SignalFile_rule.ClassLabelFile,SignalFile2],SignalFile_rule.GeneListFile,
        UserInputDef("gene_select_threshold",0.05),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("format", MUST_BE,"tdf",1),
        #Constraint("preprocess", CAN_BE_ANY_OF, ["unknown", "illumina", "agilent",
        #                             "mas5", "rma", "loess"],1),
        Constraint("logged", MUST_BE,"yes",1),
        Constraint("gene_order", MUST_BE,"no",1),
        Constraint("annotate", MUST_BE,"no",1),
        Constraint("num_features", MUST_BE,"no",1),
        Constraint("platform", MUST_BE,"no",1),
        Constraint("duplicate_probe", MUST_BE,"no",1),
        Constraint("unique_genes", MUST_BE,"no",1),
        #Consequence("preprocess",SAME_AS_CONSTRAINT,1),
        Consequence("gene_order",SET_TO_ONE_OF, ["t_test_p", "t_test_fdr"])),
         
    Module(
         "reorder_genes",  
         [SignalFile_rule.GeneListFile,SignalFile2], SignalFile2,
         Constraint("gene_order", CAN_BE_ANY_OF, ['t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list"],0),
         Constraint("format", MUST_BE,"tdf",1),
         Constraint("logged", MUST_BE,"yes",1),
         Constraint("gene_order", MUST_BE,"no",1),
         Constraint("annotate", MUST_BE,"no",1),
         Constraint("num_features", MUST_BE,"no",1),
         Constraint("platform", MUST_BE,"no",1),
         Constraint("duplicate_probe", MUST_BE,'no',1),
         Constraint("unique_genes", MUST_BE,'no',1),
         Constraint(
            "preprocess", CAN_BE_ANY_OF,
            ["unknown", "illumina", "agilent", "mas5", "rma", "loess"], 1),
         Consequence("gene_order",SET_TO_ONE_OF, ['t_test_p', "t_test_fdr",
                         'class_neighbors', "gene_list"]),
         Consequence("format",SAME_AS_CONSTRAINT,1),
         Consequence("logged",SAME_AS_CONSTRAINT,1),
         Consequence("annotate",SAME_AS_CONSTRAINT,1),
         Consequence("num_features",SAME_AS_CONSTRAINT,1),
         Consequence("platform",SAME_AS_CONSTRAINT,1),
         Consequence("duplicate_probe",SAME_AS_CONSTRAINT,1),   
         Consequence("unique_genes",SAME_AS_CONSTRAINT,1),
         Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
        ),
    Module(
         'annotate_probes',
         SignalFile2,SignalFile2,
         Constraint("format", MUST_BE,"tdf"),
         Constraint("logged", MUST_BE,"yes"),
         Constraint("annotate", MUST_BE,"no"),
         Constraint("num_features", MUST_BE,"no"),
         Constraint("platform", MUST_BE,"no"),
         Constraint("duplicate_probe", MUST_BE,'no'),
         Constraint("unique_genes", MUST_BE,'no'),
         Consequence("annotate",SET_TO,"yes"),
         Consequence("format",SAME_AS_CONSTRAINT),
         Consequence("logged",SAME_AS_CONSTRAINT),
         Consequence("num_features",SAME_AS_CONSTRAINT),
         Consequence("platform",SAME_AS_CONSTRAINT),
         Consequence("duplicate_probe",SAME_AS_CONSTRAINT),
         Consequence("unique_genes",SAME_AS_CONSTRAINT)),
    Module(
        'remove_duplicate_genes',
        SignalFile2, SignalFile2,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("annotate", MUST_BE,"yes"),
        Constraint("num_features", MUST_BE,"no"),
        Constraint("platform", MUST_BE,"no"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Constraint("unique_genes", MUST_BE,'no'),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("platform", SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
        Consequence(
            "unique_genes", SET_TO_ONE_OF,
            ['average_genes', 'high_var', 'first_gene'])),
      
    Module(
         'select_first_n_genes',
        SignalFile2,SignalFile2,
        UserInputDef("num_features_value",500),
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("num_features", MUST_BE,"no"),
        Constraint("platform", MUST_BE,"no"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("num_features",SET_TO,"yes"),
        Consequence("platform",SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT)), 
     Module(
         'add_crossplatform_probeid',
         SignalFile2,SignalFile2,
         UserInputDef("platform_name"),
         Constraint("format", MUST_BE,"tdf"),
         Constraint("logged", MUST_BE,"yes"),
         Constraint("platform", MUST_BE,"no"),
         Constraint("duplicate_probe", MUST_BE,'no'),
         Consequence("format",SAME_AS_CONSTRAINT),
         Consequence("logged",SAME_AS_CONSTRAINT),
         Consequence("platform",SET_TO,"yes"),
         Consequence("duplicate_probe",SAME_AS_CONSTRAINT)), 
     Module(
        'remove_duplicate_probes',
        SignalFile2,SignalFile2,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe",SET_TO,'high_var_probe')),
    Module(
         'select_probe_by_best_match',
        SignalFile2,SignalFile2,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe",SET_TO,'closest_probe')),
   Module( 
        'convert_signal_to_gct',
        SignalFile2,SignalFile2,
        Constraint("format", MUST_BE,"tdf"),
        Consequence("format",SET_TO,"gct")),
   Module( 
        'unlog_signal',
        SignalFile2,SignalFile2,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SET_TO,"no")),

    Module(    
        'plot_intensity_boxplot',
        SignalFile2, IntensityPlot,
        Constraint(
            "contents", CAN_BE_ANY_OF,
            ["train0", "train1", "test", "class0,class1,test",
             "class0", "class1", "class0,class1", "unspecified"]),
        Constraint("format", MUST_BE, 'tdf'),
        Constraint("logged", MUST_BE, 'yes'),
        Constraint("filter", CAN_BE_ANY_OF, ['yes',"no"]),
        Constraint(
            "preprocess", CAN_BE_ANY_OF,
            ["unknown", "illumina", "agilent", "mas5", "rma", "loess"]),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("combat_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ['yes','no']),
        Constraint(
            "gene_center", CAN_BE_ANY_OF, ['mean','median','no','unknown']),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            ["unknown", "no", "variance", "sum_of_squares"]),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"]),
        Constraint(
            "unique_genes", CAN_BE_ANY_OF,
            ["no", "average_genes", "high_var", "first_gene"]),
        Constraint("num_features", CAN_BE_ANY_OF, ["yes","no"]),
        Constraint("predataset", CAN_BE_ANY_OF, ["yes","no"]),
        Constraint("platform", CAN_BE_ANY_OF, ["yes","no"]),
        Constraint("group_fc", CAN_BE_ANY_OF, ["yes","no"]),
        Constraint("annotate", MUST_BE, 'no'),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm", SAME_AS_CONSTRAINT),
        Consequence("unique_genes", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("platform", SAME_AS_CONSTRAINT),
        Consequence("group_fc", SAME_AS_CONSTRAINT),
        ),
                    
    Module(
        'plot_actb_line',
        # Why SignalFile and not SignalFile2?
        SignalFile_rule.SignalFile, ActbPlot,
        Constraint(
            'preprocess', CAN_BE_ANY_OF,
            ["unknown", "agilent","illumina", "mas5", "rma", "loess"]),
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("missing_values", MUST_BE,"no"),
        Constraint("quantile_norm", MUST_BE,"no"),
        Constraint("combat_norm", MUST_BE,"no"),
        Constraint("shiftscale_norm", MUST_BE,"no"),
        Constraint("dwd_norm", MUST_BE,"no"),
        Constraint("bfrm_norm", MUST_BE,"no"),
        Consequence('preprocess',SAME_AS_CONSTRAINT)),
        
    Module(
        'plot_affy_affx_line',
        SignalFile2,ControlPlot,
        Constraint("annotate", MUST_BE,'no'),
        Constraint("preprocess", CAN_BE_ANY_OF, ["unknown", "agilent",
                                "mas5", "rma", "loess"]),
        Consequence("preprocess",SAME_AS_CONSTRAINT)),
    Module(
       'plot_illu_hyb_bar',
       SignalFile_rule.ControlFile,Hyb_barPlot,
       Constraint("preprocess", MUST_BE,"illumina"),
       Constraint("format", MUST_BE,"gct"),
       Constraint("logged", MUST_BE,"no")),
    Module(
        'plot_illu_biotin_line',
        SignalFile_rule.ControlFile,BiotinPlot,
        Constraint("preprocess", MUST_BE,'illumina'),
        Constraint("format", MUST_BE,"gct"),
        Constraint("logged", MUST_BE,"no")),
       
    Module(
        'plot_illu_housekeeping_line',
        SignalFile_rule.ControlFile,HousekeepingPlot,
        Constraint("preprocess", MUST_BE,'illumina'),
        Constraint("format", MUST_BE,"gct"),
        Constraint("logged", MUST_BE,"no")),
    
    ]
