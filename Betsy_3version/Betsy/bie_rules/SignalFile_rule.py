# SignalFile
#from Betsy import bie3
from Betsy.bie3 import *

ILLU_MANIFEST = [
    'HumanHT-12_V3_0_R2_11283641_A.txt',
    'HumanHT-12_V4_0_R2_15002873_B.txt',
    'HumanHT-12_V3_0_R3_11283641_A.txt',
    'HumanHT-12_V4_0_R1_15002873_B.txt',
    'HumanMI_V1_R2_XS0000122-MAP.txt',
    'HumanMI_V2_R0_XS0000124-MAP.txt',
    'HumanRef-8_V2_0_R4_11223162_A.txt',
    'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
    'HumanRef-8_V3_0_R2_11282963_A.txt',
    'HumanRef-8_V3_0_R3_11282963_A.txt',
    'HumanWG-6_V2_0_R4_11223189_A.txt',
    'HumanWG-6_V3_0_R2_11282955_A.txt',
    'HumanWG-6_V3_0_R3_11282955_A.txt',
    'MouseMI_V1_R2_XS0000127-MAP.txt',
    'MouseMI_V2_R0_XS0000129-MAP.txt',
    'MouseRef-8_V1_1_R4_11234312_A.txt',
    'MouseRef-8_V2_0_R2_11278551_A.txt',
    'MouseRef-8_V2_0_R3_11278551_A.txt',
    'MouseWG-6_V1_1_R4_11234304_A.txt',
    'MouseWG-6_V2_0_R2_11278593_A.txt',
    'MouseWG-6_V2_0_R3_11278593_A.txt',
    'RatRef-12_V1_0_R5_11222119_A.txt'
    ]
ILLU_CHIP = [
    'ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
    'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
    'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
    'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
    'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
    'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
    'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
    'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
    'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
    'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
    'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
    'ilmn_RatRef_12_V1_0_R5_11222119_A.chip'
    ]


GEOSeries = DataType("GEOSeries",
                     AttributeDef("contents",[
                         "train0", "train1","test", "class0,class1,test",
                         "class0", "class1", "class0,class1","unspecified"],
                                  'unspecified','unspecified'))

ExpressionFiles = DataType("ExpressionFiles",
                           AttributeDef("contents",
                                        ["train0","train1", "test",
                                         "class0,class1,test","class0",
                                         "class1", "class0,class1","unspecified"],
                                        'unspecified','unspecified'))

CELFiles = DataType(
    "CELFiles",
    AttributeDef("version", ["unknown", "cc", "v3", "v4"], "unknown", "v3"),
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified')
    )
RenameFile = DataType(
    'RenameFile',
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                             "class1", "class0,class1","unspecified"],
                             'unspecified','unspecified'))
AgilentFiles = DataType(
    "AgilentFiles",
    AttributeDef("contents",["train0","train1", "test",
                            "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                              'unspecified','unspecified'))

ControlFile = DataType(
    "ControlFile",
    AttributeDef(
        'preprocess',["illumina"],
        "illumina","illumina"),
    AttributeDef(
        'missing_values',["unknown", "no", "yes"],
        "no","no"),
    AttributeDef(
        "missing_algorithm",["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill"),
    AttributeDef(
        "logged",["no"],"no","no"),
    AttributeDef(
        'format',["gct"],
        "gct","gct"),
    AttributeDef("contents",["train0","train1", "test",
                            "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                              'unspecified','unspecified'))
    
    
GPRFiles = DataType(
    "GPRFiles",
    AttributeDef("contents",["train0","train1", "test",
                            "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                              'unspecified','unspecified'))

IDATFiles = DataType(
    "IDATFiles",
    AttributeDef("contents",["train0","train1", "test",
                            "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                              'unspecified','unspecified'))

ClassLabelFile = DataType(
    "ClassLabelFile",
    AttributeDef(
    "contents",["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "unspecified"],'unspecified','unspecified'),
    AttributeDef("cls_format",['cls','label','unknown'],"unknown","cls")
    )

ILLUFolder = DataType(
    "ILLUFolder", 
    AttributeDef(
        "illu_manifest",ILLU_MANIFEST,
        'HumanHT-12_V4_0_R2_15002873_B.txt','HumanHT-12_V4_0_R2_15002873_B.txt'),
    AttributeDef(
        'illu_chip',ILLU_CHIP,
        'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip','ilmn_HumanHT_12_V4_0_R1_15002873_B.chip'),
    AttributeDef('illu_bg_mode',['false', 'true'], "false", "false"),
    AttributeDef('illu_coll_mode',['none', 'max', 'median'], "none","none"),
    AttributeDef("contents",["train0","train1", "test",
                            "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                              'unspecified','unspecified')
    )


GeneListFile=DataType(
    "GeneListFile", 
    AttributeDef('cn_mean_or_median',['mean', 'median'], 'mean','mean'),
    AttributeDef('cn_ttest_or_snr',['t_test','snr'], 't_test','t_test'),
    AttributeDef('cn_filter_data',['yes','no'], 'no','no'),
    AttributeDef('gene_order',['no', "gene_list", "class_neighbors",
                          "t_test_p", "t_test_fdr"], 'gene_list',"gene_list"),
    AttributeDef("contents",["train0","train1", "test",
                            "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                              'unspecified','unspecified'))
   
SignalFile = DataType(
    "SignalFile",
    AttributeDef("format", ["unknown", "tdf", "pcl", "gct", "res", "jeffs"],
              "unknown", "tdf"),

    # Properties of the data.
    AttributeDef("preprocess",
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
              "unknown", "unknown"),
    AttributeDef("missing_values", ["unknown", "no", "yes"], "unknown", "no"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill"),
    AttributeDef("logged", ["unknown", "no", "yes"], "unknown", "yes"),
    AttributeDef("filter", ["no", "yes"], "no", "no"),
    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no"),
    
    # Other attributes.
    AttributeDef("predataset", ["no", "yes"], "no", "no"),
    # Normalizing the genes.
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "unknown", "no"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "unknown", "no"),
    AttributeDef("contents", [
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],
              "unspecified", "unspecified"),
    AttributeDef("processing_step",["postprocess","impute","merge",
                                    "normalize","processed"],"postprocess","processed")
##    AttributeDef("processing_step",["postprocess","impute","merge",
##                                    "normalize","order","annotate","filter",
##                                    "processed"],"postprocess","processed")
    )

PrettySignalFile = DataType(
    "PrettySignalFile",
    AttributeDef("format", ["tdf"],
              "tdf", "tdf"),

    # Properties of the data.
    AttributeDef("preprocess",
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
              "unknown", "unknown"),
    AttributeDef("missing_values", ["no"], "no", "no"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill"),
    AttributeDef("logged", [ "yes"], "yes", "yes"),
    AttributeDef("filter", ["no", "yes"], "no", "no"),
    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no"),
    
    # Other attributes.
    AttributeDef("predataset", ["no", "yes"], "no", "no"),
    AttributeDef("rename_sample", ["no", "yes"], "no", "no"),
    # Normalizing the genes.
    AttributeDef(
        "gene_center", [ "no", "mean", "median"],
        "no", "no"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"],
        "no", "no"),
    # Annotations.
    AttributeDef("annotate", ["no", "yes"], "no", "no"),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no"),
    AttributeDef(
        "duplicate_probe", ["no", "yes", "closest_probe", "high_var_probe"],
        "no", "no"),
    # Unclassified.
    AttributeDef("num_features", ["yes","no"], "no", "no"),
    AttributeDef(
        "gene_order",
        ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr"],
       "no", "no"),
    
    AttributeDef("platform", ["yes","no"], "no", "no"),
    AttributeDef("group_fc", ["yes","no"], "no","no"),
    AttributeDef("contents", [
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],
              "unspecified", "unspecified"),
    AttributeDef("processing_step",["order","annotate","filter",
                                    "processed"],"order","processed")
    )
    

all_modules = [
    Module(
        "download_geo", GEOSeries, ExpressionFiles,
         UserInputDef("GSEID"), UserInputDef("GPLID",""),
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
        ),
    #CELFiles
    Module(
        "extract_CEL_files", ExpressionFiles, CELFiles,
         Consequence("version", SET_TO, "unknown"),
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
        ),
    Module(
        "detect_CEL_version",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "unknown"),
        Consequence("version", BASED_ON_DATA, ["cc", "v3", "v4"]),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
        ),
    Module(
        "convert_CEL_cc_to_CEL_v3",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "cc"),
        Consequence("version", SET_TO, "v3"),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
        ),
     # IDATFiles
    Module("extract_illumina_idat_files",
            ExpressionFiles, IDATFiles,
           Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),),
    Module(
        "preprocess_illumina",
        IDATFiles, ILLUFolder,
        Consequence('illu_manifest',SET_TO_ONE_OF,ILLU_MANIFEST),
        Consequence('illu_chip',SET_TO_ONE_OF,ILLU_CHIP),
        Consequence('illu_bg_mode',SET_TO_ONE_OF,["false", "true"]),
        Consequence('illu_coll_mode',SET_TO_ONE_OF,["none", "max", "median"]),
        UserInputDef("illu_clm",''),
        UserInputDef("illu_custom_chip",''),
        UserInputDef("illu_custom_manifest",''),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
        Consequence("contents",SAME_AS_CONSTRAINT),
        ),
       Module(
        "get_illumina_signal",
         ILLUFolder, SignalFile,
         Consequence('preprocess',SET_TO,"illumina"),
         Consequence('format', SET_TO, "gct"),
         Consequence('logged', SET_TO, "no"),
         Consequence('missing_values', SET_TO, "unknown"),
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),),
    Module(
        "get_illumina_control",
         ILLUFolder,ControlFile,
         Consequence('preprocess',SET_TO,"illumina"),
         Consequence("format",SET_TO,"gct"),
         Consequence("logged",SET_TO,"no"),
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT), 
        ),
    
    # AgilentFiles
    Module(
        "extract_agilent_files", ExpressionFiles, AgilentFiles,
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT)),
    Module(
        "preprocess_agilent",
         AgilentFiles,SignalFile,
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence('logged',SET_TO,"unknown"),
        Consequence('preprocess',SET_TO,"agilent"),
        Consequence('format',SET_TO,"tdf")),

    # GPRFiles
    Module(
        "extract_gpr_files", ExpressionFiles, GPRFiles,
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT)),
    Module(
        "normalize_with_loess",
        GPRFiles,SignalFile,
        Consequence("format",SET_TO,"tdf"),
        Consequence("logged",SET_TO,"unknown"),
        Consequence("preprocess",SET_TO,"loess"),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT)),
    
    Module(
        "preprocess_rma",
        CELFiles, SignalFile,
        Constraint("version", CAN_BE_ANY_OF, ["v3", "v4"]),
        Consequence("logged", SET_TO, "yes"),
        Consequence("preprocess", SET_TO, "rma"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence("quantile_norm", SET_TO, "yes"),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
        ),
    Module(
        "preprocess_mas5",
        CELFiles, SignalFile,
        Constraint("version", CAN_BE_ANY_OF, ["v3", "v4"]),
        Consequence("logged", SET_TO, "no"),
        Consequence("preprocess", SET_TO, "mas5"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence("quantile_norm", SET_TO, "no"),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"]),
         Consequence("contents",SAME_AS_CONSTRAINT),
        ),
    ####postprocess
    Module(
        "convert_signal_to_tdf",
        SignalFile, SignalFile,
        Constraint("format", CAN_BE_ANY_OF, ['pcl', 'gct', 'res', 'jeffs']),
        Consequence("format", SET_TO, "tdf"),
        Constraint("processing_step",MUST_BE,"postprocess"),
        Consequence("processing_step",SET_TO_ONE_OF,["postprocess","impute","merge",
                                    "normalize",
                                    "processed"])
        ),
    Module(
        "check_for_log",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", BASED_ON_DATA, ["yes", "no"]),
        Constraint("processing_step",MUST_BE,"postprocess"),
        Consequence("processing_step",SET_TO_ONE_OF,["postprocess","impute","merge",
                                    "normalize",
                                    "processed"])
        ),
    Module(
        "filter_and_threshold_genes",
        SignalFile,SignalFile,
        Constraint('format',MUST_BE,"tdf"),
        Constraint('logged',MUST_BE,"no"),
        Constraint('predataset',MUST_BE,"no"),
        Consequence('format',SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('predataset',SET_TO,'yes'),
        Constraint("processing_step",MUST_BE,"postprocess"),
        Consequence("processing_step",SET_TO_ONE_OF,["postprocess","impute","merge",
                                    "normalize",
                                    "processed"])),
    Module(
        "log_signal",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "yes"),
        Constraint("processing_step",MUST_BE,"postprocess"),
        Consequence("processing_step",SET_TO_ONE_OF,["postprocess","impute","merge",
                                    "normalize",
                                    "processed"])
        ),
    #impute
    Module(
        "check_for_missing_values",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", BASED_ON_DATA, ["no", "yes"]),
        Constraint("processing_step",MUST_BE,"impute"),
        Consequence("processing_step",SET_TO_ONE_OF,["impute","merge",
                                    "normalize",
                                    "processed"])
        ),
    Module(
        "filter_genes_by_missing_values",
        SignalFile, SignalFile,
        UserInputDef("filter_genes_with_missing_values", 0.50),
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "yes"),
        Constraint("filter", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", SAME_AS_CONSTRAINT),
        Consequence("filter", SET_TO, "yes"),
        Constraint("processing_step",MUST_BE,"impute"),
        Consequence("processing_step",SET_TO_ONE_OF,["impute","merge",
                                    "normalize",
                                    "processed"])
        ),
    Module(
        "fill_missing_with_zeros",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "yes"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", SET_TO, "no"),
        Consequence(
            "missing_algorithm", SET_TO, "zero_fill", side_effect=True),
        Constraint("processing_step",MUST_BE,"impute"),
        Consequence("processing_step",SET_TO_ONE_OF,["impute","merge",
                                    "normalize",
                                    "processed"])
        ),
    Module(
        "fill_missing_with_median",
        SignalFile,SignalFile,
        Constraint('format',MUST_BE,'tdf'),
        Constraint('logged',MUST_BE,'yes'),
        Constraint('missing_values',MUST_BE,'yes'),
        Constraint("processing_step",MUST_BE,"impute"),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('missing_algorithm',SET_TO,"median_fill",side_effect=True),
        Consequence('missing_values',SET_TO,'no'),
        Consequence("processing_step",SET_TO_ONE_OF,["impute","merge",
                                    "normalize",
                                    "processed"])),
   
    #merge
    Module(  ###
        "merge_two_classes", [SignalFile, SignalFile], SignalFile,  
         Constraint("contents", MUST_BE, "class0", 0),
         Constraint("format", MUST_BE, "tdf", 0),
         Constraint("logged", MUST_BE, "yes", 0),
         Constraint("missing_values",MUST_BE,'no',0),
         Constraint("combat_norm",MUST_BE,'no',0),
         Constraint("quantile_norm",MUST_BE,"no",0),
         Constraint("dwd_norm",MUST_BE,"no",0),
         Constraint("bfrm_norm",MUST_BE,"no",0),
         Constraint("shiftscale_norm",MUST_BE,"no",0),
         Constraint("contents", MUST_BE, "class1", 1),
         Constraint("format", MUST_BE, "tdf", 1),
         Constraint("logged", MUST_BE, "yes", 1),
         Constraint("missing_values",MUST_BE,'no',1),
         Constraint("combat_norm",MUST_BE,'no',1),
         Constraint("quantile_norm",MUST_BE,"no",1),
         Constraint("dwd_norm",MUST_BE,"no",1),
         Constraint("bfrm_norm",MUST_BE,"no",1),
         Constraint("shiftscale_norm",MUST_BE,"no",1),
         Constraint("predataset",CAN_BE_ANY_OF,["yes","no"],0),###
         Constraint("predataset",CAN_BE_ANY_OF,["yes","no"],1),###
         Constraint("preprocess",CAN_BE_ANY_OF,
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],0),###
         Constraint("preprocess",CAN_BE_ANY_OF,
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],1),###
         Consequence("contents", SET_TO, "class0,class1"),
         Consequence("format", SAME_AS_CONSTRAINT, 0),
         Consequence("logged", SAME_AS_CONSTRAINT, 0),
         Consequence("missing_values",SAME_AS_CONSTRAINT,1),
         Consequence("combat_norm",SAME_AS_CONSTRAINT,1),
         Consequence("quantile_norm",SAME_AS_CONSTRAINT,1),
         Consequence("dwd_norm",SAME_AS_CONSTRAINT,1),
         Consequence("bfrm_norm",SAME_AS_CONSTRAINT,1),
         Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,1),
         Consequence("predataset",SAME_AS_CONSTRAINT,1),###
         Consequence("preprocess",SAME_AS_CONSTRAINT,1),###
         Constraint("processing_step",MUST_BE,"merge",1),
         Consequence("processing_step",SET_TO_ONE_OF,["merge",
                                    "normalize",
                                    "processed"])
        ),
        Module(
        "normalize_samples_with_quantile",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint('missing_values',MUST_BE,"no"),
        Constraint("quantile_norm", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence('missing_values',SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SET_TO, "yes"),
        Constraint("processing_step",MUST_BE,"merge"),
        Consequence("processing_step",SET_TO_ONE_OF,["merge",
                                    "normalize",
                                    "processed"])
        ),
    Module(
        "normalize_samples_with_bfrm",
        SignalFile,SignalFile,
        UserInputDef("num_factors",1),
        Constraint('format',MUST_BE,"tdf"),
        Constraint('logged',MUST_BE,"yes"),
        Constraint('missing_values',MUST_BE,"no"),
        Constraint('bfrm_norm',MUST_BE,"no"),
        Consequence('format',SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('missing_values',SAME_AS_CONSTRAINT),
        Consequence('bfrm_norm',SET_TO,'yes'),
        Constraint("processing_step",MUST_BE,"merge"),
        Consequence("processing_step",SET_TO_ONE_OF,["merge",
                                    "normalize",
                                    "processed"])),

    Module(###
        "normalize_samples_with_combat",  
        [ClassLabelFile,SignalFile],SignalFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("format",MUST_BE,"tdf",1),
        Constraint("logged",MUST_BE,"yes",1),
        Constraint("missing_values",MUST_BE,"no",1),
        Constraint("combat_norm",MUST_BE,"no",1),
        Constraint("predataset",CAN_BE_ANY_OF,["yes","no"],1),###
        Constraint("preprocess",CAN_BE_ANY_OF,
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],1),###
        Consequence("predataset",SAME_AS_CONSTRAINT,1),###
        Consequence("format",SAME_AS_CONSTRAINT,1),
        Consequence("logged",SAME_AS_CONSTRAINT,1),
        Consequence("missing_values",SAME_AS_CONSTRAINT,1),
        Consequence("preprocess",SAME_AS_CONSTRAINT,1),  ###
        Consequence("combat_norm",SET_TO,"yes"),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],0),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        Constraint("processing_step",MUST_BE,"merge",1),
        Consequence("processing_step",SET_TO_ONE_OF,["merge",
                                    "normalize",
                                    "processed"])),
           
    Module(###
        "normalize_samples_with_dwd",  
        [SignalFile,ClassLabelFile],SignalFile,
        Constraint("cls_format",MUST_BE,'cls',1),
        Constraint("format",MUST_BE,"tdf",0),
        Constraint("logged",MUST_BE,"yes",0),
        Constraint("missing_values",MUST_BE,"no",0),
        Constraint("dwd_norm",MUST_BE,"no",0),
        
        Consequence("format",SAME_AS_CONSTRAINT,0),
        Consequence("logged",SAME_AS_CONSTRAINT,0),
        Consequence("missing_values",SAME_AS_CONSTRAINT,0),
        Consequence("dwd_norm",SET_TO,"yes"),
        Constraint("preprocess",CAN_BE_ANY_OF,
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],0),###
        Consequence("preprocess",SAME_AS_CONSTRAINT,0),  ###
        Constraint("predataset",CAN_BE_ANY_OF,["yes","no"],0),###
        Consequence("predataset",SAME_AS_CONSTRAINT,0),###
        Constraint("processing_step",MUST_BE,"merge",0),
        Consequence("processing_step",SET_TO_ONE_OF,["merge",
                                    "normalize",
                                    "processed"])),

    Module(###
        "normalize_samples_with_shiftscale",  
        [ClassLabelFile,SignalFile],SignalFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("format",MUST_BE,"tdf",1),
        Constraint("logged",MUST_BE,"yes",1),
        Constraint("missing_values",MUST_BE,"no",1),
        Constraint("shiftscale_norm",MUST_BE,"no",1),
        Constraint("predataset",CAN_BE_ANY_OF,["yes","no"],1),##
        Consequence("predataset",SAME_AS_CONSTRAINT,1),###
        Constraint("preprocess",CAN_BE_ANY_OF,
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],1),###
        Consequence("preprocess",SAME_AS_CONSTRAINT,1),  ###
        Consequence("format",SAME_AS_CONSTRAINT,1),
        Consequence("logged",SAME_AS_CONSTRAINT,1),
        Consequence("missing_values",SAME_AS_CONSTRAINT,1),
        Consequence("shiftscale_norm",SET_TO,"yes"),
        Constraint("processing_step",MUST_BE,"merge",1),
        Consequence("processing_step",SET_TO_ONE_OF,["merge",
                                    "normalize",
                                    "processed"])),
    ###normalize
    Module(
        "check_gene_center",
        SignalFile,SignalFile,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_center",MUST_BE,"unknown"),
        # Why does gene_normalize matter here?
        Constraint("gene_normalize", MUST_BE, "unknown"),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("missing_values",SAME_AS_CONSTRAINT),
        Consequence("gene_center",BASED_ON_DATA,["no", "mean", "median"]),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Constraint("processing_step",MUST_BE,"normalize"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "normalize",
                                    "processed"])),
        
    Module(
        "check_gene_normalize",
        SignalFile,SignalFile,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,"unknown"),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("missing_values",SAME_AS_CONSTRAINT),
        Consequence(
            "gene_normalize", BASED_ON_DATA,
            ["no", "variance", "sum_of_squares"]),
        Constraint("processing_step",MUST_BE,"normalize"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "normalize",
                                    "processed"])),
        
    Module(   
        "convert_signal_to_pcl",
        SignalFile,SignalFile,
        Constraint("format",MUST_BE,'tdf'),
        Consequence("format",SET_TO,'pcl'),
        Constraint("processing_step",MUST_BE,"normalize"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "normalize",#"order","annotate","filter",
                                    "processed"])),
    Module(
        "center_genes",
        SignalFile,SignalFile,
        Constraint("format",MUST_BE,"pcl"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_center",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,'unknown'),
        Consequence("format",SET_TO,"tdf"),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("missing_values",SAME_AS_CONSTRAINT),
        Consequence("gene_center",SET_TO_ONE_OF,["mean", "median"]),
        Consequence("gene_normalize",SAME_AS_CONSTRAINT),
        Constraint("processing_step",MUST_BE,"normalize"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "normalize",#"order","annotate","filter",
                                    "processed"])),
    Module(
        "normalize_genes",
        SignalFile,SignalFile,
        Constraint("format",MUST_BE,"pcl"),
        Constraint("logged",MUST_BE,"yes"),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,"no"),
        Consequence("format",SET_TO,"tdf"),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("missing_values",SAME_AS_CONSTRAINT),
        Consequence("gene_normalize",SET_TO_ONE_OF,["variance", "sum_of_squares"]),
        Constraint("processing_step",MUST_BE,"normalize"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "normalize",#"order","annotate","filter",
                                    "processed"])),
    ###transfer
    Module(
        "transfer",
        SignalFile,PrettySignalFile,
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
            "gene_center", CAN_BE_ANY_OF, [ "no", "mean", "median"]),
        Constraint(
            "gene_normalize", CAN_BE_ANY_OF,
            [ "no", "variance", "sum_of_squares"]),
        Constraint(
            "missing_algorithm", CAN_BE_ANY_OF,
            ["none", "median_fill", "zero_fill"]),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
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
        Consequence("filter", SAME_AS_CONSTRAINT, 0),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("combat_norm", SAME_AS_CONSTRAINT, 0),
        Consequence("gene_center", SAME_AS_CONSTRAINT, 0),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT, 0),
        Constraint("processing_step",MUST_BE,"normalize"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                   "order","annotate","filter",
                                    "processed"])
        ),
    ## Order
    Module(  ###
        "rank_genes_by_class_neighbors",
        [ClassLabelFile,PrettySignalFile],GeneListFile,
        UserInputDef("cn_num_neighbors",50),
        UserInputDef("cn_num_perm",100),
        UserInputDef("cn_user_pval",0.5),
        UserInputDef("cn_min_threshold",10),
        UserInputDef("cn_max_threshold",16000),
        UserInputDef("cn_min_folddiff",5),
        UserInputDef("cn_abs_diff",50),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("gene_order", MUST_BE,"no",1),
        Consequence("gene_order",SET_TO,"class_neighbors"),
        Consequence("cn_mean_or_median",SET_TO_ONE_OF, ['mean', 'median']),
        Consequence("cn_ttest_or_snr",SET_TO_ONE_OF, ['t_test','snr']),
        Consequence("cn_filter_data",SET_TO_ONE_OF, ['yes','no']),
        Constraint("processing_step",MUST_BE,"order",1)),
    
    Module(###
         "rank_genes_by_sample_ttest",
         [ClassLabelFile,PrettySignalFile],GeneListFile,
        UserInputDef("gene_select_threshold",0.05),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("gene_order", MUST_BE,"no",1),
        Consequence("gene_order",SET_TO_ONE_OF, ["t_test_p", "t_test_fdr"]),
        Constraint("processing_step",MUST_BE,"order",1)),
         
    Module(  ###
         "reorder_genes",  
         [GeneListFile,PrettySignalFile], PrettySignalFile,
         Constraint("gene_order", CAN_BE_ANY_OF, ['t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list"],0),
         Constraint("gene_order", MUST_BE,"no",1),
         Constraint(
            "preprocess", CAN_BE_ANY_OF,
            ["unknown", "illumina", "agilent", "mas5", "rma", "loess"], 1),
         Constraint("gene_center",CAN_BE_ANY_OF,['median','mean','no'],1),
         Constraint("gene_normalize",CAN_BE_ANY_OF,['sum_of_squares','variance','no'],1),
         Consequence("gene_center",SAME_AS_CONSTRAINT,1),
         Consequence("gene_normalize",SAME_AS_CONSTRAINT,1),
         Consequence("gene_order",SET_TO_ONE_OF, ['t_test_p', "t_test_fdr",
                         'class_neighbors', "gene_list"]),
         Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
         Constraint("processing_step",MUST_BE,"order",1),
         Consequence("processing_step",SET_TO_ONE_OF,[
                                    "order","annotate","filter",
                                    "processed"])
        ),
    ##Annotate
    Module(
         'annotate_probes',
         PrettySignalFile,PrettySignalFile,
         Constraint("annotate", MUST_BE,"no"),
         Constraint("platform", MUST_BE,"no"),
         Consequence("annotate",SET_TO,"yes"),
         Consequence("platform",SAME_AS_CONSTRAINT),
         Constraint("processing_step",MUST_BE,"annotate"),
         Consequence("processing_step",SET_TO_ONE_OF,[
                                    "annotate","filter",
                                    "processed"])),
    Module(### 
       "relabel_samples",  
        [RenameFile,PrettySignalFile],PrettySignalFile,
         Constraint("rename_sample",MUST_BE,"no",1),
         Constraint("gene_center",CAN_BE_ANY_OF,['median','mean','no'],1),
         Constraint("gene_normalize",CAN_BE_ANY_OF,['sum_of_squares','variance','no'],1),
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],0),
         Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],1),
         Consequence("rename_sample",SET_TO,"yes"),
         Consequence("contents", SAME_AS_CONSTRAINT,0),
         Consequence("gene_center",SAME_AS_CONSTRAINT,1),
         Consequence("gene_normalize",SAME_AS_CONSTRAINT,1),
         Constraint("processing_step",MUST_BE,"annotate",1),
         Consequence("processing_step",SET_TO_ONE_OF,[
                                    "annotate","filter",
                                    "processed"])
       ),
    Module(
         'add_crossplatform_probeid',
         PrettySignalFile,PrettySignalFile,
         UserInputDef("platform_name"),
         Constraint("platform", MUST_BE,"no"),
         Consequence("platform",SET_TO,"yes"),
         Constraint("processing_step",MUST_BE,"annotate"),
         Consequence("processing_step",SET_TO_ONE_OF,[
                                    "annotate","filter",
                                    "processed"])),
    #Filter
    Module(
        'remove_duplicate_genes',
        PrettySignalFile, PrettySignalFile,
        Constraint("annotate", MUST_BE,"yes"),
        Constraint("num_features", MUST_BE,"no"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Constraint("unique_genes", MUST_BE,'no'),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
        Consequence(
            "unique_genes", SET_TO_ONE_OF,
            ['average_genes', 'high_var', 'first_gene']),
        Constraint("processing_step",MUST_BE,"filter"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "filter",
                                    "processed"])),
      
    Module(
         'select_first_n_genes',
        PrettySignalFile,PrettySignalFile,
        UserInputDef("num_features_value",500),
        Constraint("num_features", MUST_BE,"no"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("num_features",SET_TO,"yes"),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT),
         Constraint("processing_step",MUST_BE,"filter"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "filter",
                                    "processed"])), 
      
     Module(
        'remove_duplicate_probes',
        PrettySignalFile,PrettySignalFile,
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("duplicate_probe",SET_TO,'high_var_probe'),
        Constraint("processing_step",MUST_BE,"filter"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "filter",
                                    "processed"])),
    Module(
         'select_probe_by_best_match',
        PrettySignalFile,PrettySignalFile,
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("duplicate_probe",SET_TO,'closest_probe'),
         Constraint("processing_step",MUST_BE,"filter"),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "filter",
                                    "processed"])),
    
    Module(  ###
        "filter_genes_by_fold_change_across_classes",
        [ClassLabelFile,PrettySignalFile],PrettySignalFile,
        UserInputDef("group_fc_num"),
        Constraint("preprocess",
            CAN_BE_ANY_OF,
            ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],1),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("group_fc", MUST_BE,"no",1),
        Constraint("num_features", MUST_BE,"no",1),
        Constraint("duplicate_probe", MUST_BE,"no",1),
        Constraint("unique_genes", MUST_BE,"no",1),
        Constraint("gene_center",CAN_BE_ANY_OF,['median','mean','no'],1),
        Constraint("gene_normalize",CAN_BE_ANY_OF,['sum_of_squares','variance','no'],1),
        Consequence("gene_center",SAME_AS_CONSTRAINT,1),
        Consequence("gene_normalize",SAME_AS_CONSTRAINT,1),
        Consequence("num_features",SAME_AS_CONSTRAINT,1),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT,1),
        Consequence("unique_genes",SAME_AS_CONSTRAINT,1),
        Consequence("group_fc",SET_TO,"yes"),
        Consequence("preprocess",SAME_AS_CONSTRAINT,1),
        Constraint("processing_step",MUST_BE,"filter",1),
        Consequence("processing_step",SET_TO_ONE_OF,[
                                    "filter",
                                    "processed"])),
    Module(###
        "convert_label_to_cls",   
        [ClassLabelFile,SignalFile],ClassLabelFile,
        Constraint("cls_format",MUST_BE,'label',0),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("logged",MUST_BE,"yes",1),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],0),
        Constraint("contents",CAN_BE_ANY_OF,[
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        Consequence("contents", SAME_AS_CONSTRAINT,1),
        Consequence("cls_format",SET_TO,'cls'),
        )

    ]

list_files=[RenameFile,AgilentFiles,CELFiles,ControlFile,ExpressionFiles,
            GPRFiles,GEOSeries,IDATFiles,ClassLabelFile,ILLUFolder,GeneListFile,
            SignalFile,PrettySignalFile]

