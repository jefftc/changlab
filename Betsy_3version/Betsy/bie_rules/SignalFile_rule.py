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


GEOSeries = DataType("GEOSeries")
ExpressionFiles = DataType("ExpressionFiles")
CELFiles = DataType(
    "CELFiles",
    Attribute("version", ["unknown", "cc", "v3", "v4"], "unknown", "v3"),
    )
RenameFile = DataType(
    'RenameFile')
AgilentFiles = DataType(
    "AgilentFiles")

ControlFile = DataType(
    "ControlFile",
    Attribute(
        'preprocess',["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        "unknown","unknown"),
    Attribute(
        'missing_values',["unknown", "no", "yes"],
        "no","no"),
    Attribute(
        "missing_algorithm",["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill"),
    Attribute(
        "logged",["unknown", "no", "yes"],"no","no"),
    Attribute(
        'format',["unknown", "tdf", "gct", "jeffs", "pcl", "res", "xls"],
        "gct","gct"),
    )
GPRFiles = DataType(
    "GPRFiles")

IDATFiles = DataType(
    "IDATFiles")

ClassLabelFile = DataType(
    "ClassLabelFile",
    Attribute(
    "contents",["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "unspecified"],'unspecified','unspecified'),
    Attribute("cls_format",['cls','label','unknown'],"unknown","cls")
    )

ILLUFolder = DataType(
    "ILLUFolder", 
    Attribute(
        "illu_manifest",ILLU_MANIFEST,
        'HumanHT-12_V4_0_R2_15002873_B.txt','HumanHT-12_V4_0_R2_15002873_B.txt'),
    Attribute(
        'illu_chip',ILLU_CHIP,
        'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip','ilmn_HumanHT_12_V4_0_R1_15002873_B.chip'),
    Attribute('illu_bg_mode',['false', 'true'], "false", "false"),
    Attribute('illu_coll_mode',['none', 'max', 'median'], "none","none"),
    #bie3.Attribute('illu_clm',bie3.ANYATOM, "",""),
    #bie3.Attribute('illu_custom_chip',bie3.ANYATOM, "",""),
    #bie3.Attribute('illu_custom_manifest',bie3.ANYATOM, "",""),
    )


GeneListFile=DataType(
    "GeneListFile",
    #bie3.Attribute('cn_num_neighbors',bie3.ANYATOM, "",""),
    #bie3.Attribute('cn_num_perm',bie3.ANYATOM, "",""),
    #bie3.Attribute('cn_user_pval',bie3.ANYATOM, "",""),  
    Attribute('cn_mean_or_median',['mean', 'median'], 'mean','mean'),
    Attribute('cn_ttest_or_snr',['t_test','snr'], 't_test','t_test'),
    Attribute('cn_filter_data',['yes','no'], 'no','no'),
    #bie3.Attribute('cn_min_threshold',bie3.ANYATOM, "",""),
    #bie3.Attribute('cn_max_threshold',bie3.ANYATOM, "",""),
    #bie3.Attribute('cn_min_folddiff',bie3.ANYATOM, "",""),
    #bie3.Attribute('cn_abs_diff',bie3.ANYATOM, "",""),
    #bie3.Attribute('gene_select_threshold',bie3.ANYATOM,"",""),
    Attribute('gene_order',['no', "gene_list", "class_neighbors",
                          "t_test_p", "t_test_fdr"], 'gene_list',"gene_list"))
   
SignalFile = DataType(
    "SignalFile",
    Attribute("format", ["unknown", "tdf", "pcl", "gct", "res", "jeffs"],
              "unknown", "tdf"),

    # Properties of the data.
    Attribute("preprocess",
              ["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
              "unknown", "unknown"),
    Attribute("missing_values", ["unknown", "no", "yes"], "unknown", "no"),
    Attribute("missing_algorithm", ["none", "median_fill", "zero_fill"],
    ##          "none", "none"),
              "zero_fill","zero_fill"),
    Attribute("logged", ["unknown", "no", "yes"], "unknown", "yes"),
    Attribute("filter", ["no", "yes"], "no", "no"),
    # Normalization of the data.
    Attribute("dwd_norm", ["no", "yes"], "no", "no"),
    Attribute("bfrm_norm", ["no", "yes"], "no", "no"),
    Attribute("quantile_norm", ["no", "yes"], "no", "no"),
    Attribute("shiftscale_norm", ["no", "yes"], "no", "no"),
    Attribute("combat_norm", ["no", "yes"], "no", "no"),
    
    # Other attributes.
    Attribute("predataset", ["no", "yes"], "no", "no"),
    Attribute("rename_sample", ["no", "yes"], "no", "no"),
    Attribute("contents", [
        "unspecified", "train0", "train1", "test", 'class0,class1,test',
        "class0", "class1", "class0,class1"],
              "unspecified", "unspecified"),
    )
    

all_modules = [
    Module(
        "download_geo", GEOSeries, ExpressionFiles,
         UserInput("GSEID"), UserInput("GPLID","")),
    Module(
        "extract_CEL_files", ExpressionFiles, CELFiles,
        Consequence("version", SET_TO, "unknown"),
        ),
    Module(
        "detect_CEL_version",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "unknown"),
        Consequence("version", BASED_ON_DATA, ["cc", "v3", "v4"]),
        ),
    Module(
        "convert_CEL_cc_to_CEL_v3",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "cc"),
        Consequence("version", SET_TO, "v3"),
        ),
     # IDATFiles
    Module("extract_illumina_idat_files",
                ExpressionFiles, IDATFiles),
    Module(
        "preprocess_illumina",
        IDATFiles, ILLUFolder,
        Consequence('illu_manifest',SET_TO_ONE_OF,ILLU_MANIFEST),
        Consequence('illu_chip',SET_TO_ONE_OF,ILLU_CHIP),
        Consequence('illu_bg_mode',SET_TO_ONE_OF,["false", "true"]),
        Consequence('illu_coll_mode',SET_TO_ONE_OF,["none", "max", "median"]),
        UserInput("illu_clm",''),
        UserInput("illu_custom_chip",''),
        UserInput("illu_custom_manifest",'')
        ),
       Module(
        "get_illumina_signal",
         ILLUFolder, SignalFile,
         Consequence('preprocess',SET_TO,"illumina"),
         Consequence('format',SET_TO,"gct"),
         Consequence('logged',SET_TO,"no"),
         Consequence('missing_values',SET_TO,"unknown")),
    Module(
        "get_illumina_control",
         ILLUFolder,ControlFile,
         Consequence('preprocess',SET_TO,"illumina"),
         Consequence("format",SET_TO,"gct"),
         Consequence("logged",SET_TO,"no")
        ),
    
    # AgilentFiles
    Module(
        "extract_agilent_files", ExpressionFiles, AgilentFiles),
    Module(
        "preprocess_agilent",
         AgilentFiles,SignalFile,
        Consequence('logged',SET_TO,"unknown"),
        Consequence('preprocess',SET_TO,"agilent"),
        Consequence('format',SET_TO,"tdf")),

    # GPRFiles
    Module(
        "extract_gpr_files", ExpressionFiles, GPRFiles),
    Module(
        "normalize_with_loess",
        GPRFiles,SignalFile,
        Consequence("format",SET_TO,"tdf"),
        Consequence("logged",SET_TO,"unknown"),
        Consequence("preprocess",SET_TO,"loess")),
    
    Module(
        "preprocess_rma",
        CELFiles, SignalFile,
        Constraint("version", CAN_BE_ANY_OF, ["v3", "v4"]),
        Consequence("logged", SET_TO, "yes"),
        Consequence("preprocess", SET_TO, "rma"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence("quantile_norm", SET_TO, "yes"),
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
        ),
 
    Module(
        "convert_signal_to_tdf",
        SignalFile, SignalFile,
        Constraint("format", CAN_BE_ANY_OF, ['pcl', 'gct', 'res', 'jeffs']),
        Consequence("format", SET_TO, "tdf"),
        ),
    Module(
        "check_for_log",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", BASED_ON_DATA, ["yes", "no"]),
        ),
    Module(
        "log_signal",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "yes"),
        ),
    Module(
        "check_for_missing_values",
        SignalFile, SignalFile,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", BASED_ON_DATA, ["no", "yes"]),
        ),
    Module(
        "filter_genes_by_missing_values",
        SignalFile, SignalFile,
        UserInput("filter_genes_with_missing_values", 0.50),
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "yes"),
        Constraint("missing_values", MUST_BE, "yes"),
        Constraint("filter", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence("missing_values", SAME_AS_CONSTRAINT),
        Consequence("filter", SET_TO, "yes")
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
        ),
    Module(
        "fill_missing_with_median",
        SignalFile,SignalFile,
        Constraint('format',MUST_BE,'tdf'),
        Constraint('logged',MUST_BE,'yes'),
        Constraint('missing_values',MUST_BE,'yes'),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('missing_algorithm',SET_TO,"median_fill",side_effect=True),
        Consequence('missing_values',SET_TO,'no')),
   
   Module(
        "filter_and_threshold_genes",
        SignalFile,SignalFile,
        Constraint('format',MUST_BE,"tdf"),
        Constraint('logged',MUST_BE,"no"),
        Constraint('predataset',MUST_BE,"no"),
        Consequence('format',SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('predataset',SET_TO,'yes')),
    Module(
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
         Consequence("contents", SET_TO, "class0,class1"),
         Consequence("format", SAME_AS_CONSTRAINT, 0),
         Consequence("logged", SAME_AS_CONSTRAINT, 0),
         Consequence("missing_values",SAME_AS_CONSTRAINT,1),
         Consequence("combat_norm",SAME_AS_CONSTRAINT,1),
         Consequence("quantile_norm",SAME_AS_CONSTRAINT,1),
         Consequence("dwd_norm",SAME_AS_CONSTRAINT,1),
         Consequence("bfrm_norm",SAME_AS_CONSTRAINT,1),
         Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,1)
        ),
    Module(
       "relabel_samples",
        [RenameFile,SignalFile],SignalFile,
         Constraint('format',MUST_BE,'tdf',1),
         Constraint("rename_sample",MUST_BE,"no",1),
         Constraint("logged",MUST_BE,"yes",1),
         Constraint("missing_values",MUST_BE,'no',1),
         Constraint("combat_norm",MUST_BE,'no',1),
         Constraint("quantile_norm",MUST_BE,"no",1),
         Constraint("dwd_norm",MUST_BE,"no",1),
         Constraint("bfrm_norm",MUST_BE,"no",1),
         Constraint("shiftscale_norm",MUST_BE,"no",1),
         Consequence("format",SAME_AS_CONSTRAINT,1),
         Consequence("rename_sample",SET_TO,"yes"),
         Consequence("logged",SAME_AS_CONSTRAINT,1),
         Consequence("missing_values",SAME_AS_CONSTRAINT,1),
         Consequence("combat_norm",SAME_AS_CONSTRAINT,1),
         Consequence("quantile_norm",SAME_AS_CONSTRAINT,1),
         Consequence("dwd_norm",SAME_AS_CONSTRAINT,1),
         Consequence("bfrm_norm",SAME_AS_CONSTRAINT,1),
         Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,1)),

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
        ),
    Module(
        "normalize_samples_with_bfrm",
        SignalFile,SignalFile,
        UserInput("num_factors",1),
        Constraint('format',MUST_BE,"tdf"),
        Constraint('logged',MUST_BE,"yes"),
        Constraint('missing_values',MUST_BE,"no"),
        Constraint('bfrm_norm',MUST_BE,"no"),
        Consequence('format',SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('missing_values',SAME_AS_CONSTRAINT),
        Consequence('bfrm_norm',SET_TO,'yes')),

    Module(
        "normalize_samples_with_combat",
        [ClassLabelFile,SignalFile],SignalFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("format",MUST_BE,"tdf",1),
        Constraint("logged",MUST_BE,"yes",1),
        Constraint("missing_values",MUST_BE,"no",1),
        Constraint("combat_norm",MUST_BE,"no",1),
        Consequence("format",SAME_AS_CONSTRAINT,1),
        Consequence("logged",SAME_AS_CONSTRAINT,1),
        Consequence("missing_values",SAME_AS_CONSTRAINT,1),
        Consequence("combat_norm",SET_TO,"yes")),
           
    Module(
        "normalize_samples_with_dwd",
        [ClassLabelFile,SignalFile],SignalFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("format",MUST_BE,"tdf",1),
        Constraint("logged",MUST_BE,"yes",1),
        Constraint("missing_values",MUST_BE,"no",1),
        Constraint("dwd_norm",MUST_BE,"no",1),
        Consequence("format",SAME_AS_CONSTRAINT,1),
        Consequence("logged",SAME_AS_CONSTRAINT,1),
        Consequence("missing_values",SAME_AS_CONSTRAINT,1),
        Consequence("dwd_norm",SET_TO,"yes")),

    Module(
        "normalize_samples_with_shiftscale",
        [ClassLabelFile,SignalFile],SignalFile,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("format",MUST_BE,"tdf",1),
        Constraint("logged",MUST_BE,"yes",1),
        Constraint("missing_values",MUST_BE,"no",1),
        Constraint("shiftscale_norm",MUST_BE,"no",1),
        Consequence("format",SAME_AS_CONSTRAINT,1),
        Consequence("logged",SAME_AS_CONSTRAINT,1),
        Consequence("missing_values",SAME_AS_CONSTRAINT,1),
        Consequence("shiftscale_norm",SET_TO,"yes")),
    Module(
        "convert_label_to_cls",
        [ClassLabelFile,SignalFile],ClassLabelFile,
        Constraint("cls_format",MUST_BE,'label',0),
        Constraint("format",MUST_BE,'tdf',1),
        Constraint("logged",MUST_BE,"yes",1),
        Consequence("cls_format",SET_TO,'cls'))

    ]

list_files=[RenameFile,AgilentFiles,CELFiles,ControlFile,ExpressionFiles,
            GPRFiles,GEOSeries,IDATFiles,ClassLabelFile,ILLUFolder,GeneListFile,
            SignalFile]

