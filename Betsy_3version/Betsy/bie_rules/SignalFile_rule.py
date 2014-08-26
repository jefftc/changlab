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
CONTENTS = ["train0", "train1","test", "class0,class1,test",
             "class0", "class1", "class0,class1","unspecified",
            "diff_class0","diff_class1","diff_class0,diff_class1","diff_unspecified"]

PREPROCESS = ["unknown", "illumina", "agilent", "mas5", "rma", "loess",
              "tcga","rsem"]

GEOSeries = DataType("GEOSeries",
                     AttributeDef("contents",CONTENTS,
                                  'unspecified','unspecified',help="contents"),
                     help="GEOID to download from the GEO database")
GEOfamily = DataType("GEOfamily",
                     AttributeDef("contents",CONTENTS,
                                  'unspecified','unspecified',help="contents"),
                     help="GEO fmaily soft file download from the GEO database")
TCGAID = DataType("TCGAID",
                     AttributeDef("contents",CONTENTS,
                                  'unspecified','unspecified',help="contents"),
                     help="TCGA ID to download from TCGA database")
TCGAFile = DataType("TCGAFile",
                     AttributeDef("contents",CONTENTS,
                                  'unspecified','unspecified',help="contents"),
                     AttributeDef("data",['RSEM_genes','RSEM_exons',
                                          'humanmethylation450','mirnaseq',
                                          'rppa','clinical'],'RSEM_genes',
                                          'RSEM_genes',help="TCGA data type"),
                     help="TCGA file download from TCGA database")


ExpressionFiles = DataType("ExpressionFiles",
                           AttributeDef("contents",
                                        CONTENTS,
                                        'unspecified','unspecified',help="contents"),
                           help="Expression file folder, can be CELFiles, IDATFiles,"\
                                 "AgilentFile,GPRFiles")

CELFiles = DataType(
    "CELFiles",
    AttributeDef("version", ["unknown", "cc", "v3_v4"], "unknown", "v3_v4",help="cel file version"),
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="A folder of cel files")

RenameFile = DataType(
    'RenameFile',
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    AttributeDef("labels_from",["title","description"],'title','title',help="labels from title or description"),
    help="A file used to rename the sample name in the gene expression file.")

AgilentFiles = DataType(
    "AgilentFiles",
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="A folder of agilent files.")

ControlFile = DataType(
    "ControlFile",
    AttributeDef(
        'preprocess',["illumina"],
        "illumina","illumina",help="preprocess for ControlFile"),
    AttributeDef(
        'missing_values',["unknown", "no", "yes"],
        "no","no",help="missing values yes or not"),
    AttributeDef(
        "missing_algorithm",["none", "median_fill", "zero_fill"],
        "zero_fill", "zero_fill",help="missing algorithm"),
    AttributeDef(
        "logged",["no"],"no","no",help="logged yes or not"),
    AttributeDef(
        'format',["gct"],
        "gct","gct",help="file format"),
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="The control file in the ILLUFolder")
    
    
GPRFiles = DataType(
    "GPRFiles",
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="A folder of GPR files.")

IDATFiles = DataType(
    "IDATFiles",
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="A folder of IDAFiles.")

ClassLabelFile = DataType(
    "ClassLabelFile",
    AttributeDef(
    "contents",CONTENTS,'unspecified','unspecified',help="contents"),
    AttributeDef("cls_format",['cls','label','unknown'],"unknown","cls",help="cls format for ClassLabelFile"),
    help="The Class label file, can be cls format or label format")

ILLUFolder = DataType(
    "ILLUFolder", 
    AttributeDef(
        "illu_manifest",ILLU_MANIFEST,
        'HumanHT-12_V4_0_R2_15002873_B.txt','HumanHT-12_V4_0_R2_15002873_B.txt',help="illumina manifest"),
    AttributeDef(
        'illu_chip',ILLU_CHIP,
        'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip','ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',help="illumina chip type"),
    AttributeDef('illu_bg_mode',['false', 'true'], "false", "false",help="illumina background mode"),
    AttributeDef('illu_coll_mode',['none', 'max', 'median'], "none","none",help="illumina coll mode"),
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="A folder generated from preprocess_illumina," \
          "it contains SignalFile_Postprocess and ControlFile.")


GeneListFile=DataType(
    "GeneListFile", 
    AttributeDef('cn_mean_or_median',['mean', 'median'], 'mean','mean',help="class neighbors mean or median"),
    AttributeDef('cn_ttest_or_snr',['t_test','snr'], 't_test','t_test',help="class neighbors ttest or snr"),
    AttributeDef('cn_filter_data',['yes','no'], 'no','no',help="class neighbors filter data or not"),
    AttributeDef('gene_order',['no', "gene_list", "class_neighbors",
                               "t_test_p", "t_test_fdr",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change'],
                 't_test_p',"t_test_p",help="gene order method"),
    AttributeDef("contents",CONTENTS,'unspecified','unspecified',help="contents"),
    help="A file contains a list of genes.")
   
SignalFile_Postprocess = DataType(
    "SignalFile_Postprocess",
    AttributeDef("format", ["unknown", "tdf", "pcl", "gct", "res", "jeffs"],
              "unknown", "tdf",help="file format"),
    # Properties of the data.
    AttributeDef("preprocess", PREPROCESS, "unknown", "unknown",help="preprocess method"),
    AttributeDef("logged", ["unknown", "no", "yes"], "unknown", "yes",help="logged or not"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    help="The input SignalFile which care the format,preprocess," \
          "logged,predataset,contents.")

SignalFile_Impute = DataType(
    "SignalFile_Impute",
    # Properties of the data.
    AttributeDef("preprocess", PREPROCESS, "unknown", "unknown",help="preprocess method"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef("missing_values", ["unknown", "no", "yes"], "unknown", "no",
                 help="missing values unknown,yes or not"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    help="The SignalFile after SignalFile_Postprocess, care missing_values,missing_algorithm"\
          "and filter.")


SignalFile_Merge = DataType(
    "SignalFile_Merge",
    # Properties of the data.
    AttributeDef("preprocess",
              PREPROCESS,
              "unknown", "unknown",help="preprocess method"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no",help="dwd normalization"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no",help="bfrm normalization"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no",help="quantile normalization"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no",help="shiftscale normalization"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no",help="combat normalization"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    help="The SignalFile after SiganlFile_Impute, care dwd_norm,bfrm_norm,"\
          "quantile_norm,shiftscale_norm,combat_norm.")

SignalFile_Normalize = DataType(
    "SignalFile_Normalize",
    # Properties of the data.
    AttributeDef("preprocess", PREPROCESS, "unknown", "unknown",help="preprocess method"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("format", ["tdf", "pcl"], "tdf", "tdf",help="file format"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no",help="dwd normalization"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no",help="bfrm normalization"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no",help="quantile normalization"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no",help="shiftscale normalization"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no",help="combat normalization"),
    AttributeDef(
        "gene_center", ["unknown", "no", "mean", "median"],
        "unknown", "no",help="gene center method"),
    AttributeDef(
        "gene_normalize", ["unknown", "no", "variance", "sum_of_squares"],
        "unknown", "no",help="gene normalize method"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    help="The SignalFile after SiganlFile_Merge, care gene_center,"\
          "gene_normalize.")    
    
    
    

SignalFile_Order = DataType(
    "SignalFile_Order",
    # Properties of the data.
    AttributeDef("preprocess", PREPROCESS, "unknown", "unknown",help="preprocess method"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no",help="dwd normalization"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no",help="bfrm normalization"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no",help="quantile normalization"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no",help="shiftscale normalization"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no",help="combat normalization"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef(
        "gene_center", [ "no", "mean", "median"],
        "no", "no",help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"],
        "no", "no",help="gene normalize method"),
    AttributeDef(
        "gene_order",
        ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr",
         'diff_ttest','diff_sam','diff_ebayes','diff_fold_change'],
       "no", "no",help="gene order method"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    help="The SignalFile after SiganlFile_Normalize, care gene_order.")

SignalFile_Annotate= DataType( 
    "SignalFile_Annotate",
    # Properties of the data.
    AttributeDef("preprocess", PREPROCESS, "unknown", "unknown",help="preprocess method"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no",help="dwd normalization"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no",help="bfrm normalization"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no",help="quantile normalization"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no",help="shiftscale normalization"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no",help="combat normalization"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef(
        "gene_center", [ "no", "mean", "median"],
        "no", "no",help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"],
        "no", "no",help="gene normalize method"),
    AttributeDef(
        "gene_order",
        ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr",
         'diff_ttest','diff_sam','diff_ebayes','diff_fold_change'],
       "no", "no",help="gene order method"),
    AttributeDef("annotate", ["no", "yes"], "no", "no",help="annotate file or not"),
    AttributeDef("rename_sample", ["no", "yes"], "no", "no",help="rename sample or not"),
    AttributeDef("platform", ["yes","no"], "no", "no",help="add platform or not"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    help="The SignalFile after SiganlFile_Order, care annotate,rename_sample,platform.")

SignalFile_Filter= DataType( 
    "SignalFile_Filter",
    # Properties of the data.
    AttributeDef("preprocess",PREPROCESS,"unknown", "unknown",help="preprocess method"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no",help="dwd normalization"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no",help="bfrm normalization"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no",help="quantile normalization"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no",help="shiftscale normalization"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no",help="combat normalization"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef(
        "gene_center", [ "no", "mean", "median"],
        "no", "no",help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"],
        "no", "no",help="gene normalize method"),
    AttributeDef(
        "gene_order",
        ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr",
         'diff_ttest','diff_sam','diff_ebayes','diff_fold_change'],
       "no", "no",help="gene order method"),
    AttributeDef("annotate", ["no", "yes"], "no", "no",help="annotate file or not"),
    AttributeDef("rename_sample", ["no", "yes"], "no", "no",help="rename sample or not"),
    AttributeDef("platform", ["yes","no"], "no", "no",help="add platform or not"),
    AttributeDef("num_features", ["yes","no"], "no", "no",help="select a num of features or not"),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no",help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "closest_probe", "high_var_probe"],
        "no", "no",help="method to remove duplicated probes"),
    AttributeDef("group_fc", ["yes","no"], "no","no",help="group fold change or not"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    AttributeDef("logged", [ "no", "yes"], "yes", "yes",help="logged or not"),
    AttributeDef("format", [ "tdf", "gct"], "tdf", "tdf",help="file format"),
    help="The SignalFile after SiganlFile_Annotate, care num_features,"\
          "unique_genes,duplicate_probe,group_fc,logged,format.")

SignalFile= DataType( 
    "SignalFile",
    # Properties of the data.
    AttributeDef("preprocess",PREPROCESS,"unknown", "unknown",help="preprocess method"),
    AttributeDef("missing_algorithm", ["none", "median_fill", "zero_fill"],
              "zero_fill","zero_fill",help="missing algorithm"),
    AttributeDef("filter", ["no", "yes"], "no", "no",help="filter missing or not"),
    # Normalization of the data.
    AttributeDef("dwd_norm", ["no", "yes"], "no", "no",help="dwd normalization"),
    AttributeDef("bfrm_norm", ["no", "yes"], "no", "no",help="bfrm normalization"),
    AttributeDef("quantile_norm", ["no", "yes"], "no", "no",help="quantile normalization"),
    AttributeDef("shiftscale_norm", ["no", "yes"], "no", "no",help="shiftscale normalization"),
    AttributeDef("combat_norm", ["no", "yes"], "no", "no",help="combat normalization"),
    AttributeDef("predataset", ["no", "yes"], "no", "no",help="predataset or not"),
    AttributeDef(
        "gene_center", [ "no", "mean", "median"],
        "no", "no",help="gene center method"),
    AttributeDef(
        "gene_normalize", [ "no", "variance", "sum_of_squares"],
        "no", "no",help="gene normalize method"),
    AttributeDef(
        "gene_order",
        ["no", "class_neighbors", "gene_list", "t_test_p", "t_test_fdr",
         'diff_ttest','diff_sam','diff_ebayes','diff_fold_change'],
       "no", "no",help="gene order method"),
    AttributeDef("annotate", ["no", "yes"], "no", "no",help="annotate file or not"),
    AttributeDef("rename_sample", ["no", "yes"], "no", "no",help="rename sample or not"),
    AttributeDef("platform", ["yes","no"], "no", "no",help="add platform or not"),
    AttributeDef("num_features", ["yes","no"], "no", "no",help="select a num of features or not"),
    AttributeDef(
        "unique_genes", ["no", "average_genes", "high_var", "first_gene"],
        "no", "no",help="method to get unique genes"),
    AttributeDef(
        "duplicate_probe", ["no", "closest_probe", "high_var_probe"],
        "no", "no",help="method to remove duplicated probes"),
    AttributeDef("group_fc", ["yes","no"], "no","no",help="group fold change or not"),
    AttributeDef("contents", CONTENTS,"unspecified", "unspecified",help="contents"),
    AttributeDef("logged", [ "no", "yes"], "yes", "yes",help="logged or not"),
    AttributeDef("format", [ "tdf", "gct"], "tdf", "tdf",help="file format"),
    help="The SignalFile after SiganlFile_Filter, the attributes are the same as SignalFile_Filter.")

    

all_modules = [
    #TCGA Files
    Module('download_tcga', TCGAID, TCGAFile,
           OptionDef("disease",help="tcga disease type"),
           OptionDef("date","",help="date for tcga disease"),
           Constraint("contents",CAN_BE_ANY_OF,CONTENTS,),
           Consequence("contents",SAME_AS_CONSTRAINT),
           Consequence("data",SET_TO_ONE_OF,['RSEM_genes','RSEM_exons',
                                          'humanmethylation450','mirnaseq',
                                          'rppa','clinical']),
           help="download data from tcga website according to TCGAID"),
    
    Module('preprocess_tcga',TCGAFile,SignalFile_Postprocess,
        Constraint("contents",CAN_BE_ANY_OF, CONTENTS,),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence('logged',SET_TO,"unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence('preprocess',SET_TO,"tcga"),
        Consequence('format',SET_TO,"tdf"),
        help="preprocess tcga file, generate to SignalFile_Postprocess"),
           
    Module(
        "download_geo", GEOSeries, ExpressionFiles,
         OptionDef("GSEID",help="GSEID to download"),
         OptionDef("GPLID","",help="GPDID to download"),
         Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="download GEO data from geo website according to GSEID and GPLID"
        ),
    #CELFiles
    Module(
        "extract_CEL_files", ExpressionFiles, CELFiles,
         Consequence("version", SET_TO, "unknown"),
         Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="extract CEL files from Expression Files folder"
        ),
    Module(
        "detect_CEL_version",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "unknown"),
        Consequence("version", BASED_ON_DATA, ["cc", "v3_v4"]),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="detect the version of cel files"
        ),
    Module(
        "convert_CEL_to_v3",
        CELFiles, CELFiles,
        Constraint("version", MUST_BE, "cc"),
        Consequence("version", SET_TO, "v3_v4"),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="convert cel files to v3 version"
        ),
     # IDATFiles
    Module("extract_illumina_idat_files",
            ExpressionFiles, IDATFiles,
           Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
           help="extract idat files from Expression File folder"),
    Module(
        "preprocess_illumina",
        IDATFiles, ILLUFolder,
        Consequence('illu_manifest',SET_TO_ONE_OF,ILLU_MANIFEST),
        Consequence('illu_chip',SET_TO_ONE_OF,ILLU_CHIP),
        Consequence('illu_bg_mode',SET_TO_ONE_OF,["false", "true"]),
        Consequence('illu_coll_mode',SET_TO_ONE_OF,["none", "max", "median"]),
        OptionDef("illu_clm",'',help="illumina clm"),
        OptionDef("illu_custom_chip",'',help="illumina custom chip name"),
        OptionDef("illu_custom_manifest",'',help='illumina custrom manifest file'),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="preprocess idat files,generate SignalFile_Postprocess"
        ),
       Module(
        "get_illumina_signal",
         ILLUFolder, SignalFile_Postprocess,
         Consequence('preprocess',SET_TO,"illumina"),
         Consequence('format', SET_TO, "gct"),
         Consequence('logged', SET_TO, "no"),
         Consequence('predataset', SET_TO, "no"),
         Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help="extract the SignalFile_Postprocess from ILLUFolder"),
    Module(
        "get_illumina_control",
         ILLUFolder,ControlFile,
         Consequence('preprocess',SET_TO,"illumina"),
         Consequence("format",SET_TO,"gct"),
         Consequence("logged",SET_TO,"no"),
         Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="extract illumina ControlFile from ILLUFolder"
        ),
    
    # AgilentFiles
    Module(
        "extract_agilent_files", ExpressionFiles, AgilentFiles,
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="extract agilent files from ExpressionFiles"),
    Module(
        "preprocess_agilent",
         AgilentFiles,SignalFile_Postprocess,
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence('logged',SET_TO,"unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence('preprocess',SET_TO,"agilent"),
        Consequence('format',SET_TO,"tdf"),
        help="preprocess agilent, generate SignalFile_Postprocess"),

    # GPRFiles
    Module(
        "extract_gpr_files", ExpressionFiles, GPRFiles,
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="extract gpr files from ExpressionFiles"),
    Module(
        "normalize_with_loess",
        GPRFiles,SignalFile_Postprocess,
        Consequence("format",SET_TO,"tdf"),
        Consequence("logged",SET_TO,"unknown"),
        Consequence('predataset', SET_TO, "no"),
        Consequence("preprocess",SET_TO,"loess"),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="normalize GPRFiles,generate SignalFile_Postprocess"),
    
    Module(
        "preprocess_rma",
        CELFiles, SignalFile_Postprocess,
        Constraint("version", MUST_BE,'v3_v4'),
        Consequence("logged", SET_TO, "yes"),
        Consequence("preprocess", SET_TO, "rma"),
        Consequence("format", SET_TO, "jeffs"),
        Consequence('predataset', SET_TO, "no"),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="preprocess CELFiles with rma method,generate SignalFile_Postprocess"
        ),
    Module(
        "preprocess_mas5",
        CELFiles, SignalFile_Postprocess,
        Constraint("version", MUST_BE,'v3_v4'),
        Consequence("logged", SET_TO, "no"),
        Consequence("preprocess", SET_TO, "mas5"),
        Consequence('predataset', SET_TO, "no"),
        Consequence("format", SET_TO, "jeffs"),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="preprocess CELFiles with mas5 method,generate SignalFile_Postprocess"
        ),
    ####postprocess
    Module(
        "convert_signal_to_tdf",
        SignalFile_Postprocess, SignalFile_Postprocess,
        Constraint("format", CAN_BE_ANY_OF, ["unknown", "pcl", "gct", "res", "jeffs"]),
        Consequence("format", SET_TO, "tdf"),
        help="convert SignalFile_Postprocess to tdf format"
        ),
    Module(
        "check_for_log",
        SignalFile_Postprocess, SignalFile_Postprocess,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", BASED_ON_DATA, ["yes", "no"]),
        help="check for log in SignalFile_Postprocess"
        ),
    Module(
        "filter_and_threshold_genes",
        SignalFile_Postprocess,SignalFile_Postprocess,
        Constraint('format',MUST_BE,"tdf"),
        Constraint('logged',MUST_BE,"no"),
        Constraint('predataset',MUST_BE,"no"),
        Consequence('format',SAME_AS_CONSTRAINT),
        Consequence('logged',SAME_AS_CONSTRAINT),
        Consequence('predataset',SET_TO,'yes'),
        help="filter genes by a threshold using genepattern module"),
    Module(
        "log_signal",
        SignalFile_Postprocess, SignalFile_Postprocess,
        Constraint("format", MUST_BE, "tdf"),
        Constraint("logged", MUST_BE, "no"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("logged", SET_TO, "yes"),
        help="log SignalFile_Postprocess"),
    
    #impute
    Module(
        "convert_postprocess_impute",
        SignalFile_Postprocess, SignalFile_Impute,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF,CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SET_TO, 'zero_fill'),
        Consequence("missing_values", SET_TO, "unknown"),
        Consequence("filter", SET_TO, 'no'),
        help="convert SignalFile_Postprocess to SignalFile_Impute"),
        
    Module(
        "check_for_missing_values",
        SignalFile_Impute, SignalFile_Impute,
        Constraint("missing_values", MUST_BE, "unknown"),
        Consequence("missing_values", BASED_ON_DATA, ["no", "yes"]),
        help="check missing values in SignalFile_Impute"
        ),
    
    Module(
        "filter_genes_by_missing_values",
        SignalFile_Impute, SignalFile_Impute,
        OptionDef("filter_value", 0.50,
                     help="filter by missing values in percentage, etc.(0-1)"),
        Constraint("missing_values", MUST_BE, "yes"),
        Constraint("filter", MUST_BE, "no"),
        Consequence("missing_values", SAME_AS_CONSTRAINT),
        Consequence("filter", SET_TO, "yes"),
        help="filter genes by missing values in SignalFile_Impute"
        ),
    Module(
        "fill_missing_with_zeros",
        SignalFile_Impute, SignalFile_Impute,
        Constraint("missing_values", MUST_BE, "yes"),
        Consequence("missing_values", SET_TO, "no"),
        Consequence(
            "missing_algorithm", SET_TO, "zero_fill", side_effect=True),
        help="fill missing values in SignalFile_Impute with zeros"
        ),
    Module(
        "fill_missing_with_median",
        SignalFile_Impute,SignalFile_Impute,
        Constraint('missing_algorithm',MUST_BE,'zero_fill'),
        Constraint('missing_values',MUST_BE,'yes'),
        Consequence('missing_algorithm',SET_TO,"median_fill",side_effect=True),
        Consequence('missing_values',SET_TO,'no'),
        help="fill missing values in SignalFile_Impute with median"),
    #merge
    Module(
        "convert_impute_merge",
        SignalFile_Impute, SignalFile_Merge,
        Constraint("preprocess", CAN_BE_ANY_OF, ["unknown", "illumina", "agilent", "mas5", "loess","tcga","rsem"]),
        Constraint("contents", CAN_BE_ANY_OF,CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SET_TO,'no'),
        Consequence("dwd_norm", SET_TO,'no'),
        Consequence("combat_norm", SET_TO,'no'),
        Consequence("bfrm_norm",SET_TO,'no'),
        Consequence("shiftscale_norm", SET_TO,'no'),
        help="convert SignalFile_Impute to SignalFile_Merge"),
    
    Module(
        "convert_impute_merge_rma",
        SignalFile_Impute, SignalFile_Merge,
        Constraint("preprocess",
            MUST_BE,"rma"),
        Constraint("contents", CAN_BE_ANY_OF,CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_values",MUST_BE,"no"),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SET_TO,'yes'),
        Consequence("dwd_norm", SET_TO,'no'),
        Consequence("combat_norm", SET_TO,'no'),
        Consequence("bfrm_norm",SET_TO,'no'),
        Consequence("shiftscale_norm", SET_TO,'no'),
        help="convert SignalFile_Impute to SignalFile_Merge with preprocess=rma"),
    
    Module(  #did not consider the diff_expr case 
        "merge_two_classes_rma", [SignalFile_Merge, SignalFile_Merge], SignalFile_Merge, 
         Constraint("contents", MUST_BE, "class0", 0),
         Constraint("preprocess",MUST_BE,'rma',0),
         Constraint("combat_norm",MUST_BE,'no',0),
         Constraint("quantile_norm",MUST_BE,'yes',0), 
         Constraint("dwd_norm",MUST_BE,"no",0),
         Constraint("bfrm_norm",MUST_BE,"no",0),
         Constraint("shiftscale_norm",MUST_BE,"no",0),
         Constraint("preprocess",MUST_BE,'rma',1),
         Constraint("contents", MUST_BE, "class1", 1),
         Constraint("combat_norm",MUST_BE,'no',1),
         Constraint("quantile_norm",MUST_BE,"yes",1),
         Constraint("dwd_norm",MUST_BE,"no",1),
         Constraint("bfrm_norm",MUST_BE,"no",1),
         Constraint("shiftscale_norm",MUST_BE,"no",1),
         Consequence("contents", SET_TO, "class0,class1"),
         Consequence("preprocess",SAME_AS_CONSTRAINT,0),
         Consequence("combat_norm",SAME_AS_CONSTRAINT,0),
         Consequence("quantile_norm",SAME_AS_CONSTRAINT,0),
         Consequence("dwd_norm",SAME_AS_CONSTRAINT,0),
         Consequence("bfrm_norm",SAME_AS_CONSTRAINT,0),
         Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,0),
         DefaultAttributesFrom(0),
         DefaultAttributesFrom(1),
        help="merge two classes SignalFile_Merge with preprocess=rma,generate SignalFile_Merge"
        ),

    Module(  #did not consider the diff_expr case 
        "merge_two_classes", [SignalFile_Merge, SignalFile_Merge], SignalFile_Merge,
         Constraint("preprocess",CAN_BE_ANY_OF, ["unknown", "illumina",
                                                 "agilent", "mas5", "loess","tcga"]),
         Constraint("contents", MUST_BE, "class0", 0),
         Constraint("combat_norm",MUST_BE,'no',0),
         Constraint("quantile_norm",MUST_BE,'no',0), 
         Constraint("dwd_norm",MUST_BE,"no",0),
         Constraint("bfrm_norm",MUST_BE,"no",0),
         Constraint("shiftscale_norm",MUST_BE,"no",0),
         Constraint("contents", MUST_BE, "class1", 1),
         Constraint("combat_norm",MUST_BE,'no',1),
         Constraint("quantile_norm",MUST_BE,"no",1),
         Constraint("dwd_norm",MUST_BE,"no",1),
         Constraint("bfrm_norm",MUST_BE,"no",1),
         Constraint("shiftscale_norm",MUST_BE,"no",1),
         Consequence("contents", SET_TO, "class0,class1"),
         Consequence("preprocess",SAME_AS_CONSTRAINT,0),
         Consequence("combat_norm",SAME_AS_CONSTRAINT,0),
         Consequence("quantile_norm",SAME_AS_CONSTRAINT,0),
         Consequence("dwd_norm",SAME_AS_CONSTRAINT,0),
         Consequence("bfrm_norm",SAME_AS_CONSTRAINT,0),
         Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,0),
         DefaultAttributesFrom(0),
         DefaultAttributesFrom(1),
         help="merge two classes SignalFile_Merge,generate SignalFile_Merge"
         ),
    
        Module(
        "normalize_samples_with_quantile",
        SignalFile_Merge, SignalFile_Merge,
        Constraint("combat_norm", MUST_BE, "no"),
        Constraint("shiftscale_norm", MUST_BE, "no"),
        Constraint("bfrm_norm", MUST_BE, "no"),
        Constraint("dwd_norm", MUST_BE, "no"),
        Constraint("quantile_norm", MUST_BE, "no"),
        Consequence("quantile_norm", SET_TO, "yes"),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        help="nommalize SignalFile_Merge with quantile method"
        ),
    Module(
        "normalize_samples_with_bfrm",
        SignalFile_Merge,SignalFile_Merge,
        OptionDef("num_factors",1,help="num factors for bfrm normalization"),
        Constraint('bfrm_norm',MUST_BE,"no"),
        Constraint('combat_norm',MUST_BE,"no"),
        Constraint('shiftscale_norm',MUST_BE,"no"),
        Constraint('dwd_norm',MUST_BE,"no"),
        Consequence('bfrm_norm',SET_TO,'yes'),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        help="nommalize SignalFile_Merge with bfrm method"
        ),

    Module(
        "normalize_samples_with_combat",  
        [ClassLabelFile,SignalFile_Merge],SignalFile_Merge,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("combat_norm",MUST_BE,"no",1),
        Constraint('shiftscale_norm',MUST_BE,"no",1),
        Constraint('bfrm_norm',MUST_BE,"no",1),
        Constraint('dwd_norm',MUST_BE,"no",1),
        Consequence("combat_norm",SET_TO,"yes"),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT,1),
        Consequence("dwd_norm",SAME_AS_CONSTRAINT,1),
        Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,1),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with combat method"),
           
    Module(
        "normalize_samples_with_dwd",  
        [ClassLabelFile,SignalFile_Merge],SignalFile_Merge,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS,0),
        Constraint("dwd_norm",MUST_BE,"no",1),
        Consequence("dwd_norm",SET_TO,"yes"),
        Constraint("bfrm_norm",MUST_BE,"no",1),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT,1),
        Constraint("combat_norm",MUST_BE,"no",1),
        Consequence("combat_norm",SAME_AS_CONSTRAINT,1),
        Constraint("shiftscale_norm",MUST_BE,"no",1),
        Consequence("shiftscale_norm",SAME_AS_CONSTRAINT,1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with dwd method"),

    Module(
        "normalize_samples_with_shiftscale",  
        [ClassLabelFile,SignalFile_Merge],SignalFile_Merge,
        Constraint("cls_format",MUST_BE,'cls',0),
        Constraint("contents",CAN_BE_ANY_OF,CONTENTS,0),
        Constraint("shiftscale_norm",MUST_BE,"no",1),
        Consequence("shiftscale_norm",SET_TO,"yes"),
        Constraint("bfrm_norm",MUST_BE,"no",1),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT,1),
        Constraint("combat_norm",MUST_BE,"no",1),
        Consequence("combat_norm",SAME_AS_CONSTRAINT,1),
        Constraint("dwd_norm",MUST_BE,"no",1),
        Consequence("dwd_norm",SAME_AS_CONSTRAINT,1),
        Constraint("contents",SAME_AS,0,1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        DefaultAttributesFrom(1),
        help="nommalize SignalFile_Merge with shiftscale method"),
    ###normalize
    Module(
        "convert_merge_normalize",
        SignalFile_Merge, SignalFile_Normalize,
        Constraint("preprocess",CAN_BE_ANY_OF, PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF,CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("combat_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SET_TO,'unknown'),
        Consequence("gene_normalize", SET_TO,'unknown'),
        Consequence("format", SET_TO,'tdf'),
        help="convert SignalFile_Merge to SignalFile_Normalize"
        ),
    Module(
        "check_gene_center",
        SignalFile_Normalize,SignalFile_Normalize,
        Constraint("format",MUST_BE,'tdf'),
        Constraint("gene_center",MUST_BE,"unknown"),
        Constraint("gene_normalize", MUST_BE, "unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence("gene_center",BASED_ON_DATA,["no", "mean", "median"]),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        help="check SignalFile_Normalize if it is gene center or not"),
        
    Module(
        "check_gene_normalize",
        SignalFile_Normalize,SignalFile_Normalize,
        Constraint("format",MUST_BE,'tdf'),
        Constraint("gene_normalize",MUST_BE,"unknown"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Consequence(
            "gene_normalize", BASED_ON_DATA,
            ["no", "variance", "sum_of_squares"]),
        help="check SignalFile_Normalize if gene normalize or nto"),
        
    Module(   
        "convert_signal_to_pcl",
        SignalFile_Normalize,SignalFile_Normalize,
        Constraint("format",MUST_BE,'tdf'),
        Consequence("format",SET_TO,'pcl'),
        help="convert SignalFile_Normalize from tdf format to pcl format"),
    
    Module(
        "center_genes",
        SignalFile_Normalize,SignalFile_Normalize,
        Constraint("format",MUST_BE,"pcl"),
        Constraint("gene_center",MUST_BE,"no"),
        Constraint("gene_normalize",MUST_BE,'unknown'),
        Consequence("format",SET_TO,"tdf"),
        Consequence("gene_center",SET_TO_ONE_OF,["mean", "median"]),
        Consequence("gene_normalize",SAME_AS_CONSTRAINT),
        help="center genes in SignalFile_Normalize"),
    
    Module(
        "normalize_genes",
        SignalFile_Normalize,SignalFile_Normalize,
        Constraint("format",MUST_BE,"pcl"),
        Constraint("gene_normalize",MUST_BE,"no"),
        Consequence("format",SET_TO,"tdf"),
        Consequence("gene_normalize",SET_TO_ONE_OF,["variance", "sum_of_squares"]),
        help="normalize genes in SignalFile_Normalize"),
    ##Order
    Module(
        "convert_normalize_order",
        SignalFile_Normalize, SignalFile_Order,
        Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Constraint("format",MUST_BE,'tdf'),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("combat_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median","mean"]),
        Constraint("gene_normalize", CAN_BE_ANY_OF, ["no", "variance","sum_of_squares"]),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order", SET_TO,'no'),
        help="convert SignalFile_Normalize to SignalFile_Order"
        ),


    Module(  
        "rank_genes_by_class_neighbors",
        [ClassLabelFile,SignalFile_Order],GeneListFile,
        OptionDef("cn_num_neighbors",50,help='number of neighbors for class neighbors method'),
        OptionDef("cn_num_perm",100,help='number of permutation for class neighbors method'),
        OptionDef("cn_user_pval",0.5,help='number of user p value for class neighbors method'),
        OptionDef("cn_min_threshold",10,help='min threshold for class neighbors method'),
        OptionDef("cn_max_threshold",16000,help='max threshold for class neighbors method'),
        OptionDef("cn_min_folddiff",5,help='min fold diff for class neighbors method'),
        OptionDef("cn_abs_diff",50,help='abs diff for class neighbors method'),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint(
            "contents", CAN_BE_ANY_OF, CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Constraint("gene_order", MUST_BE,"no",1),
        Consequence("gene_order",SET_TO,"class_neighbors"),
        Consequence("cn_mean_or_median",SET_TO_ONE_OF, ['mean', 'median']),
        Consequence("cn_ttest_or_snr",SET_TO_ONE_OF, ['t_test','snr']),
        Consequence("cn_filter_data",SET_TO_ONE_OF, ['yes','no']),
        Consequence("contents",SAME_AS_CONSTRAINT,0),
        help="rank the genes in SignalFile_Order by class neighbors method"
        ),
    
    Module(
        "rank_genes_by_sample_ttest",
        [ClassLabelFile, SignalFile_Order], GeneListFile,
        OptionDef("gene_select_threshold", 0.05,help="threshold for sample ttest"),
        Constraint("cls_format", MUST_BE, 'cls', 0),
        Constraint("gene_order", MUST_BE, "no", 1),
        Constraint(
            "contents", CAN_BE_ANY_OF, CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Consequence("gene_order",SET_TO_ONE_OF, ["t_test_p", "t_test_fdr"]),
        Consequence("contents",SAME_AS_CONSTRAINT,0),
        help="rank the genes in SignalFile_Order by ttest method"
       ),
         
    Module(  
         "reorder_genes",  
         [GeneListFile,SignalFile_Order], SignalFile_Order,
         Constraint("gene_order", CAN_BE_ANY_OF, ['t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list"],0),
         Constraint("gene_order", MUST_BE,"no",1),
         Constraint(
            "preprocess", CAN_BE_ANY_OF, PREPROCESS, 1),
         Constraint(
            "contents", CAN_BE_ANY_OF, CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Constraint("gene_center",CAN_BE_ANY_OF,['median','mean','no'],1),
         Constraint("gene_normalize",CAN_BE_ANY_OF,['sum_of_squares','variance','no'],1),
         Consequence("gene_center",SAME_AS_CONSTRAINT,1),
         Consequence("gene_normalize",SAME_AS_CONSTRAINT,1),
         Consequence("gene_order",SAME_AS_CONSTRAINT,0),
         Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         DefaultAttributesFrom(1),
         help="rank the genes in SignalFile_Order by genes in GeneListFile"
        ),
    Module(  
         "reorder_genes_with_diff",  
         [GeneListFile,SignalFile_Order], SignalFile_Order,
         Constraint("gene_order", CAN_BE_ANY_OF, ['diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change'],0),
         Constraint("gene_order", MUST_BE,"no",1),
         Constraint(
            "preprocess", CAN_BE_ANY_OF, PREPROCESS, 1),
         Constraint(
            "contents", CAN_BE_ANY_OF, CONTENTS,0),
         Constraint("contents",CAN_BE_ANY_OF,CONTENTS,1),
         Constraint("gene_center",CAN_BE_ANY_OF,['median','mean','no'],1),
         Constraint("gene_normalize",CAN_BE_ANY_OF,['sum_of_squares','variance','no'],1),
         Consequence("gene_center",SAME_AS_CONSTRAINT,1),
         Consequence("gene_normalize",SAME_AS_CONSTRAINT,1),
         Consequence("gene_order",SAME_AS_CONSTRAINT,0),
         Consequence("preprocess", SAME_AS_CONSTRAINT, 1),
         Consequence("contents",SAME_AS_CONSTRAINT,1),
         
         DefaultAttributesFrom(1),
         help="reorder SignalFile_Order with genes in GeneListFile"
        ),
    ##Annotate
    Module(
        "convert_order_annotate",
        SignalFile_Order, SignalFile_Annotate,
        Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("combat_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median","mean"]),
        Constraint("gene_normalize", CAN_BE_ANY_OF, ["no", "variance","sum_of_squares"]),
        Constraint("gene_order", CAN_BE_ANY_OF,["no",'t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change']),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order",SAME_AS_CONSTRAINT),
        Consequence("annotate",SET_TO,"no"),
        Consequence("rename_sample",SET_TO,"no"),
        Consequence("platform",SET_TO,"no"),
        help="convert SignalFile_Order to SignalFile_Annotate"
        ),
    Module(
         'annotate_probes',
         SignalFile_Annotate,SignalFile_Annotate,
         Constraint("annotate", MUST_BE,"no"),
         Constraint("platform", MUST_BE,"no"),
         Consequence("annotate",SET_TO,"yes"),
         Consequence("platform",SAME_AS_CONSTRAINT),
         help="annotate SignalFile_Annotate"),
    
    Module(
         'download_GEO_family_soft',
         GEOSeries,GEOfamily,
         OptionDef("GSEID",help='GSEID for download family_soft file'),
         Constraint("contents",CAN_BE_ANY_OF, CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help="download geo family soft file"),
    
    Module(
         'convert_family_soft_to_rename',
         GEOfamily,RenameFile,
         OptionDef("GSEID",help='GSEID for download family_soft file'),
         Constraint("contents",CAN_BE_ANY_OF, CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("labels_from",SET_TO_ONE_OF,["title","description"]),
         help="convert famliy soft file to RenameFile"),
    
    Module( 
       "relabel_samples",  
        [RenameFile,SignalFile_Annotate],SignalFile_Annotate,
         Constraint("rename_sample",MUST_BE,"no",1),
         Constraint("labels_from",CAN_BE_ANY_OF,["title","description"],0),
         Constraint("contents",CAN_BE_ANY_OF,CONTENTS,0),
         Constraint("contents",SAME_AS,0,1),
         Consequence("rename_sample",SET_TO,"yes"),
         Consequence("contents", SAME_AS_CONSTRAINT,0),
         Constraint('annotate',MUST_BE,'no',1),
         Consequence("annotate",SAME_AS_CONSTRAINT,1),
         Constraint('platform',MUST_BE,'no',1),
         Consequence("platform",SAME_AS_CONSTRAINT,1),
         DefaultAttributesFrom(1),
         help="relabel the sample names in SignalFile_Annotate given RenameFile"
       ),
    Module(
         'add_crossplatform_probeid',
         SignalFile_Annotate,SignalFile_Annotate,
         OptionDef("platform_name",help="given the new platform name to add to the file"),
         Constraint("platform", MUST_BE,"no"),
         Consequence("platform",SET_TO,"yes"),
         help="add a cross platform to SignalFile_Annotate"),
    
    #Filter
    Module(
        "convert_annotate_filter",
        SignalFile_Annotate, SignalFile_Filter,
        Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("combat_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median","mean"]),
        Constraint("gene_normalize", CAN_BE_ANY_OF, ["no", "variance","sum_of_squares"]),
        Constraint("gene_order", CAN_BE_ANY_OF,["no",'t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change']),
        Constraint("annotate",CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("rename_sample",CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("platform",CAN_BE_ANY_OF, ["no", "yes"]),
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order",SAME_AS_CONSTRAINT),
        Consequence("annotate",SAME_AS_CONSTRAINT),
        Consequence("rename_sample",SAME_AS_CONSTRAINT),
        Consequence("platform",SAME_AS_CONSTRAINT),
        Consequence("logged",SET_TO,"yes"),
        Consequence("format",SET_TO,"tdf"),
        Consequence("num_features",SET_TO,"no"),
        Consequence("unique_genes",SET_TO,"no"),
        Consequence("duplicate_probe",SET_TO,"no"),
        Consequence("group_fc",SET_TO,"no"),
        help="convert SignalFile_Annotate to SignalFile_Filter"
        ),
    Module(
        'remove_duplicate_genes',
        SignalFile_Filter, SignalFile_Filter,
        Constraint("annotate", MUST_BE,"yes"),
        Constraint("num_features", MUST_BE,"no"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Constraint("unique_genes", MUST_BE,'no'),
        Consequence("annotate", SAME_AS_CONSTRAINT),
        Consequence("num_features", SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe", SAME_AS_CONSTRAINT),
        Constraint("format", MUST_BE,"tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE,"yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Consequence(
            "unique_genes", SET_TO_ONE_OF,
            ['average_genes', 'high_var', 'first_gene']),
        help="remove duplicate genes in SignalFile_Filter"),
      
    Module(
         'select_first_n_genes',
        SignalFile_Filter,SignalFile_Filter,
        OptionDef("num_features_value",500,help="num of features to be selected in the SignalFile"),
        Constraint("format", MUST_BE,"tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE,"yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("num_features", MUST_BE,"no"),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("num_features",SET_TO,"yes"),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT),
         help="select first n genes in SignalFile_Filter"), 
      
     Module(
        'remove_duplicate_probes',
        SignalFile_Filter,SignalFile_Filter,
        Constraint("format", MUST_BE,"tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE,"yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("duplicate_probe",SET_TO,'high_var_probe'),
        Constraint("platform", MUST_BE,"yes"),
        Consequence("platform",SAME_AS_CONSTRAINT),
        help="remove duplciate probes in SignalFile_Filter by high_var_probe method"),
    Module(
         'select_probe_by_best_match',
        SignalFile_Filter,SignalFile_Filter,
        Constraint("format", MUST_BE,"tdf"),
        Consequence("format", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE,"yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("duplicate_probe", MUST_BE,'no'),
        Consequence("duplicate_probe",SET_TO,'closest_probe'),
        Constraint("platform", MUST_BE,"yes"),
        Consequence("platform",SAME_AS_CONSTRAINT),
         help="remove duplciate probes in SignalFile_Filter by closest_probe method"),
    
    Module(
        "filter_genes_by_fold_change_across_classes",
        [ClassLabelFile,SignalFile_Filter],SignalFile_Filter,
        OptionDef("group_fc_num",help="group fold change number"),
        Constraint("cls_format", MUST_BE,'cls',0),
        Constraint("format", MUST_BE,"tdf",1),
        Consequence("format", SAME_AS_CONSTRAINT,1),
        Constraint("logged", MUST_BE,"yes",1),
        Consequence("logged", SAME_AS_CONSTRAINT,1),
        Constraint("group_fc", MUST_BE,"no",1),
        Constraint("num_features", MUST_BE,"no",1),
        Constraint("duplicate_probe", MUST_BE,"no",1),
        Constraint("unique_genes", MUST_BE,"no",1),
        Constraint("contents",CAN_BE_ANY_OF, CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Consequence("contents",SAME_AS_CONSTRAINT,0),
        Consequence("num_features",SAME_AS_CONSTRAINT,1),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT,1),
        Consequence("unique_genes",SAME_AS_CONSTRAINT,1),
        Consequence("group_fc",SET_TO,"yes"),
        DefaultAttributesFrom(1),
        help="filter genes in SignalFile_Filter by fold change in different classes"),
    Module(   
        "convert_signal_to_gct",
        SignalFile_Filter,SignalFile_Filter,
        Constraint("format",MUST_BE,'tdf'),
        Consequence("format",SET_TO,'gct'),
        help="convert SignalFile_Filter in tdf format to gct format"),
    Module( 
        'unlog_signal',
        SignalFile_Filter,SignalFile_Filter,
        Constraint("format", MUST_BE,"tdf"),
        Constraint("logged", MUST_BE,"yes"),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("logged",SET_TO,"no"),
        help="unlog SignalFile_Filter"),
    
    Module(
        "convert_label_to_cls",   
        [ClassLabelFile,SignalFile_Merge],ClassLabelFile,
        Constraint("cls_format",MUST_BE,'label',0),
        Constraint("contents",CAN_BE_ANY_OF, CONTENTS,0),
        Constraint("contents",SAME_AS,0,1),
        Consequence("contents", SAME_AS_CONSTRAINT,0),
        Consequence("cls_format",SET_TO,'cls'),
        help="convert ClassLabelFile with label format to cls format"
        ),
    Module( 
        'transfer',
        SignalFile_Filter,SignalFile,
        Constraint("preprocess", CAN_BE_ANY_OF, PREPROCESS),
        Constraint("contents", CAN_BE_ANY_OF, CONTENTS),
        Constraint("predataset", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("filter",CAN_BE_ANY_OF,["no", "yes"]),
        Constraint("missing_algorithm",CAN_BE_ANY_OF,['none','zero_fill','median_fill']),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Consequence("predataset", SAME_AS_CONSTRAINT),
        Consequence("missing_algorithm",SAME_AS_CONSTRAINT),
        Consequence("filter", SAME_AS_CONSTRAINT),
        Constraint("quantile_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("combat_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("shiftscale_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("bfrm_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("dwd_norm", CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("gene_center", CAN_BE_ANY_OF, ["no", "median","mean"]),
        Constraint("gene_normalize", CAN_BE_ANY_OF, ["no", "variance","sum_of_squares"]),
        Constraint("gene_order", CAN_BE_ANY_OF,["no",'t_test_p', "t_test_fdr",
                                   'class_neighbors', "gene_list",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change']),
        Constraint("annotate",CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("rename_sample",CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("platform",CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("logged",CAN_BE_ANY_OF, ["no", "yes"]),
        Constraint("format",CAN_BE_ANY_OF,['tdf','gct']),
        Constraint("num_features",CAN_BE_ANY_OF,['yes',"no"]),
        Constraint("unique_genes",CAN_BE_ANY_OF,["no", "average_genes", "high_var", "first_gene"]),
        Constraint("duplicate_probe",CAN_BE_ANY_OF,["no", "closest_probe", "high_var_probe"]),
        Constraint("group_fc",CAN_BE_ANY_OF,['yes',"no"]),
        
        Consequence("dwd_norm", SAME_AS_CONSTRAINT),
        Consequence("combat_norm", SAME_AS_CONSTRAINT),
        Consequence("bfrm_norm",SAME_AS_CONSTRAINT),
        Consequence("shiftscale_norm", SAME_AS_CONSTRAINT),
        Consequence("quantile_norm", SAME_AS_CONSTRAINT),
        Consequence("gene_center", SAME_AS_CONSTRAINT),
        Consequence("gene_normalize", SAME_AS_CONSTRAINT),
        Consequence("gene_order",SAME_AS_CONSTRAINT),
        Consequence("annotate",SAME_AS_CONSTRAINT),
        Consequence("rename_sample",SAME_AS_CONSTRAINT),
        Consequence("platform",SAME_AS_CONSTRAINT),
        Consequence("logged",SAME_AS_CONSTRAINT),
        Consequence("format",SAME_AS_CONSTRAINT),
        Consequence("num_features",SAME_AS_CONSTRAINT),
        Consequence("unique_genes",SAME_AS_CONSTRAINT),
        Consequence("duplicate_probe",SAME_AS_CONSTRAINT),
        Consequence("group_fc",SAME_AS_CONSTRAINT),
        help="transfer SignalFile_Filter to SignalFile"
        )
    ]

list_files=[RenameFile,AgilentFiles,CELFiles,ControlFile,ExpressionFiles,
            GPRFiles,GEOSeries,IDATFiles,ClassLabelFile,ILLUFolder,GeneListFile,
           SignalFile,SignalFile_Postprocess,SignalFile_Impute, SignalFile_Merge,
            SignalFile_Normalize,SignalFile_Order,SignalFile_Annotate,SignalFile_Filter,GEOfamily,
            TCGAID,TCGAFile]

