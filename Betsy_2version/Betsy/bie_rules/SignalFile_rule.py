# SignalFile
from Betsy import bie
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

RenameFile = bie.DataType(
    'RenameFile',
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),)
AgilentFiles = bie.DataType(
    "AgilentFiles",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),)
CELFiles = bie.DataType(
    "CELFiles",
    bie.Attribute(
        version=["unknown", "cc", "v3", "v4"],
        DEFAULT="unknown"),
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    )
ControlFile = bie.DataType(
    "ControlFile",
    bie.Attribute(
        preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        DEFAULT="unknown"),
    bie.Attribute(
        missing_values=["unknown", "no", "yes"],
        DEFAULT="no"),
    bie.Attribute(
        missing_algorithm=["none", "median_fill", "zero_fill"],
        DEFAULT="none", OPTIONAL=True),
    bie.Attribute(
        logged=["unknown", "no", "yes"],
        DEFAULT="no"),
    bie.Attribute(
        format=["unknown", "tdf", "gct", "jeffs", "pcl", "res", "xls"],
        DEFAULT="gct"),
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    )
ExpressionFiles = bie.DataType(
    "ExpressionFiles",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    )
GPRFiles = bie.DataType(
    "GPRFiles",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),)
GEOSeries = bie.DataType(
    "GEOSeries",
    bie.Attribute(GSEID=bie.ANYATOM, DEFAULT=bie.ANYATOM),
    bie.Attribute(GPLID=bie.ANYATOM, DEFAULT=bie.ANYATOM,OPTIONAL=True),
    )
IDATFiles = bie.DataType(
    "IDATFiles",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),)
ClassLabelFile = bie.DataType(
    "ClassLabelFile",
    bie.Attribute(
    contents=["train0", "train1", "test", "class0,class1,test",
                  "class0", "class1", "class0,class1",
                  "no"],DEFAULT='no'),
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(cls_format=['cls','label','unknown'],DEFAULT='unknown')
    )

ILLUFolder = bie.DataType(
    "ILLUFolder", 
    bie.Attribute(
        illu_manifest=ILLU_MANIFEST,
        DEFAULT='HumanHT-12_V4_0_R2_15002873_B.txt'),
    bie.Attribute(
        illu_chip=ILLU_CHIP,
        DEFAULT='ilmn_HumanHT_12_V4_0_R1_15002873_B.chip'),
    bie.Attribute(illu_bg_mode=['false', 'true'], DEFAULT="false"),

    bie.Attribute(illu_coll_mode=['none', 'max', 'median'], DEFAULT="none"),
    bie.Attribute(illu_clm=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(illu_custom_chip=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(illu_custom_manifest=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),)

GeneListFile=bie.DataType(
    "GeneListFile",
    bie.Attribute(cn_num_neighbors=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(cn_num_perm=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(cn_user_pval=bie.ANYATOM, DEFAULT=""),  
    bie.Attribute(cn_mean_or_median=['mean', 'median'], DEFAULT='mean'),
    bie.Attribute(cn_ttest_or_snr=['t_test','snr'], DEFAULT='t_test'),
    bie.Attribute(cn_filter_data=['yes','no'], DEFAULT='no'),
    bie.Attribute(cn_min_threshold=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(cn_max_threshold=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(cn_min_folddiff=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(cn_abs_diff=bie.ANYATOM, DEFAULT=""),
    bie.Attribute(gene_select_threshold=bie.ANYATOM,DEFAULT=""),
    bie.Attribute(gene_order=['no', "gene_list", "class_neighbors",
                          "t_test_p", "t_test_fdr"], DEFAULT='gene_list'),
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),)
SignalFile = bie.DataType(
    "SignalFile",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    # Properties of the format.
    bie.Attribute(
        format=["unknown", "tdf", "gct", "jeffs", "pcl", "res", "xls"],
        DEFAULT="unknown"),
    # Properties of the data.
    bie.Attribute(
        preprocess=["unknown", "illumina", "agilent", "mas5", "rma", "loess"],
        DEFAULT="unknown"),
    bie.Attribute(
        missing_values=["unknown", "no", "yes"],
        DEFAULT="unknown"),
    bie.Attribute(
        missing_algorithm=["none", "median_fill", "zero_fill"],
        DEFAULT="none", OPTIONAL=True),
    bie.Attribute(
        logged=["unknown", "no", "yes"],
        DEFAULT="unknown"),

    # Normalizing the data.  Very difficult to check normalization.
    # If you're not sure if the data is normalized, then the answer is
    # "no".
    bie.Attribute(dwd_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(bfrm_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(quantile_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(shiftscale_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(combat_norm=["no", "yes"], DEFAULT="no"),
    bie.Attribute(predataset=["no", "yes"], DEFAULT="no"),
    bie.Attribute(filter=bie.ANYATOM, DEFAULT="no"),
    bie.Attribute(rename_sample=["no", "yes"], DEFAULT="no"),
    bie.Attribute(
        contents=["train0", "train1", "test", 'class0,class1,test',
                  "class0", "class1", "class0,class1",
                  "no"],
        DEFAULT="no"))


list_files=[RenameFile,AgilentFiles,CELFiles,ControlFile,ExpressionFiles,
            GPRFiles,GEOSeries,IDATFiles,ClassLabelFile,ILLUFolder,GeneListFile,
            SignalFile]
            
all_modules = [
    # GSEID
    bie.Module("download_geo", GEOSeries, ExpressionFiles),
    bie.Module("extract_CEL_files", ExpressionFiles, CELFiles(version="unknown")),

    # CELFiles
    bie.QueryModule(
        "detect_CEL_version",
        CELFiles(version="unknown"),
        CELFiles(version=["cc", "v3", "v4"])),
    bie.Module(
        "convert_CEL_cc_to_CEL_v3",
        CELFiles(version="cc"),
        CELFiles(version="v3")),
    bie.Module(
        "preprocess_rma",
        CELFiles(version=["v3", "v4"]),
        SignalFile(
            logged="yes", preprocess="rma", format="jeffs",
            missing_values="no")
        ),
    bie.Module(
        "preprocess_mas5",
        CELFiles(version=["v3", "v4"]),
        SignalFile(
            logged="no", preprocess="mas5", format="jeffs",
            missing_values="no")),
    # IDATFiles
    bie.Module("extract_illumina_idat_files", ExpressionFiles, IDATFiles),
    bie.Module(
        "preprocess_illumina",
        IDATFiles,
        ILLUFolder(
            illu_manifest=ILLU_MANIFEST, illu_chip=ILLU_CHIP,
            illu_bg_mode=["false", "true"],
            illu_coll_mode=["none", "max", "median"],
            illu_clm=bie.ANYATOM, illu_custom_chip=bie.ANYATOM,
            illu_custom_manifest=bie.ANYATOM)),
    bie.Module(
        "get_illumina_signal",
        ILLUFolder,
        SignalFile(preprocess="illumina", format="gct",logged="no")
        ),
    bie.Module(
        "get_illumina_control",
        ILLUFolder,
        ControlFile(preprocess="illumina", format="gct", logged="no")
        ),
    
    # AgilentFiles
    bie.Module(
        "extract_agilent_files", ExpressionFiles, AgilentFiles),
    bie.Module(
        "preprocess_agilent",
        AgilentFiles,
        SignalFile(logged="unknown", preprocess="agilent", format="tdf")),

    # GPRFiles
    bie.Module(
        "extract_gpr_files", ExpressionFiles, GPRFiles),
    bie.Module(
        "normalize_with_loess",
        GPRFiles,
        SignalFile(format="tdf", logged="unknown", preprocess="loess")
        ),
    
    # SignalFile
    bie.Module(
        "convert_signal_to_tdf",
        SignalFile(
            format=['pcl', 'res', 'gct', 'jeffs', 'unknown', 'xls']),
        SignalFile(format='tdf')),
    bie.QueryModule(
        "check_for_log",
        SignalFile(format="tdf", logged='unknown'),
        SignalFile(format="tdf", logged=['yes', "no"])),
    bie.Module(
        "log_signal",
        SignalFile(logged='no', format="tdf"),
        SignalFile(logged='yes', format='tdf')),
    bie.Module(
        "fill_missing_with_median",
        SignalFile(
            format="tdf", logged="yes", missing_algorithm="none",
            missing_values="yes"),
        SignalFile(
            format="tdf", logged="yes", missing_algorithm="median_fill",
            missing_values="no")),
    bie.Module(
        "fill_missing_with_zeros",
        SignalFile(
            format="tdf", logged="yes", missing_algorithm="none",
            missing_values="yes"),
        SignalFile(
            format="tdf", logged="yes", missing_algorithm="zero_fill",
            missing_values="no")),
    bie.QueryModule(
        "check_for_missing_values",
        SignalFile(format="tdf", missing_values="unknown",logged="yes"),
        SignalFile(format="tdf", missing_values=["no", "yes"],logged="yes")),
    bie.Module(
        "filter_genes_by_missing_values",
        SignalFile(
            format='tdf', logged="yes", missing_values="yes", filter="no"),
        SignalFile(
            format='tdf', logged="yes", missing_values="yes", filter=bie.ANYATOM)),
    bie.Module(
        "filter_and_threshold_genes",
        SignalFile(format="tdf", logged="no", predataset="no"),
        SignalFile(format="tdf", logged="no", predataset="yes")),
    bie.Module(
        "relabel_samples",
        [RenameFile,
         SignalFile(format='tdf', rename_sample="no",logged="yes",
                    missing_values='no',
                   combat_norm='no',quantile_norm="no",
                   dwd_norm="no",bfrm_norm="no",shiftscale_norm="no")],
        SignalFile(format='tdf', rename_sample="yes",logged="yes",
                   missing_values='no',
                   combat_norm='no',quantile_norm="no",
                   dwd_norm="no",bfrm_norm="no",shiftscale_norm="no")),

    bie.Module(
           "merge_two_classes",
           [SignalFile(contents="class0",format='tdf',logged="yes",
                       missing_values='no',
                       combat_norm='no',quantile_norm="no",
                       dwd_norm="no",bfrm_norm="no",shiftscale_norm="no"),
            SignalFile(contents="class1",format='tdf',logged="yes",
                       missing_values='no',
                       combat_norm='no',quantile_norm="no",
                       dwd_norm="no",bfrm_norm="no",shiftscale_norm="no")],
            SignalFile(contents="class0,class1",format='tdf',logged="yes",
                       missing_values='no',
                       combat_norm='no',quantile_norm="no",
                       dwd_norm="no",bfrm_norm="no",shiftscale_norm="no")),
    
    # Sample normalization.
    bie.Module(
        "normalize_samples_with_quantile",
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            quantile_norm="no"),
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            quantile_norm="yes")),
    bie.Module(
        "normalize_samples_with_combat",
        [ClassLabelFile(cls_format='cls'),
            SignalFile(
            format="tdf", logged="yes", missing_values="no",
            combat_norm="no")],
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            combat_norm="yes")),
    bie.Module(
        "normalize_samples_with_dwd",
        [ClassLabelFile(cls_format='cls'),
         SignalFile(
            format="tdf", logged="yes", missing_values="no",
            dwd_norm="no")],
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            dwd_norm="yes")),
    bie.Module(
        "normalize_samples_with_bfrm",
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            bfrm_norm="no"),
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            bfrm_norm="yes")),
    bie.Module(
        "normalize_samples_with_shiftscale",
        [ClassLabelFile(cls_format='cls'),
            SignalFile(
            format="tdf", logged="yes", missing_values="no",
            shiftscale_norm="no")],
        SignalFile(
            format="tdf", logged="yes", missing_values="no",
            shiftscale_norm="yes")),
    bie.Module(
        "convert_label_to_cls",
        [ClassLabelFile(cls_format='label'),
         SignalFile(format='tdf',logged='yes')],
         ClassLabelFile(cls_format='cls'))
    ]
