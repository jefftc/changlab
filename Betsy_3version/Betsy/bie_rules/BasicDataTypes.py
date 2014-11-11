#BasicDataTypes

from Betsy.bie3 import *
import Database


RenameFile = DataType(
    'RenameFile',
    AttributeDef("contents",Database.CONTENTS,'unspecified','unspecified',help="contents"),
    AttributeDef("labels_from",["title","description"],'title','title',help="labels from title or description"),
    help="A file used to rename the sample name in the gene expression file.")

ExpressionFiles = DataType("ExpressionFiles",
                           AttributeDef("contents",Database.CONTENTS,
                                        'unspecified','unspecified',help="contents"),
                           AttributeDef("filetype",['unknown','cel','gpr','idat','agilent','matrix'],
                                        'unknown','unknown',help="filetype"),
                           help="Expression file folder, can be CELFiles, IDATFiles,"\
                                 "AgilentFile,GPRFiles")


GeneListFile=DataType(
    "GeneListFile", 
    AttributeDef('cn_mean_or_median',['mean', 'median'], 'mean','mean',help="class neighbors mean or median"),
    AttributeDef('cn_ttest_or_snr',['t_test','snr'], 't_test','t_test',help="class neighbors ttest or snr"),
    AttributeDef('cn_filter_data',['yes','no'], 'no','no',help="class neighbors filter data or not"),
    AttributeDef('gene_order',['no', "gene_list", "class_neighbors",
                               "t_test_p", "t_test_fdr",'diff_ttest','diff_sam',
                               'diff_ebayes','diff_fold_change'],
                 't_test_p',"t_test_p",help="gene order method"),
    AttributeDef("contents",Database.CONTENTS,'unspecified','unspecified',help="contents"),
    help="A file contains a list of genes.")


list_files = [RenameFile,ExpressionFiles,GeneListFile]

all_modules = [
 Module(
        "download_geo_supplement", Database.GEOSeries, ExpressionFiles,
         OptionDef("GSEID",help="GSEID to download"),
         OptionDef("GPLID","",help="GPDID to download"),
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("filetype",SET_TO,'unknown'),
        help="download GEO data from geo website according to GSEID and GPLID"
        ),
Module("download_geo_seriesmatrix", Database.GEOSeries, Database.MatrixFile,
         OptionDef("GSEID",help="GSEID to download"),
         OptionDef("GPLID","",help="GPDID to download"),
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
        help="download GEO data from geo website according to GSEID and GPLID"
        ),

    Module(
         'convert_family_soft_to_rename',
         Database.GEOfamily,RenameFile,
         OptionDef("GSEID",help='GSEID for download family_soft file'),
         Constraint("contents",CAN_BE_ANY_OF, Database.CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         Consequence("labels_from",SET_TO_ONE_OF,["title","description"]),
         help="convert famliy soft file to RenameFile"),

     Module(
        "check_geo_file_type", [ExpressionFiles,Database.MatrixFile], ExpressionFiles,
         Constraint("contents",CAN_BE_ANY_OF,Database.CONTENTS),
         Constraint("contents",SAME_AS,0,1),
         Consequence("contents",SAME_AS_CONSTRAINT,0),
         Constraint("filetype",MUST_BE,'unknown',0),
         Consequence("filetype",BASED_ON_DATA,['matrix','cel','gpr','idat','agilent']),
        help="check the file type downloaded from geo database"
        ),

]
