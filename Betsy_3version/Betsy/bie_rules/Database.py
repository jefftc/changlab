#Database

from Betsy.bie3 import *
CONTENTS = ["train0", "train1","test", "class0,class1,test",
             "class0", "class1", "class0,class1","unspecified",
            "diff_class0","diff_class1","diff_class0,diff_class1","diff_unspecified"]

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

list_files = [GEOSeries,GEOfamily,TCGAID,TCGAFile]

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
    
    Module(
         'download_GEO_family_soft',
         GEOSeries,GEOfamily,
         OptionDef("GSEID",help='GSEID for download family_soft file'),
         Constraint("contents",CAN_BE_ANY_OF, CONTENTS),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help="download geo family soft file"),
    
    
    ]
