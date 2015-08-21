from Betsy.bie3 import *
import BasicDataTypes as BDT
import GeneExpProcessing

ConsensusClusteringFolder = DataType('ConsensusClusteringFolder',
                    AttributeDef("contents",BDT.CONTENTS,
                               'unspecified','unspecified',help="contents"),
                    AttributeDef("Consensus_algorithm",["Hierarchical", "SOM","NMF","KMeans"],
                                 "Hierarchical","Hierarchical",help="clustering algorithm in ConsensusClustering"),
                    AttributeDef("clusterby",['columns', 'rows'],
                        "columns","columns",help="cluster by columns or rows"),
                    AttributeDef("cc_distance",['Euclidean', 'Pearson'],
                        "Euclidean","Euclidean",help="distance measure in consensusClustering"),
                    AttributeDef("merge_type",['average', 'complete','single'],
                        "average","average",help="merge type for consensus clustering"),
                     AttributeDef("normalize_type",['row-wise', 'column-wise','both','none'],
                        "row-wise","row-wise",help="normalize type for consensus clustering"),
                    AttributeDef("create_heatmap",['no', 'yes'],
                        "no","no",help="create heatmap or not for consensus clustering"),
                     AttributeDef("cc_resample",['subsample', 'reatures','nosampling'],
                        "subsample","subsample",help="resampling scheme(one of 'subsample[ratio]','reatures[nfeat]','nosampling'"), 
                    help="Consensus Clustering File")
list_files = [ConsensusClusteringFolder]
all_modules=[
    ModuleNode(
        'consensusClustering',
         GeneExpProcessing.SignalFile,ConsensusClusteringFolder,
         OptionDef('cc_kmax',5,help="number of kmax cluster(must be >1)"),
         OptionDef('cc_resampling_iter',20,help="number of resampling iterations"),
         OptionDef('cc_seed_value','12345',help="random number generator seed"),
         OptionDef('cc_decent_iter','2000',help="Number of SOM/NMF iterations"),
         OptionDef('cc_norm_iter','0',help="Number of row/column normalization iterations (supercedes normalize.type)"),
         OptionDef('cc_heatmap_size','2',help="point size of a consensus matrix's heat map (between 1 and 20)"),
         Constraint("contents",CAN_BE_ANY_OF,BDT.CONTENTS),
         Constraint("format",MUST_BE,'gct'),
         Consequence("contents",SAME_AS_CONSTRAINT),
         help="consensus clustering")]
    

