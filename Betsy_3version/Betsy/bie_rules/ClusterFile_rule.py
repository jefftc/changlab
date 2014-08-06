#ClusterFile
from Betsy.bie3 import *
import SignalFile_rule
ClusterFile = DataType(
    "ClusterFile",
    AttributeDef("cluster_alg",['som','pca','kmeans','hierarchical'],'kmeans','kmeans',
                 help="cluster algorithm"),
    AttributeDef("distance",['correlation','euclidean'],'correlation','correlation',
                 help="distance for cluster algorithm"),
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                               'unspecified','unspecified',help='contents'),
    help="Cluster file")

Heatmap = DataType(
    "Heatmap",
    AttributeDef("cluster_alg",['som','pca','kmeans','hierarchical','no_cluster_alg'],
              'no_cluster_alg','no_cluster_alg',help="cluster algorithm"),
    AttributeDef("distance",['correlation','euclidean'],'euclidean','euclidean',
                 help="distance for cluster algorithm"),
    AttributeDef('color',['red_green', 'blue_yellow'],'red_green','red_green',
                 help='color to plot the heatmap'),
    AttributeDef("contents",SignalFile_rule.CONTENTS,
                               'unspecified','unspecified',help='contents'),
    help="Heatmap file")

list_files = [ClusterFile,Heatmap]

all_modules = [
    Module(
        'cluster_genes_by_som',
        SignalFile_rule.SignalFile,ClusterFile,
        Consequence("cluster_alg",SET_TO,'som'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="cluster genes by som method"
        ),
    Module(
        'cluster_genes_by_pca',
        SignalFile_rule.SignalFile,ClusterFile,
        Consequence("cluster_alg",SET_TO,'pca'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="cluster genes by pca method"),
    Module(
        'cluster_genes_by_kmeans',
         SignalFile_rule.SignalFile,ClusterFile,
        UserInputDef("k_value",5,help='k value for k means'),
        Consequence("cluster_alg",SET_TO,'kmeans'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="cluster genes by kmeans method"),
    Module(
        'cluster_genes_by_hierarchical',
        SignalFile_rule.SignalFile,ClusterFile,
        Consequence("cluster_alg",SET_TO,'hierarchical'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("contents",SAME_AS_CONSTRAINT),
        help="cluster genes by hierarchical method"
        ),
    Module(
        'plot_signal_heatmap',
        SignalFile_rule.SignalFile,Heatmap,
        UserInputDef('hm_width',20,help="width in heatmap plot"),
        UserInputDef('hm_height',1,help="heigth in heatmap plot"),
        Constraint("format",MUST_BE,'tdf'),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("cluster_alg",SET_TO,'no_cluster_alg'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence("color",SET_TO_ONE_OF,['red_green', 'blue_yellow']),
        help="plot heatmap for signal file"),
    Module(
        'plot_cluster_heatmap',
        ClusterFile,Heatmap,
        UserInputDef('hm_width',20,help="width in heatmap plot"),
        UserInputDef('hm_height',1,help="heigth in heatmap plot"),
        Constraint("cluster_alg",CAN_BE_ANY_OF,['hierarchical','pca','som','kmeans']),
        Constraint("distance",CAN_BE_ANY_OF,['correlation','euclidean']),
        Constraint("contents",CAN_BE_ANY_OF,SignalFile_rule.CONTENTS),
        Consequence("cluster_alg",SAME_AS_CONSTRAINT),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence("distance",SAME_AS_CONSTRAINT),
        help="plot heatmap for cluster file"),
            
    ]
