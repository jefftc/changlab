#ClusterFile
from Betsy.bie3 import *
import SignalFile_rule
ClusterFile = DataType(
    "ClusterFile",
    AttributeDef("cluster_alg",['som','pca','kmeans','hierarchical'],'kmeans','kmeans'),
    AttributeDef("distance",['correlation','euclidean'],'correlation','correlation'),
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified'))

Heatmap = DataType(
    "Heatmap",
    AttributeDef("cluster_alg",['som','pca','kmeans','hierarchical','no_cluster_alg'],
              'no_cluster_alg','no_cluster_alg'),
    AttributeDef("distance",['correlation','euclidean'],'euclidean','euclidean'),
    AttributeDef('hm_width',['yes','no'],'yes','yes'),
    AttributeDef('hm_height',['yes','no'],'yes','yes'),
    AttributeDef('color',['red_green', 'blue_yellow'],'red_green','red_green'),
    AttributeDef("contents",["train0","train1", "test",
                             "class0,class1,test","class0",
                              "class1", "class0,class1","unspecified"],
                               'unspecified','unspecified'))

list_files = [ClusterFile,Heatmap]

all_modules = [
    Module(
        'cluster_genes_by_som',
        SignalFile_rule.SignalFile,ClusterFile,
        Consequence("cluster_alg",SET_TO,'som'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        ),
    Module(
        'cluster_genes_by_pca',
        SignalFile_rule.SignalFile,ClusterFile,
        Consequence("cluster_alg",SET_TO,'pca'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
    Module(
        'cluster_genes_by_kmeans',
         SignalFile_rule.SignalFile,ClusterFile,
        UserInputDef("k_value",5),
        Consequence("cluster_alg",SET_TO,'kmeans'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
    Module(
        'cluster_genes_by_hierarchical',
        SignalFile_rule.SignalFile,ClusterFile,
        Consequence("cluster_alg",SET_TO,'hierarchical'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
    Module(
        'plot_signal_heatmap',
        SignalFile_rule.SignalFile,Heatmap,
        UserInputDef('hm_width_value',20),
        UserInputDef('hm_height_value',1),
        Constraint("format",MUST_BE,'tdf'),
        Consequence("cluster_alg",SET_TO,'no_cluster_alg'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Consequence("hm_width",SET_TO,'yes'),
        Consequence("hm_height",SET_TO,'yes'),
        Consequence("color",SET_TO_ONE_OF,['red_green', 'blue_yellow'])),
    Module(
        'plot_signal_heatmap',
        ClusterFile,Heatmap,
        UserInputDef("k_value",5),
        Constraint("cluster_alg",CAN_BE_ANY_OF,['hierarchical','pca','som','kmeans']),
        Constraint("distance",CAN_BE_ANY_OF,['correlation','euclidean']),
        Consequence("cluster_alg",SAME_AS_CONSTRAINT),
        Consequence("distance",SAME_AS_CONSTRAINT)),
            
    ]
