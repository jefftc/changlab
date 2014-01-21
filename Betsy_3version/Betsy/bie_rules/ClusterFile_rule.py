#ClusterFile
from Betsy.bie3 import *
import SignalFile2_rule
ClusterFile = DataType(
    "ClusterFile",
    Attribute("cluster_alg",['som','pca','kmeans','hierarchical'],'kmeans','kmeans'),
    Attribute("distance",['correlation','euclidean'],'correlation','correlation'),
    Attribute('k',['yes','no'],'no','no'))


Heatmap = DataType(
    "Heatmap",
    Attribute("cluster_alg",['som','pca','kmeans','hierarchical','no_cluster_alg'],
              'no_cluster_alg','no_cluster_alg'),
    Attribute("distance",['correlation','euclidean'],'euclidean','euclidean'),
    Attribute('k',['yes','no'],'no','no'),
    Attribute('hm_width',['yes','no'],'no','no'),
    Attribute('hm_height',['yes','no'],'no','no'),
    Attribute('color',['red_green', 'blue_yellow'],'red_green','red_green'))

list_files = [ClusterFile,Heatmap]

all_modules = [
    Module(
        'cluster_genes_by_som',
        SignalFile2_rule.SignalFile2,ClusterFile,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Consequence("cluster_alg",SET_TO,'som'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])
        ),
    Module(
        'cluster_genes_by_pca',
        SignalFile2_rule.SignalFile2,ClusterFile,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Consequence("cluster_alg",SET_TO,'pca'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
    Module(
        'cluster_genes_by_kmeans',
        SignalFile2_rule.SignalFile2,ClusterFile,
        UserInput("k_value",5),
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Consequence("cluster_alg",SET_TO,'kmeans'),
        Consequence("k",SET_TO,'yes'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
    Module(
        'cluster_genes_by_hierarchical',
        SignalFile2_rule.SignalFile2,ClusterFile,
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Consequence("cluster_alg",SET_TO,'hierarchical'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
    Module(
        'plot_signal_heatmap',
        SignalFile2_rule.SignalFile2,Heatmap,
        UserInput('hm_width_value',20),
        UserInput('hm_height_value',1),
        Constraint("format",MUST_BE,"tdf"),
        Constraint("logged",MUST_BE,"yes"),
        Consequence("cluster_alg",SET_TO,'no_cluster_alg'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Consequence("hm_width",SET_TO,'yes'),
        Consequence("hm_height",SET_TO,'yes'),
        Consequence("color",SET_TO_ONE_OF,['red_green', 'blue_yellow'])),
    Module(
        'plot_signal_heatmap',
        ClusterFile,Heatmap,
        UserInput("k_value",5),
        Consequence("cluster_alg",SET_TO_ONE_OF,['hierarchical','pca','som','kmeans']),
        Consequence("hm_width",SET_TO,"yes"),
        Consequence("hm_height",SET_TO,'yes'),
        Consequence("color",SET_TO_ONE_OF,['red_green', 'blue_yellow']),
        Consequence("k",SET_TO,'yes'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean'])),
            
    ]
