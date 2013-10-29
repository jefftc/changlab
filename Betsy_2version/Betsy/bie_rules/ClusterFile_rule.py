#ClusterFile
from Betsy import bie
import SignalFile2_rule
ClusterFile = bie.DataType(
    "ClusterFile",
    # Properties of the format.
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(cluster_alg=['som','pca','kmeans','hierarchical'],DEFAULT='kmeans'),
    bie.Attribute(distance=['correlation','euclidean'],DEFAULT='euclidean'),
    bie.Attribute(k=bie.ANYATOM,DEFAULT='5'))


Heatmap = bie.DataType(
    "Heatmap",
    bie.Attribute(filename=bie.ANYATOM, DEFAULT="", OPTIONAL=True),
    bie.Attribute(cluster_alg=['som','pca','kmeans','hierarchical','no_cluster_alg'],
              DEFAULT='no_cluster_alg'),
    bie.Attribute(distance=['correlation','euclidean'],DEFAULT='euclidean'),
    bie.Attribute(k=bie.ANYATOM,DEFAULT='5'),
    bie.Attribute(hm_width=bie.ANYATOM,DEFAULT="20"),
    bie.Attribute(hm_height=bie.ANYATOM,DEFAULT='1'),
    bie.Attribute(color=['red_green', 'blue_yellow'],DEFAULT='red_green'))

list_files = [ClusterFile,Heatmap]

all_modules = [
    bie.Module(
        'cluster_genes_by_som',
        SignalFile2_rule.SignalFile2(format="tdf", logged="yes"),
        ClusterFile(cluster_alg='som',
                    k=bie.ANYATOM,distance=['correlation','euclidean'])),
    bie.Module(
        'cluster_genes_by_pca',
        SignalFile2_rule.SignalFile2(format="tdf", logged="yes"),
        ClusterFile(cluster_alg='pca',
                    k=bie.ANYATOM,distance=['correlation','euclidean'])),
    bie.Module(
        'cluster_genes_by_kmeans',
        SignalFile2_rule.SignalFile2(format="tdf", logged="yes"),
        ClusterFile(cluster_alg='kmeans',
                    k=bie.ANYATOM,distance=['correlation','euclidean'])),
    bie.Module(
        'cluster_genes_by_hierarchical',
        SignalFile2_rule.SignalFile2(format="tdf", logged="yes"),
        ClusterFile(cluster_alg='hierarchical',
                    k=bie.ANYATOM,distance=['correlation','euclidean'])),
    bie.Module(
        'plot_signal_heatmap',
        SignalFile2_rule.SignalFile2(format="tdf", logged="yes"),
        Heatmap(cluster_alg='no_cluster_alg',
                    hm_width=bie.ANYATOM, hm_height=bie.ANYATOM,
                    color=['red_green', 'blue_yellow'])),
    bie.Module(
        'plot_signal_heatmap',
        ClusterFile,
        Heatmap(cluster_alg=['hierarchical','pca','som','kmeans'],
                hm_width=bie.ANYATOM, hm_height=bie.ANYATOM,
                color=['red_green', 'blue_yellow'],
                k=bie.ANYATOM,distance=['correlation','euclidean']))
    ]
