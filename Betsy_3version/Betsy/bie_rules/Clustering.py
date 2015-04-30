#Clustering
from Betsy.bie3 import *
import GeneExpProcessing
import Database
import Heatmap
import BasicDataTypes

ClusterFile = DataType("ClusterFile", AttributeDef("cluster_alg", [
    'som', 'pca', 'kmeans', 'hierarchical'
], 'kmeans', 'kmeans',
                                                   help="cluster algorithm"),
                       AttributeDef("distance", ['correlation', 'euclidean'],
                                    'correlation', 'correlation',
                                    help="distance for cluster algorithm"),
                       AttributeDef("contents", Database.CONTENTS,
                                    'unspecified', 'unspecified',
                                    help='contents'),
                       help="Cluster file")

ClusterReportFile = DataType('ClusterReportFile',
                             help="Report file for cluster report")

list_files = [ClusterFile, ClusterReportFile]

all_modules = [
    Module('cluster_genes_by_som', GeneExpProcessing.SignalFile, ClusterFile,
           Consequence("cluster_alg", SET_TO, 'som'),
           Consequence("distance", SET_TO_ONE_OF, ['correlation', 'euclidean']),
           Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
           Consequence("contents", SAME_AS_CONSTRAINT),
           help="cluster genes by som method"),
    Module('cluster_genes_by_pca', GeneExpProcessing.SignalFile, ClusterFile,
           Consequence("cluster_alg", SET_TO, 'pca'),
           Consequence("distance", SET_TO_ONE_OF, ['correlation', 'euclidean']),
           Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
           Consequence("contents", SAME_AS_CONSTRAINT),
           help="cluster genes by pca method"),
    Module('cluster_genes_by_kmeans', GeneExpProcessing.SignalFile,
           ClusterFile, OptionDef("k_value", 5,
                                  help='k value for k means'),
           Consequence("cluster_alg", SET_TO, 'kmeans'),
           Consequence("distance", SET_TO_ONE_OF, ['correlation', 'euclidean']),
           Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
           Consequence("contents", SAME_AS_CONSTRAINT),
           help="cluster genes by kmeans method"),
    Module('cluster_genes_by_hierarchical', GeneExpProcessing.SignalFile,
           ClusterFile, Consequence("cluster_alg", SET_TO, 'hierarchical'),
           Consequence("distance", SET_TO_ONE_OF, ['correlation', 'euclidean']),
           Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
           Consequence("contents", SAME_AS_CONSTRAINT),
           help="cluster genes by hierarchical method"),
    Module('plot_cluster_heatmap', ClusterFile, Heatmap.Heatmap,
           OptionDef('hm_width', 20,
                     help="width in heatmap plot"),
           OptionDef('hm_height', 20,
                     help="heigth in heatmap plot"),
           Constraint("cluster_alg", CAN_BE_ANY_OF, ['hierarchical', 'pca',
                                                     'som', 'kmeans']),
           Constraint("distance", CAN_BE_ANY_OF, ['correlation', 'euclidean']),
           Constraint("contents", CAN_BE_ANY_OF, Database.CONTENTS),
           Consequence("cluster_alg", SAME_AS_CONSTRAINT),
           Consequence("contents", SAME_AS_CONSTRAINT),
           Consequence("distance", SAME_AS_CONSTRAINT),
           help="plot heatmap for cluster file"),
    Module(
        'make_cluster_report', [ClusterFile,
                                Heatmap.Heatmap], ClusterReportFile,
        Constraint("cluster_alg", CAN_BE_ANY_OF, ['som', 'pca', 'kmeans',
                                                  'hierarchical'], 0),
        Constraint("distance", CAN_BE_ANY_OF, ['correlation', 'euclidean'], 0),
        Constraint("cluster_alg", SAME_AS, 0, 1),
        Constraint("distance", SAME_AS, 0, 1),
        help="make cluster report"),
]
