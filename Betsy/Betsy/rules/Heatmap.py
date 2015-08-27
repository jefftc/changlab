#Heatmap

from Betsy.bie3 import *
import BasicDataTypes as BDT
import GeneExpProcessing


Heatmap = DataType(
    "Heatmap",
    AttributeDef("cluster_alg",['som','pca','kmeans','hierarchical','no_cluster_alg'],
              'no_cluster_alg','no_cluster_alg',help="cluster algorithm"),
    AttributeDef("distance",['correlation','euclidean'],'euclidean','euclidean',
                 help="distance for cluster algorithm"),
    AttributeDef('color',['red_green', 'blue_yellow'],'red_green','red_green',
                 help='color to plot the heatmap'),
    AttributeDef("contents",BDT.CONTENTS,
                               'unspecified','unspecified',help='contents'),
    help="Heatmap file")

all_data_types = [Heatmap]

all_modules = [
 ModuleNode(
        'plot_signal_heatmap',
        GeneExpProcessing.SignalFile,Heatmap,
        OptionDef('hm_width',20,help="width in heatmap plot"),
        OptionDef('hm_height',20,help="heigth in heatmap plot"),
        Constraint("format",MUST_BE,'tdf'),
        Constraint("contents",CAN_BE_ANY_OF,BDT.CONTENTS),
        Consequence("cluster_alg",SET_TO,'no_cluster_alg'),
        Consequence("distance",SET_TO_ONE_OF,['correlation','euclidean']),
        Consequence("contents",SAME_AS_CONSTRAINT),
        Consequence("color",SET_TO_ONE_OF,['red_green', 'blue_yellow']),
        help="plot heatmap for signal file"),]
