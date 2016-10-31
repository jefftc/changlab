from Betsy.bie3 import *
import GeneExpProcessing as GXP
import BasicDataTypes as BDT

ClusterFolder = DataType(
    "ClusterFolder",
    AttributeDef(
        "cluster_alg", ["hierarchical", "kmeans", "som"],
        "kmeans", "kmeans",
        help="cluster algorithm"),
    #AttributeDef(
    #    "distance", ["correlation", "euclidean"], "correlation", "correlation",
    #    help="distance for cluster algorithm"),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    help="Folder of files that contain the result of a cluster analysis.  "
    "Files formats are defined by the TreeView program.  "
    "Possible files are: 1) signal.cdt.  Gene expression values.  "
    "2) array_tree.atr  3) gene_tree.gtr  "
    "4) array_cluster.kag  5) gene_cluster.kgg")


#ClusterReportFile = DataType(
#    "ClusterReportFile",
#    # Not sure what this is.
#    help="Report file for cluster report",
#    )


Heatmap = DataType(
    "Heatmap",
    AttributeDef(
        "cluster_alg", ["none", "som", "kmeans", "hierarchical"],
        "none", "none",
        help="Clustering algorithm."),
    AttributeDef(
        "contents", BDT.CONTENTS, "unspecified", "unspecified",
        help="contents"),
    help="A PNG file containing a heatmap.",
    )


all_data_types = [
    ClusterFolder,
    #ClusterReportFile,
    Heatmap,
    ]


all_modules = [
    ModuleNode(
        "cluster_genes_by_hierarchical",
        GXP.SignalFile, ClusterFolder,
        OptionDef(
            "cluster_genes", default="yes",
            help="Whether to cluster genes.  Options: yes, no"),
        OptionDef(
            "cluster_arrays", default="yes",
            help="Whether to cluster arrays.  Options: yes, no"),
        OptionDef(
            "distance_measure", default="pearson",
            help="Which distance metric to use.  "
            "Options: uncent-cor, pearson, abs-uncent-cor, "
            "abs-pearson, spearman, kendall, euclidean, city-block."),
        OptionDef(
            "linkage", default="complete",
            help="How to link branches.  "
            "Options: complete, single, centroid, average."),
        Consequence("cluster_alg", SET_TO, "hierarchical"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="cluster genes by hierarchical method"
        ),
    ModuleNode(
        "cluster_genes_by_kmeans",
        GXP.SignalFile, ClusterFolder,
        OptionDef(
            "cluster_genes", default="yes",
            help="Whether to cluster genes.  Options: yes, no"),
        OptionDef(
            "cluster_arrays", default="yes",
            help="Whether to cluster arrays.  Options: yes, no"),
        OptionDef(
            "distance_measure", default="pearson",
            help="Which distance metric to use.  "
            "Options: uncent-cor, pearson, abs-uncent-cor, "
            "abs-pearson, spearman, kendall, euclidean, city-block."),
        OptionDef("kmeans_k", 5, help="Number of clusters."),
        
        Consequence("cluster_alg", SET_TO, "kmeans"),
        #Consequence("distance", SET_TO_ONE_OF, ["correlation", "euclidean"]),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="cluster genes by kmeans method"),
    # XXX
    ModuleNode(
        "cluster_genes_by_som",
        GXP.SignalFile, ClusterFolder,
        OptionDef(
            "cluster_genes", default="yes",
            help="Whether to cluster genes.  Options: yes, no"),
        OptionDef(
            "cluster_arrays", default="yes",
            help="Whether to cluster arrays.  Options: yes, no"),
        OptionDef(
            "distance_measure", default="pearson",
            help="Which distance metric to use.  "
            "Options: uncent-cor, pearson, abs-uncent-cor, "
            "abs-pearson, spearman, kendall, euclidean, city-block."),
        OptionDef(
            "som_rows", default="2",
            help="How many rows for the SOM grid."),
        OptionDef(
            "som_cols", default="2",
            help="How many columns for the SOM grid."),
        
        Consequence("cluster_alg", SET_TO, "som"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="cluster genes by som method"
        ),
    #ModuleNode(
    #    "cluster_genes_by_pca",
    #    GXP.SignalFile, ClusterFolder,
    #    Consequence("cluster_alg", SET_TO, "pca"),
    #    #Consequence("distance", SET_TO_ONE_OF, ["correlation", "euclidean"]),
    #    Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    #    Consequence("contents", SAME_AS_CONSTRAINT),
    #    help="cluster genes by pca method"),
   
    #ModuleNode(
    #    "make_cluster_report",
    #    [ClusterFolder, Heatmap], ClusterReportFile,
    #    Constraint(
    #        "cluster_alg", CAN_BE_ANY_OF,
    #        ["som", "pca", "kmeans", "hierarchical"], 0),
    #    #Constraint("distance", CAN_BE_ANY_OF, ["correlation", "euclidean"], 0),
    #    Constraint("cluster_alg", SAME_AS, 0, 1),
    #    #Constraint("distance", SAME_AS, 0, 1),
    #    help="make cluster report"),

    ModuleNode(
        "plot_signal_heatmap",
        GXP.SignalFile, Heatmap,
        OptionDef(
            "hm_width", default="20",
            help="Width (in pixels) of each column in the heatmap."),
        OptionDef(
            "hm_height", default="20", 
            help="Height (in pixels) of each row in the heatmap."),
        OptionDef(
            "hm_color", default="brewer-rdylbu-div",
            help="Distance metric for clustering.  "
            "Options: red, white, red-green, blue-yellow, red-green-soft, "
            "red-blue-soft, matlab, bild, genepattern, genespring, yahoo, "
            "brewer-prgn-div, brewer-rdbu-div, brewer-rdylbu-div, "
            "brewer-rdylgn-div, brewer-spectral-div, brewer-blues-seq, "
            "brewer-greens-seq, brewer-reds-seq, brewer-ylorbr-seq, "
            "or brewer-qual-set1"),

        OptionDef(
            "hm_colorbar", default="no",
            help="Whether to show a colorbar with the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_colorbar_horizontal", default="no",
            help="Draw a horizontal colorbar.  Options: yes, no."),
        OptionDef(
            "hm_colorbar_height", default="1.0",
            help="Scale the height of the colorbar by this amount."),
        OptionDef(
            "hm_colorbar_width", default="1.0",
            help="Scale the width of the colorbar by this amount."),
        OptionDef(
            "hm_colorbar_font", default="1.0",
            help="Scale the font of the colorbar by this amount."),

        OptionDef(
            "hm_label_genes", default="no",
            help="Whether to label the genes in the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_scale_gene_labels", default="1.0",
            help="Change the size of the gene labels.  "
            "1.0 means default size, 1.5 will draw 50% bigger."),
        OptionDef(
            "hm_label_arrays", default="no",
            help="Whether to label the arrays in the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_scale_array_labels", default="1.0",
            help="Change the size of the array labels.  "
            "1.0 means default size, 1.5 will draw 50% bigger."),
        #OptionDef(
        #    "cluster_alg", default="none",
        #    help="How to cluster the data.  "
        #    "Options: none (default), som, pca, kmeans, hierarchical."),
        #OptionDef(
        #    "distance", default="correlation",
        #    help="Distance metric for clustering.  "
        #    "Options: correlation (default), euclidean."),
        Constraint("format", MUST_BE, "tdf"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("cluster_alg", MUST_BE, "none"),
        Consequence("cluster_alg", SET_TO, "none"),
        #Consequence("distance", SET_TO_ONE_OF, ["correlation", "euclidean"]),
        #Consequence("color", SET_TO_ONE_OF, ["red_green", "blue_yellow"]),
        help="Plot heatmap for signal file",
        ),

    ModuleNode(
        "plot_cluster_heatmap",
        ClusterFolder, Heatmap,
        OptionDef(
            "hm_width", default="20",
            help="Width (in pixels) of each column in the heatmap."),
        OptionDef(
            "hm_height", default="20", 
            help="Height (in pixels) of each row in the heatmap."),
        OptionDef(
            "hm_color", default="brewer-rdylbu-div",
            help="Distance metric for clustering.  "
            "Options: red, white, red-green, blue-yellow, red-green-soft, "
            "red-blue-soft, matlab, bild, genepattern, genespring, yahoo, "
            "brewer-prgn-div, brewer-rdbu-div, brewer-rdylbu-div, "
            "brewer-rdylgn-div, brewer-spectral-div, brewer-blues-seq, "
            "brewer-greens-seq, brewer-reds-seq, brewer-ylorbr-seq, "
            "or brewer-qual-set1"),

        OptionDef(
            "hm_colorbar", default="no",
            help="Whether to show a colorbar with the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_colorbar_horizontal", default="no",
            help="Draw a horizontal colorbar.  Options: yes, no."),
        OptionDef(
            "hm_colorbar_height", default="1.0",
            help="Scale the height of the colorbar by this amount."),
        OptionDef(
            "hm_colorbar_width", default="1.0",
            help="Scale the width of the colorbar by this amount."),
        OptionDef(
            "hm_colorbar_font", default="1.0",
            help="Scale the font of the colorbar by this amount."),
        
        OptionDef(
            "hm_label_genes", default="no",
            help="Whether to label the genes in the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_scale_gene_labels", default="1.0",
            help="Change the size of the gene labels.  "
            "1.0 means default size, 1.5 will draw 50% bigger."),
        OptionDef(
            "hm_label_arrays", default="no",
            help="Whether to label the arrays in the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_scale_array_labels", default="1.0",
            help="Change the size of the array labels.  "
            "1.0 means default size, 1.5 will draw 50% bigger."),

        OptionDef(
            "hm_show_gene_tree", default="yes",
            help="Whether to show the dendrogram for the genes.  "
            "Options: yes, no"),
        OptionDef(
            "hm_show_array_tree", default="yes",
            help="Whether to show the dendrogram for the genes.  "
            "Options: yes, no"),
        OptionDef(
            "hm_show_gene_cluster", default="yes",
            help="Whether to show the dendrogram for the genes.  "
            "Options: yes, no"),
        OptionDef(
            "hm_show_array_cluster", default="yes",
            help="Whether to show the dendrogram for the genes.  "
            "Options: yes, no"),
        
        Constraint(
            "cluster_alg", CAN_BE_ANY_OF, ["hierarchical", "som", "kmeans"]),
        Consequence("cluster_alg", SAME_AS_CONSTRAINT),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        #Constraint("distance", CAN_BE_ANY_OF, ["correlation", "euclidean"]),
        #Consequence("distance", SAME_AS_CONSTRAINT),
        help="Plot a heatmap from a clustered signal file.",
        ),
    ]
