# DataTypes:
# SignalDistributionBoxplot
# ActbPlot
# HousekeepingPlot
# ControlPlot
# IlluHybridizationProbePlot
# BiotinPlot
# PCAPlot
# Heatmap
#
# ClusterFolder
#
# Modules:
# plot_intensity_boxplot
# plot_actb_line
# plot_illu_housekeeping_line
# plot_affy_affx_line
# plot_illu_hyb_bar
# plot_illu_biotin_line
# plot_sample_pca
#
# plot_signal_heatmap
# plot_cluster_heatmap
# 
# cluster_genes_by_kmeans
# cluster_genes_by_hierarchical


from Betsy.bie3 import *
import SignalFile
import BasicDataTypes as BDT
YESNO = BDT.YESNO


SignalDistributionBoxplot = DataType(
    'SignalDistributionBoxplot',
    SignalFile.ATTR_CONTENTS,
    SignalFile.ATTR_PREPROCESS,
    SignalFile.ATTR_LOGGED,
    AttributeDef(
        "relabel_sample", YESNO, "no", "no", help="rename sample or not"),
    help="Boxplot of distribution of signal files."
    )

ActbPlot = DataType(
    "ActbPlot",
    SignalFile.ATTR_CONTENTS,
    SignalFile.ATTR_PREPROCESS,
    SignalFile.ATTR_LOGGED,
    AttributeDef(
        "relabel_sample", YESNO, "no", "no", help="rename sample or not"),
    help="Gene expression of beta-actin.")

HousekeepingPlot = DataType(
    'HousekeepingPlot',
    SignalFile.ATTR_CONTENTS,
    AttributeDef(
        "relabel_sample", YESNO, "no", "no", help="rename sample or not"),
    help="Plot of Illumina's Housekeeping control probes.")

ControlPlot = DataType(
    "ControlPlot",
    *(SignalFile.UNPROC_ATTRIBUTES),
    help="control plot file.  WHAT IS THIS?")

#Hyb_barPlot = DataType(
IlluHybridizationProbePlot = DataType(
    'IlluHybridizationProbePlot',
    SignalFile.ATTR_CONTENTS,
    help="Expression of high, medium, and low hybridization controls "
    "for Illumina microarrays.")

BiotinPlot = DataType(
    'BiotinPlot',
    SignalFile.ATTR_CONTENTS,
    AttributeDef(
        "relabel_sample", YESNO, "no", "no", help="rename sample or not"),
    help="Biotin plot file")

PCAPlot = DataType(
    "PCAPlot",
    # Technically, don't need:
    # annotate
    # platform
    # Keep anyway for simplicity.
    # 
    # relabel_sample is helpful if we're labeling the points.
    *(SignalFile.UNPROC_ATTRIBUTES +
      SignalFile.IMPUTE_ATTRIBUTES +
      SignalFile.MERGE_ATTRIBUTES +
      SignalFile.NORMALIZE_ATTRIBUTES +
      SignalFile.ORDER_ATTRIBUTES +
      SignalFile.ANNOTATE_ATTRIBUTES +
      SignalFile.FILTER_ATTRIBUTES),
    help="PNG file with a PCA plot."
    )

Heatmap = DataType(
    "Heatmap",
    # Technically, don't need:
    # annotate
    # platform
    # Keep anyway for simplicity.
    # 
    # relabel_sample is helpful if we're labeling the points.
    AttributeDef(
        "cluster_alg", ["none", "som", "kmeans", "hierarchical"],
        "none", "hierarchical",
        help="Clustering algorithm."),
    *(SignalFile.UNPROC_ATTRIBUTES +
      SignalFile.IMPUTE_ATTRIBUTES +
      SignalFile.MERGE_ATTRIBUTES +
      SignalFile.NORMALIZE_ATTRIBUTES +
      SignalFile.ORDER_ATTRIBUTES +
      SignalFile.ANNOTATE_ATTRIBUTES +
      SignalFile.FILTER_ATTRIBUTES),
    help="A PNG file containing a heatmap."
    )

ClusterFolder = DataType(
    "ClusterFolder",
    AttributeDef(
        "cluster_alg", ["hierarchical", "kmeans", "som"],
        "kmeans", "kmeans",
        help="cluster algorithm"),
    *(SignalFile.UNPROC_ATTRIBUTES +
      SignalFile.IMPUTE_ATTRIBUTES +
      SignalFile.MERGE_ATTRIBUTES +
      SignalFile.NORMALIZE_ATTRIBUTES +
      SignalFile.ORDER_ATTRIBUTES +
      SignalFile.ANNOTATE_ATTRIBUTES +
      SignalFile.FILTER_ATTRIBUTES),
    help="Folder of files that contain the result of a cluster analysis.  "
    "Files formats are defined by the TreeView program.  "
    "Possible files are: 1) signal.cdt.  Gene expression values.  "
    "2) array_tree.atr  3) gene_tree.gtr  "
    "4) array_cluster.kag  5) gene_cluster.kgg")

all_data_types = [
    SignalDistributionBoxplot,
    ActbPlot,
    HousekeepingPlot,
    ControlPlot,
    IlluHybridizationProbePlot,
    BiotinPlot,
    PCAPlot,
    Heatmap,
    ClusterFolder,
    ]

all_modules = [
    ModuleNode(
        'plot_intensity_boxplot',
        SignalFile.SignalFile, SignalDistributionBoxplot,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO),
        Consequence("relabel_sample", SAME_AS_CONSTRAINT),
        help="plot intensity boxplot"
        ),
    ModuleNode(
        'plot_actb_line',
        #SignalFile._SignalFile_Impute, ActbPlot,
        #Constraint("missing_values", MUST_BE, "no"),
        SignalFile.SignalFile, ActbPlot,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("logged", MUST_BE, "yes"),
        Consequence("logged", SAME_AS_CONSTRAINT),
        Constraint("annotate", MUST_BE, "yes"),
        Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO),
        Consequence("relabel_sample", SAME_AS_CONSTRAINT),
        help="Plot the expression values of ACTB."),
    ModuleNode(
        'plot_illu_housekeeping_line',
        SignalFile.IlluminaControlFile, HousekeepingPlot,
        #Constraint("preprocess", MUST_BE, 'illumina'),
        #Constraint("format", MUST_BE, "gct"),
        #Constraint("logged", MUST_BE, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO),
        Consequence("relabel_sample", SAME_AS_CONSTRAINT),
        help="plot illumina housekeeping line"),
    ModuleNode(
        'plot_affy_affx_line',
        SignalFile.SignalFile, ControlPlot,
        Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
        Consequence("preprocess", SAME_AS_CONSTRAINT),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot affy affx line"),
    ModuleNode(
        "plot_illu_hyb_bar",
        SignalFile.IlluminaControlFile, IlluHybridizationProbePlot,
        #Constraint("preprocess", MUST_BE, "illumina"),
        #Constraint("format", MUST_BE, "gct"),
        #Constraint("logged", MUST_BE, "no"),
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        help="plot illumina affx line"),
    ModuleNode(
        'plot_illu_biotin_line',
        SignalFile.IlluminaControlFile, BiotinPlot,
        Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
        Consequence("contents", SAME_AS_CONSTRAINT),
        Constraint("relabel_sample", CAN_BE_ANY_OF, YESNO),
        Consequence("relabel_sample", SAME_AS_CONSTRAINT),
        help="plot illumina biotin line"),

    ModuleNode(
        "plot_sample_pca",
        SignalFile.SignalFile, PCAPlot,
        #OptionDef(
        #    "pca_num_genes", default="500",
        #    help="Number of genes to use when calculating "
        #    "principal components.",
        #    ),
        *(
            SignalFile.CONVERT_UNPROC +
            SignalFile.CONVERT_IMPUTE +
            SignalFile.CONVERT_MERGE +
            SignalFile.CONVERT_NORMALIZE +
            SignalFile.CONVERT_ORDER +
            SignalFile.CONVERT_ANNOTATE +
            SignalFile.CONVERT_FILTER
            ),
        help="Make a PCA plot."
        ),
##     ModuleNode(
##         "plot_signal_heatmap",
##         SignalFile.SignalFile, Heatmap,
##         OptionDef(
##             "hm_width", default="20",
##             help="Width (in pixels) of each column in the heatmap."),
##         OptionDef(
##             "hm_height", default="20", 
##             help="Height (in pixels) of each row in the heatmap."),
##         OptionDef(
##             "hm_color", default="brewer-rdylbu-div",
##             help="Distance metric for clustering.  "
##             "Options: red, white, red-green, blue-yellow, red-green-soft, "
##             "red-blue-soft, matlab, bild, genepattern, genespring, yahoo, "
##             "brewer-prgn-div, brewer-rdbu-div, brewer-rdylbu-div, "
##             "brewer-rdylgn-div, brewer-spectral-div, brewer-blues-seq, "
##             "brewer-greens-seq, brewer-reds-seq, brewer-ylorbr-seq, "
##             "or brewer-qual-set1"),

##         OptionDef(
##             "hm_colorbar", default="no",
##             help="Whether to show a colorbar with the heatmap.  "
##             "Options: yes, no"),
##         OptionDef(
##             "hm_colorbar_horizontal", default="no",
##             help="Draw a horizontal colorbar.  Options: yes, no."),
##         OptionDef(
##             "hm_colorbar_height", default="1.0",
##             help="Scale the height of the colorbar by this amount."),
##         OptionDef(
##             "hm_colorbar_width", default="1.0",
##             help="Scale the width of the colorbar by this amount."),
##         OptionDef(
##             "hm_colorbar_font", default="1.0",
##             help="Scale the font of the colorbar by this amount."),

##         OptionDef(
##             "hm_label_genes", default="no",
##             help="Whether to label the genes in the heatmap.  "
##             "Options: yes, no"),
##         OptionDef(
##             "hm_scale_gene_labels", default="1.0",
##             help="Change the size of the gene labels.  "
##             "1.0 means default size, 1.5 will draw 50% bigger."),
##         OptionDef(
##             "hm_label_arrays", default="no",
##             help="Whether to label the arrays in the heatmap.  "
##             "Options: yes, no"),
##         OptionDef(
##             "hm_scale_array_labels", default="1.0",
##             help="Change the size of the array labels.  "
##             "1.0 means default size, 1.5 will draw 50% bigger."),
##         #OptionDef(
##         #    "cluster_alg", default="none",
##         #    help="How to cluster the data.  "
##         #    "Options: none (default), som, pca, kmeans, hierarchical."),
##         #OptionDef(
##         #    "distance", default="correlation",
##         #    help="Distance metric for clustering.  "
##         #    "Options: correlation (default), euclidean."),
##         Constraint("format", MUST_BE, "tdf"),
##         Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
##         Consequence("contents", SAME_AS_CONSTRAINT),
##         Constraint("preprocess", CAN_BE_ANY_OF, BDT.PREPROCESS),
##         Consequence("preprocess", SAME_AS_CONSTRAINT),
##         #Constraint("cluster_alg", MUST_BE, "none"),
##         Consequence("cluster_alg", SET_TO, "none"),
##         Constraint("filter_and_threshold", CAN_BE_ANY_OF, YESNO),
##         Consequence("filter_and_threshold", SAME_AS_CONSTRAINT),
##         help="Plot heatmap for signal file",
##         ),
    ModuleNode(
        "plot_cluster_heatmap",
        ClusterFolder, Heatmap,
        OptionDef(
            "hm_width", default="",
            help="Width (in pixels) of each column in the heatmap."),
        OptionDef(
            "hm_height", default="", 
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
            "hm_label_genes", default="",
            help="Whether to label the genes in the heatmap.  "
            "Options: yes, no"),
        OptionDef(
            "hm_scale_gene_labels", default="1.0",
            help="Change the size of the gene labels.  "
            "1.0 means default size, 1.5 will draw 50% bigger."),
        OptionDef(
            "hm_label_arrays", default="",
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
            "cluster_alg", CAN_BE_ANY_OF, ["som", "kmeans", "hierarchical"]),
        Consequence("cluster_alg", SAME_AS_CONSTRAINT),
        *(
            SignalFile.CONVERT_UNPROC +
            SignalFile.CONVERT_IMPUTE +
            SignalFile.CONVERT_MERGE +
            SignalFile.CONVERT_NORMALIZE +
            SignalFile.CONVERT_ORDER +
            SignalFile.CONVERT_ANNOTATE +
            SignalFile.CONVERT_FILTER
            ),
        help="Plot a heatmap from a clustered signal file."
        ),

    ModuleNode(
        "cluster_genes_by_hierarchical",
        SignalFile.SignalFile, ClusterFolder,
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
        *(
            SignalFile.CONVERT_UNPROC +
            SignalFile.CONVERT_IMPUTE +
            SignalFile.CONVERT_MERGE +
            SignalFile.CONVERT_NORMALIZE +
            SignalFile.CONVERT_ORDER +
            SignalFile.CONVERT_ANNOTATE +
            SignalFile.CONVERT_FILTER
            ),
        help="cluster genes by hierarchical method"
        ),
    ModuleNode(
        "cluster_genes_by_kmeans",
        SignalFile.SignalFile, ClusterFolder,
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
        *(
            SignalFile.CONVERT_UNPROC +
            SignalFile.CONVERT_IMPUTE +
            SignalFile.CONVERT_MERGE +
            SignalFile.CONVERT_NORMALIZE +
            SignalFile.CONVERT_ORDER +
            SignalFile.CONVERT_ANNOTATE +
            SignalFile.CONVERT_FILTER
            ),
        help="cluster genes by kmeans method"),
    #ModuleNode(
    #    "cluster_genes_by_som",
    #    SignalFile.SignalFile, ClusterFolder,
    #    OptionDef(
    #        "cluster_genes", default="yes",
    #        help="Whether to cluster genes.  Options: yes, no"),
    #    OptionDef(
    #        "cluster_arrays", default="yes",
    #        help="Whether to cluster arrays.  Options: yes, no"),
    #    OptionDef(
    #        "distance_measure", default="pearson",
    #        help="Which distance metric to use.  "
    #        "Options: uncent-cor, pearson, abs-uncent-cor, "
    #        "abs-pearson, spearman, kendall, euclidean, city-block."),
    #    OptionDef(
    #        "som_rows", default="2",
    #        help="How many rows for the SOM grid."),
    #    OptionDef(
    #        "som_cols", default="2",
    #        help="How many columns for the SOM grid."),
    #    
    #    Consequence("cluster_alg", SET_TO, "som"),
    #    Constraint("contents", CAN_BE_ANY_OF, BDT.CONTENTS),
    #    Consequence("contents", SAME_AS_CONSTRAINT),
    #    help="cluster genes by som method"
    #    ),
    ]
