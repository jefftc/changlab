#!/usr/bin/env python

# To do:
# - Allow choice of gene label.
# - Make transparent background.
# - Add borders between elements in the plot.
# - What if there are nan's in the probability file?

"""

Classes:
ClusterData

PlotLayout              Holds layouts for each individual pieces of the plot.
HeatmapLayout
DendrogramLayout
GeneDendrogramLayout
ArrayDendrogramLayout
GeneClusterLayout
ArrayClusterLayout
GeneLabelLayout
ArrayLabelLayout


Functions:
process_data_set
make_layout
calc_coords_for_layout

convert_to_pcl
filter_matrix
normalize_matrix
cluster_matrix
pretty_scale_matrix
find_data_files
read_data_set
write_data_set

plot
plot_matrix
plot_dendrogram
plot_gene_clusters
plot_array_clusters
plot_gene_labels
plot_array_labels

_cluster
_cleanup_cluster
_get_gene_ids        Return unique IDs for the genes.
_get_array_ids       Return unique IDs for the arrays.
_get_gene_labels     Return pretty labels for the genes.
_get_array_labels    Return pretty labels for the arrays.
_parse_gene_names
_get_color
_exists_nz

"""

import os, sys

## Detect if I'm being run from within GenePattern.  If so, then add
## the current directory to the library path.
#if os.path.split(__file__)[0].endswith("gp_pybinreg"):
#    sys.path.append(os.path.split(__file__)[0])

MIN_FONTSIZE = 6
MAX_MEGAPIXELS = 256  # No more than 256 megapixel plot.

class ClusterData:
    def __init__(
        self, gene_tree, array_tree, gene_tree_cluster, array_tree_cluster,
        gene_cluster, array_cluster):
        self.gene_tree = gene_tree
        self.array_tree = array_tree
        self.gene_tree_cluster = gene_tree_cluster
        self.array_tree_cluster = array_tree_cluster
        self.gene_cluster = gene_cluster
        self.array_cluster = array_cluster

class PlotLayout:
    def __init__(self, heatmap, gene_dendrogram, array_dendrogram,
                 gene_cluster, array_cluster, gene_label, array_label):
        self.heatmap = heatmap
        self.gene_dendrogram = gene_dendrogram
        self.array_dendrogram = array_dendrogram
        self.gene_cluster = gene_cluster
        self.array_cluster = array_cluster
        self.gene_label = gene_label
        self.array_label = array_label

class PlotCoords:
    def __init__(self, hm_x, hm_y, gd_x, gd_y, ad_x, ad_y,
                 gc_x, gc_y, ac_x, ac_y, gl_x, gl_y, al_x, al_y):
        self.hm_x, self.hm_y = hm_x, hm_y
        self.gd_x, self.gd_y = gd_x, gd_y
        self.ad_x, self.ad_y = ad_x, ad_y
        self.gc_x, self.gc_y = gc_x, gc_y
        self.ac_x, self.ac_y = ac_x, ac_y
        self.gl_x, self.gl_y = gl_x, gl_y
        self.al_x, self.al_y = al_x, al_y

class HeatmapLayout:
    def __init__(self, nrow, ncol, boxwidth, boxheight, color_fn):
        self.nrow = nrow
        self.ncol = ncol
        self.boxwidth = boxwidth
        self.boxheight = boxheight
        self.color_fn = color_fn
    def width(self):
        return self.size()[0]
    def height(self):
        return self.size()[1]
    def size(self):
        height = self.boxheight * self.nrow
        width = self.boxwidth * self.ncol
        return width, height
    def coord(self, row, col):
        x = col * self.boxwidth
        y = row * self.boxheight
        return x, y, self.boxwidth, self.boxheight
    def color(self, x):
        # x is from [0, 1].  find the nearest color.
        assert x >= 0 and x <= 1, "x out of range: %g" % x
        return _get_color(x, self.color_fn)

class DendrogramLayout:
    def __init__(
        self, num_items, num_other_items,
        pixels_per_item, pixels_per_other_item,
        size_scale, thickness_scale, tree, tree_cluster, color_fn):
        # This dendrogram is measured in 2 dimensions: the dimension
        # that spans across the branches, and the dimension that spans
        # across the phylogenetic distance.  color_fn is for the
        # clusters.
        import math
        
        self.num_items = num_items
        self.pixels_per_item = pixels_per_item
        self.tree = tree
        self.tree_cluster = tree_cluster
        self.color_fn = color_fn
        self.max_cluster = None
        if self.tree_cluster:
            self.max_cluster = max(self.tree_cluster.values())
        self._item_size = num_items * pixels_per_item
        
        # Should be the same height as the heatmap.  The width should
        # be 0.625x the height (perfect ratio).
        RATIO = 1.0 / 1.6 / 2.0
        # Both dendrograms should have symmetric sizes.  So base the
        # RATIO on the smaller of the two dimensions.
        x1 = num_items * pixels_per_item
        x2 = num_other_items * pixels_per_other_item
        x = min(x1, x2)
        self._dist_size = int(math.ceil(x * RATIO * size_scale))

        # Convert the distances from clustering to percentages.  The
        # percentages indicate how far across the plot to place the
        # node.  0% is the furthest distance, while 100% is the
        # closest.
        
        # These are actually similarity metrics, so put 1.0 at 100%,
        # and the lowest value given in the tree (can be negative) at
        # 0%.
        lowest = highest = None
        for node in tree:
            left, right, distance = node
            if lowest is None or distance < lowest:
                lowest = distance
            if highest is None or distance > highest:
                highest = distance
        assert highest <= 1.0
        
        # Set the closest to always be 1.0.
        highest = 1.0
        # Set a small border at the end for the root.
        self.ROOT_SIZE = 0.15 * (highest - lowest)
        lowest -= self.ROOT_SIZE
        self.lowest, self.highest = lowest, highest

        min_ppi = min(pixels_per_item, pixels_per_other_item)
        x = int(math.ceil(min_ppi*0.20 * thickness_scale))
        x = min(max(x, 1), min_ppi)
        self.LINEWIDTH = x
    def vthicken(self, x, y, width, height):
        import math
        np = self.LINEWIDTH - width
        if np <= 0:
            return x, y, width, height
        hnp = int(math.floor(np/2.0))
        return x-hnp, y, width+np, height
    def hthicken(self, x, y, width, height):
        import math
        np = self.LINEWIDTH - height
        if np <= 0:
            return x, y, width, height
        hnp = int(math.floor(np/2.0))
        return x, y-hnp, width, height+np
    def item_size(self):
        return self._item_size
    def dist_size(self):
        return self._dist_size
    def color(self, id):
        c = 0, 0, 0
        n = None
        if self.tree_cluster:
            n = self.tree_cluster[id]
        if n is not None and self.max_cluster:
            p = float(n) / self.max_cluster
            c = _get_color(p, self.color_fn)
        return c
    def item_coord(self, item):
        x = int(item * self.pixels_per_item + self.pixels_per_item/2.0)
        return x
    def dist_coord(self, distance):
        assert distance >= self.lowest and distance <= self.highest
        # Convert the distance to a percentage.
        perc = (distance - self.lowest) / (self.highest - self.lowest)
        x = int(perc * self.dist_size())
        return x
        
class GeneDendrogramLayout(DendrogramLayout):
    def __init__(self, nrow, ncol, boxwidth, boxheight,
                 size_scale, thickness_scale, tree, tree_cluster, color_fn):
        DendrogramLayout.__init__(
            self, nrow, ncol, boxheight, boxwidth, size_scale, thickness_scale,
            tree, tree_cluster, color_fn)
    def width(self):
        return self.size()[0]
    def height(self):
        return self.size()[1]
    def size(self):
        height = self.item_size()
        width = self.dist_size()
        return width, height
    def coord(self, row, distance):
        x = self.dist_coord(distance)
        y = self.item_coord(row)
        return x, y
    def lines(self, node_num, node_dist, left_num, left_dist,
              right_num, right_dist):
        import math
        node_x, node_y = self.coord(node_num, node_dist)
        left_x, left_y = self.coord(left_num, left_dist)
        right_x, right_y = self.coord(right_num, right_dist)

        # The right node is on top of the left node.
        lines = []
        # Add two left lines, then two right lines.
        x = node_x, left_y, left_x-node_x+1, 1
        lines.append(self.vthicken(node_x, node_y, 1, left_y-node_y+1))
        lines.append(self.hthicken(node_x, left_y, left_x-node_x+1, 1))
        lines.append(self.vthicken(node_x, right_y, 1, node_y-right_y+1))
        lines.append(self.hthicken(node_x, right_y, right_x-node_x+1, 1))
        return lines
    def root(self, node_num, node_dist):
        root_x, root_y = self.coord(node_num, self.lowest)
        node_x, node_y = self.coord(node_num, node_dist)
        x = root_x, root_y, node_x-root_x+1, 1
        return self.hthicken(*x)

class ArrayDendrogramLayout(DendrogramLayout):
    # Bug: This doesn't work if the boxheight is changed (e.g. -y 2).
    def __init__(self, nrow, ncol, boxwidth, boxheight,
                 size_scale, thickness_scale, tree, tree_cluster, color_fn):
        DendrogramLayout.__init__(
            self, ncol, nrow, boxwidth, boxheight, size_scale, thickness_scale,
            tree, tree_cluster, color_fn)
    def width(self):
        return self.size()[0]
    def height(self):
        return self.size()[1]
    def size(self):
        height = self.dist_size()
        width = self.item_size()
        return width, height
    def coord(self, row, distance):
        x = self.item_coord(row)
        y = self.dist_coord(distance)
        return x, y
    def lines(self, node_num, node_dist, left_num, left_dist,
              right_num, right_dist):
        node_x, node_y = self.coord(node_num, node_dist)
        left_x, left_y = self.coord(left_num, left_dist)
        right_x, right_y = self.coord(right_num, right_dist)

        lines = []
        # Add two left lines, then two right lines.
        lines.append(self.hthicken(left_x, node_y, right_x-left_x+1, 1))
        lines.append(self.vthicken(left_x, node_y, 1, left_y-node_y+1))
        lines.append(self.hthicken(node_x, node_y, right_x-node_x+1, 1))
        lines.append(self.vthicken(right_x, node_y, 1, right_y-node_y+1))
        return lines
    def root(self, node_num, node_dist):
        root_x, root_y = self.coord(node_num, self.lowest)
        node_x, node_y = self.coord(node_num, node_dist)
        x = root_x, root_y, 1, node_y-root_y+1
        return self.vthicken(*x)

class GeneClusterLayout:
    def __init__(self, num_items, item_width, item_height):
        self.num_items = num_items
        self.item_width = item_width
        self.item_height = item_height
    def width(self):
        return self.item_width
    def height(self):
        return self.item_height * self.num_items
    def size(self):
        return self.width(), self.height()
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self.num_items
        x = 0
        y = num * self.item_height
        return x, y, self.item_width, self.item_height

class ArrayClusterLayout:
    def __init__(self, num_items, item_width, item_height):
        self.num_items = num_items
        self.item_width = item_width
        self.item_height = item_height
    def width(self):
        return self.item_width * self.num_items
    def height(self):
        return self.item_height
    def size(self):
        return self.width(), self.height()
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self.num_items
        x = num * self.item_width
        y = 0
        return x, y, self.item_width, self.item_height

class GeneLabelLayout:
    def __init__(self, item_height, item_widths):
        self._item_height = item_height
        self._item_widths = item_widths
        num_items = len(item_widths)
        self._width = max(item_widths)
        self._height = item_height * num_items
        self._num_items = num_items
    def item_height(self):
        return self._item_height
    def width(self):
        return self._width
    def height(self):
        return self._height
    def size(self):
        return self.width(), self.height()
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self._num_items
        x = 0
        y = num * self._item_height
        return x, y, self._item_widths[num], self._item_height

class ArrayLabelLayout:
    def __init__(self, item_height, item_widths):
        # item_height refers to the text, not rotated.
        self._item_height = item_height
        self._item_widths = item_widths
        num_items = len(item_widths)
        # _width and _height refer to the layout object.
        self._width = item_height * num_items
        self._height = max(item_widths)
        self._num_items = num_items
    def item_height(self):
        return self._item_height
    def width(self):
        return self._width
    def height(self):
        return self._height
    def size(self):
        return self.width(), self.height()
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self._num_items
        x = num * self._item_height
        y = self._height-self._item_widths[num]
        return x, y, self._item_height, self._item_widths[num]

def process_data_set(
    MATRIX, cluster, cluster_data, jobname,
    gene_indexes, gene_names, gene_file, num_genes_var,
    array_indexes, array_file,
    log_transform, gene_center, gene_normalize, array_center, array_normalize,
    cluster_genes, cluster_arrays, cluster_alg, distance, method,
    gene_k, array_k, kmeans_k, scale, gain, autoscale):
    import arrayio

    assert MATRIX.nrow() > 0, "Matrix has no genes."
    MATRIX = filter_matrix(
        MATRIX, gene_indexes, gene_names, gene_file, num_genes_var,
        array_indexes, array_file)
    assert MATRIX.nrow() > 0, "Filtered matrix has no genes."
    MATRIX, cluster_data = normalize_matrix(
        MATRIX, cluster, cluster_data, log_transform,
        gene_center, gene_normalize, array_center, array_normalize)
    MATRIX, cluster_data = cluster_matrix(
        MATRIX, cluster, cluster_data, cluster_genes, cluster_arrays,
        cluster_alg, distance, method, gene_k, array_k, kmeans_k)
    if jobname:
        write_data_set(MATRIX, cluster_data, jobname)
    # Scale after the clustering, so it doesn't affect the clustering
    # results.
    MATRIX = pretty_scale_matrix(MATRIX, scale, gain, autoscale)
        
    return MATRIX, cluster_data

def make_layout(
    MATRIX, cluster_data, plotlib, boxwidth, boxheight, color_scheme, 
    cluster_genes, gene_tree_scale, gene_tree_thickness,
    cluster_arrays, array_tree_scale, array_tree_thickness,
    cluster_alg, label_genes, label_arrays):
    from genomicode import colorfns

    # Choose the color scheme.
    scheme2fn = {
        "red" : colorfns.red_shade,
        "red-green" : colorfns.rg_array_colors,
        "blue-yellow" : colorfns.by_array_colors,
        "matlab" : colorfns.matlab_colors,
        "bild" : colorfns.bild_colors,
        "genespring" : colorfns.genespring_colors,
        "yahoo" : colorfns.yahoo_weather_colors,
        }
    assert color_scheme in scheme2fn, "Unknown color scheme: %s" % color_scheme
    color_fn = scheme2fn[color_scheme]

    # Make the layout for the heatmap.
    hm_layout = HeatmapLayout(
        MATRIX.nrow(), MATRIX.ncol(), boxwidth, boxheight, color_fn)

    # Make layouts for the dendrograms.
    gd_layout = ad_layout = None
    if(cluster_genes and cluster_data.gene_tree and gene_tree_scale > 0 and 
       cluster_alg == "hierarchical" and MATRIX.ncol() > 1):
        # Only add the dendrogram if hierarchical clustering was
        # requested.  If clustering not done, then the matrix file
        # will not have the GID annotations, and there will be no way
        # to match up the genes with the clusters.
        assert gene_tree_scale > 0
        assert gene_tree_thickness > 0
        gd_layout = GeneDendrogramLayout(
            MATRIX.nrow(), MATRIX.ncol(), boxwidth, boxheight,
            gene_tree_scale, gene_tree_thickness, 
            cluster_data.gene_tree, cluster_data.gene_tree_cluster,
            colorfns.matlab_colors)
    if(cluster_arrays and cluster_data.array_tree and array_tree_scale > 0 and
       cluster_alg == "hierarchical" and MATRIX.nrow() > 1):
        assert array_tree_scale > 0
        assert array_tree_thickness > 0
        ad_layout = ArrayDendrogramLayout(
            MATRIX.nrow(), MATRIX.ncol(), boxwidth, boxheight,
            array_tree_scale, array_tree_thickness, 
            cluster_data.array_tree, cluster_data.array_tree_cluster, 
            colorfns.matlab_colors)

    # Make layouts for the clusters.
    # Can plot these (k-means) clusters if either kmeans or
    # hierarchical clustering was requested.  Unlike hierarchical
    # clustering, plotting this does not require any annotations in
    # the matrix file.
    gc_layout = ac_layout = None
    if cluster_data.gene_cluster:
        gc_layout = GeneClusterLayout(MATRIX.nrow(), boxwidth, boxheight)
    if cluster_data.array_cluster:
        ac_layout = ArrayClusterLayout(MATRIX.ncol(), boxwidth, boxheight)

    # Make the layout for the gene or array labels.
    gl_layout = al_layout = None
    gene_labels = array_labels = None
    gl_fontsize = al_fontsize = None
    if label_genes:
        gene_labels = _get_gene_labels(MATRIX)
        gl_fontsize = plotlib.fit_fontsize_to_height(boxheight)
        widths = [plotlib.get_text_size(x, gl_fontsize)[0]
                  for x in gene_labels]
        gl_layout = GeneLabelLayout(boxheight, widths)
    if label_arrays:
        array_labels = _get_array_labels(MATRIX)
        al_fontsize = plotlib.fit_fontsize_to_height(boxwidth)
        widths = [plotlib.get_text_size(x, al_fontsize)[0]
                  for x in array_labels]
        al_layout = ArrayLabelLayout(boxwidth, widths)

    x = PlotLayout(
        hm_layout, gd_layout, ad_layout, gc_layout, ac_layout,
        gl_layout, al_layout)
    return x

def calc_coords_for_layout(layout):
    x = y = 0
    def _safe_size(layout):
        if layout is None:
            return 0, 0
        return layout.size()

    hm_width, hm_height = _safe_size(layout.heatmap)
    gd_width, gd_height = _safe_size(layout.gene_dendrogram)
    ad_width, ad_height = _safe_size(layout.array_dendrogram)
    gc_width, gc_height = _safe_size(layout.gene_cluster)
    ac_width, ac_height = _safe_size(layout.array_cluster)
    gl_width, gl_height = _safe_size(layout.gene_label)
    al_width, al_height = _safe_size(layout.array_label)
    
    # On X-axis: gene dendrogram, cluster, label, then heatmap.
    gd_x, gd_y = x, y+ad_height+ac_height+al_height
    gc_x, gc_y = gd_x+gd_width, gd_y
    gl_x, gl_y = gc_x+gc_width, gd_y
    
    # On Y-axis: array dendrogram, cluster, label, then heatmap.
    ad_x, ad_y = x+gd_width+gc_width+gl_width, y
    ac_x, ac_y = ad_x, ad_y+ad_height
    al_x, al_y = ad_x, ac_y+ac_height

    # Put heatmap at the end.
    hm_x, hm_y = ad_x, gd_y

    x = PlotCoords(
        hm_x, hm_y, gd_x, gd_y, ad_x, ad_y,
        gc_x, gc_y, ac_x, ac_y, gl_x, gl_y, al_x, al_y)
    return x

def _choose_gene_id(MATRIX):
    # Given a user-specified matrix, try to pick a good unique ID for
    # the genes.
    import arrayio
    
    headers = MATRIX.row_names()

    # Prioritize some potential ones.  Don't use the standard headers,
    # e.g. arrayio.ROW_ID, so that we can preserve the user's header.
    IDS = ["Probe.Set.ID"]
    for id in IDS:
        if id in headers:
            return id

    # If no known headers are found, then choose a standard one.
    IDS = [arrayio.AFFY_PROBESET_ID, arrayio.GENE_ID, arrayio.ROW_ID]
    for id in IDS:
        if id in headers:
            return id

    # If no standard ones are found, then just arbitrarily use the
    # first column.
    if headers:
        return headers[0]
    
    raise AssertionError, "I could not find an ID for the matrix."

def _choose_gene_label(MATRIX):
    import arrayio

    names = MATRIX.row_names()

    # Prioritize some potential ones.
    IDS = [
        "DESCRIPTION",   # For GCT files.
        "NAME",
        arrayio.GENE_SYMBOL, "Gene.Symbol",
        arrayio.GENE_DESCRIPTION, "Description",
        arrayio.GENE_ID, "LocusLink",
        arrayio.AFFY_PROBESET_ID, "Probe.Set.ID",
        arrayio.ROW_ID
        ]
    for id in IDS:
        if id in names:
            return id
    if names:
        return names[0]
    raise AssertionError, "I could not find an ID for the matrix."

def convert_to_pcl(MATRIX, label_name=None):
    # Convert the matrix to PCL format.
    # Row names   <ID>  NAME
    # Col names   
    import arrayio
    from genomicode import Matrix

    # Select from the row names an ID and a NAME.
    id_name = _choose_gene_id(MATRIX)
    name_name = _choose_gene_label(MATRIX)

    # Should not use "GID" as column name for PCL file.  When
    # clustering, cluster will add another "GID" column, and then
    # there will be two columns called "GID".  Rename this to
    # something else, if necessary.
    pretty_id_name = id_name
    if pretty_id_name == "GID":
        pretty_id_name = "GID.OLD"
    if pretty_id_name == "NAME":
        # GCT files uses "NAME" for ID, which conflicts with PCL definition.
        pretty_id_name = "ID.NAME"
    pretty_name_name = "NAME"

    SAMPLE_NAME = arrayio.tab_delimited_format.SAMPLE_NAME 
    row_order = [pretty_id_name, pretty_name_name]
    col_order = [SAMPLE_NAME]
    row_names = {}
    col_names = {}
    synonyms = {}
    
    row_names[pretty_id_name] = MATRIX.row_names(id_name)
    row_names[pretty_name_name] = MATRIX.row_names(name_name)
    col_names[SAMPLE_NAME] = MATRIX.col_names(arrayio.COL_ID)
    synonyms[arrayio.ROW_ID] = pretty_id_name
    synonyms[arrayio.COL_ID] = SAMPLE_NAME

    x = Matrix.InMemoryMatrix(
        MATRIX.slice(), row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order)
    pcl_matrix = Matrix.add_synonyms(x, synonyms)
    assert arrayio.pcl_format.is_matrix(pcl_matrix)
    return pcl_matrix

def read_filecol(filecol):
    from genomicode import iofns

    # filecol is either <filename> or <filename>,<col>.  commas
    # are not allowed in the filenames.
    filename, colnum = filecol, 0
    if filecol.find(",") >= 0:
        x = filecol.split(",")
        assert len(x) == 2, "File should be specified: <filename>,<col>"
        filename, colnum = x
        colnum = int(colnum)
    assert os.path.exists(filename), "could not find file %s" % filename
    data = iofns.split_tdf(open(filename).read())
    names = [x[colnum].strip() for x in data]
    names = [x for x in names if x]
    return names

def _parse_gene_names(gene_name_list):
    # This can a list of comma separated genes, e.g.
    # ["E2F1", "E2F2,E2F3"]
    # Need to separate them out.
    gene_names = []
    for x in gene_name_list:
        x = x.split(",")
        gene_names.extend(x)
    return gene_names

def filter_matrix(
    MATRIX, gene_indexes, gene_names, gene_filecol, num_genes_var,
    array_indexes, array_filecol):
    # Filter the genes, maintaining the order specified in the input.
    from genomicode import pcafns
    from genomicode import parsefns
    
    # User provides indexes as 1-based inclusive.  Convert to 0-based
    # exclusive.
    if gene_indexes is not None:
        I = []
        for s, e in parsefns.parse_ranges(gene_indexes):
            assert s >= 1
            s, e = s-1, e
            if e > MATRIX.nrow():
                e = MATRIX.nrow()
            I.extend(range(s, e))
        MATRIX = MATRIX.matrix(I, None)

    # Use the names specified by the user.
    gene_names = _parse_gene_names(gene_names)

    # Read gene names from a file.
    if gene_filecol:
        x = read_filecol(gene_filecol)
        gene_names.extend(x)

    if array_indexes is not None:
        I = []
        for s, e in parsefns.parse_ranges(array_indexes):
            assert s >= 1
            s, e = s-1, e
            if e > MATRIX.ncol():
                e = MATRIX.ncol()
            I.extend(range(s, e))
        MATRIX = MATRIX.matrix(None, I)

    # Read array names from a file.
    array_names = []
    if array_filecol:
        x = read_filecol(array_filecol)
        array_names.extend(x)
    #print "GENES", gene_names
    #print "ARRAYS", array_names
    #print MATRIX.row_names()

    # Specify genes by the annotations.  Microarray formats do not
    # name genes, but give them annotations.  Specify arrays by their
    # name.
    if gene_names and array_names:
        MATRIX = MATRIX.matrix(row=gene_names, col_name=array_names)
    elif gene_names:
        MATRIX = MATRIX.matrix(row=gene_names)
    elif array_names:
        MATRIX = MATRIX.matrix(col_name=array_names)

    # Now select the genes based on variance.
    if num_genes_var and MATRIX.nrow() > num_genes_var:
        I = pcafns.select_genes_var(MATRIX._X, num_genes_var)
        MATRIX = MATRIX.matrix(I, None)
        
    return MATRIX

def normalize_matrix(
    MATRIX, cluster, cluster_data, log_transform, gene_center, gene_normalize,
    array_center, array_normalize):
    # log_transform    boolean
    # gene_center      None, "mean", or "median"
    # gene_normalize   None, "ss", or "var"
    # array_center     None, "mean", or "median"
    # array_normalize  None, "ss", or "var"

    import arrayio

    # If no normalization requested, then just return the matrix.
    if (not log_transform and
        not gene_center and not gene_normalize and
        not array_center and not array_normalize):
        return MATRIX, cluster_data

    args = []
    if log_transform:
        args.append("-l")
    if gene_center == "mean":
        args.append("-cg a")
    elif gene_center == "median":
        args.append("-cg m")
    if gene_normalize == "ss":
        args.append("-ng")
    if array_center == "mean":
        args.append("-ca a")
    elif array_center == "median":
        args.append("-ca m")
    if array_normalize == "ss":
        args.append("-na")

    # No clustering.
    args.append("-g 0")
    args.append("-e 0")

    filestem = _cluster(MATRIX, cluster=cluster, *args)
    files = find_data_files(filestem)
    assert "nrm" in files, "No normalization file produced."
    MATRIX, cluster_data = read_data_set(filestem, cluster_data)
    _cleanup_cluster(filestem)

    if gene_normalize == "var":
        normalize_genes_var(MATRIX)
    if array_normalize == "var":
        normalize_arrays_var(MATRIX)

    return MATRIX, cluster_data

def normalize_genes_var(MATRIX):
    from genomicode import jmath
    
    # Normalize the genes in place.
    X = MATRIX._X
    for i in range(len(X)):
        X_i = X[i]
        m = jmath.mean(X_i)
        # Subtract the mean.
        X_i = [x-m for x in X_i]
        # Normalize to stddev of 1.
        s = jmath.stddev(X_i)
        if s != 0:
            X_i = [x/s for x in X_i]
        # Add the mean back.
        X_i = [x+m for x in X_i]
        X[i] = X_i

def normalize_arrays_var(MATRIX):
    from genomicode import jmath

    # Normalize the arrays in place.
    X = MATRIX._X
    if not X or not X[0]:
        return
    for i in range(len(X[0])):
        X_i = [x[i] for x in X]
        m = jmath.mean(X_i)
        # Subtract the mean.
        X_i = [x-m for x in X_i]
        # Normalize to stddev of 1.
        s = jmath.stddev(X_i)
        if s != 0:
            X_i = [x/s for x in X_i]
        # Add the mean back.
        X_i = [x+m for x in X_i]
        for j in range(len(X)):
            X[j][i] = X_i[j]

def cluster_matrix(
    MATRIX, cluster, cluster_data,
    cluster_genes, cluster_arrays, algorithm, distance, method,
    gene_k, array_k, kmeans_k):
    import arrayio
    from genomicode import clusterio

    assert algorithm in ["hierarchical", "kmeans"]

    dist2id = {
        "uncent-cor"  : 1,  "pearson"    : 2,  "abs-uncent-cor" : 3,
        "abs-pearson" : 4,  "spearman"   : 5,  "kendall"        : 6,
        "euclidean"   : 7,  "city-block" : 8,
        }
    method2id = {
        "complete" : "m", "single" : "s", "centroid" : "c", "average" : "a",
        }

    # Skip if all conditions are true:
    # - not clustering genes
    # - not clustering arrays
    # - not cutting gene tree (and gene tree already exists)
    # - not cutting array tree (and array tree already exists)
    if (not cluster_genes and not cluster_arrays and
        not (gene_k and cluster_data.gene_tree) and 
        not (array_k and cluster_data.array_tree)):
        return MATRIX, cluster_data

    # If not clustering and just re-cutting the tree, then don't
    # bother regenerating the clusters.
    if cluster_genes or cluster_arrays:
        args = []

        id = dist2id[distance]
        if cluster_genes:
            args.append("-g %s" % id)
        else:
            args.append("-g 0")
        if cluster_arrays:
            args.append("-e %s" % id)
        else:
            args.append("-e 0")

        id = method2id[method]
        args.append("-m %s" % id)

        if algorithm == "kmeans":
            args.append("-k %d" % kmeans_k)

        filestem = _cluster(MATRIX, cluster=cluster, *args)
        files = find_data_files(filestem)
        assert "cdt" in files, "No cdt file produced."
        MATRIX, cluster_data = read_data_set(filestem, cluster_data)
        _cleanup_cluster(filestem)
    
    # Cluster the hierarchical trees, if necessary.
    gene_tree_cluster = array_tree_cluster = None
    # If I haven't reclustered the data, then the old tree is still
    # valid.
    if not cluster_genes:
        gene_tree_cluster = cluster_data.gene_tree_cluster
    if not cluster_arrays:
        array_tree_cluster = cluster_data.array_tree_cluster
    if cluster_data.gene_tree and gene_k:
        assert gene_k <= MATRIX.nrow(), "more gene clusters than genes"
        gene_tree_cluster = clusterio.cut_dendrogram(
            cluster_data.gene_tree, gene_k)
    if cluster_data.array_tree and array_k:
        assert array_k <= MATRIX.ncol(), "more array clusters than arrays"
        array_tree_cluster = clusterio.cut_dendrogram(
            cluster_data.array_tree, array_k)
    cluster_data.gene_tree_cluster = gene_tree_cluster
    cluster_data.array_tree_cluster = array_tree_cluster

    return MATRIX, cluster_data

def pretty_scale_matrix(MATRIX, scale, gain, autoscale):
    # Find a good default gain value.  Values should be scaled from
    # [-1, 1].
    # Will change the MATRIX variable.

    nrow, ncol = MATRIX.dim()
    X = MATRIX._X

    # Apply the scale specified by the user.
    for i in range(nrow):
        for j in range(ncol):
            X[i][j] = X[i][j] + scale

    # Choose a default gain based on the maximum expression level.
    defgain = 1.0
    if autoscale:
        x_max = None
        for i in range(nrow):
            for j in range(ncol):
                if x_max is None or abs(X[i][j]) > x_max:
                    x_max = abs(X[i][j])
        defgain = 1.0/x_max
        defgain = defgain * 2.0    # By default, automatically multiply by 2.0.

    # Apply the gain specified by the user.
    for i in range(nrow):
        for j in range(ncol):
            # The gain is scaled from the default gain.
            X[i][j] = X[i][j] * defgain * gain

    #for x in X:
    #   print "\t".join(map(str, x))

    # Finally, rescale to [0, 1].
    for i in range(nrow):
        for j in range(ncol):
            x = X[i][j]
            x = (x + 1) * 0.5
            x = max(min(x, 1), 0)
            X[i][j] = x

    return MATRIX

def _guess_filestem(file_or_job):
    # Examples:
    # file_or_job          stem
    # test                 test
    # GSE5451.l2.mas5.gtr  GSE5451.l2.mas5
    # GSE1456.mas5.gz      GSE1456
    # GSE1456.mas5         GSE1456
    # out.dat              out
    # out.pcl              out
    # out.txt              out
    # 
    # Rule:
    # - If the file doesn't exist, then use the whole thing as the
    #   stem.
    # - If there's a .gz, then chop it off.
    # - If there's an extension, then chop it off.
    
    if not _exists_nz(file_or_job):
        return file_or_job
    
    stem = file_or_job

    # Chop off the .gz at the end.
    COMPRESSION_EXTS = [".gz", ".bz2", ".zip"]
    s, e = os.path.splitext(stem)
    for ext in COMPRESSION_EXTS:
        if e.lower() == ext:
            stem = s
            break

    # Chop off one more extension, if it exists.
    stem, e = os.path.splitext(stem)

    return stem

def find_data_files(file_or_stem):
    # Return a dictionary of extension -> filename.

    # Needs to take a stem, because it can be hard to find the file
    # names based on the CDT file, because the stem is different.
    # Example files:
    # <stem>.gtc
    # <stem>_K_A5.kag
    # <stem>_K_G5.kgg
    # <stem>_K_G5_A5.cdt    Hard to cut into a stem.
    fullstem = _guess_filestem(file_or_stem)

    # Bug: should do a case insensitive search.
    path, stem = os.path.split(fullstem)
    if not path:
        path = "."
    EXTENSIONS = [
        "nrm", "pcl", "cdt", "gtr", "atr", "gtc", "atc", "kgg", "kag"]
    ext2file = {}
    for file in os.listdir(path):
        if not file.startswith(stem):
            continue
        f, e = os.path.splitext(file)
        if e.startswith("."):
            e = e[1:]
        if e not in EXTENSIONS:
            continue
        
        # Do some checking to make sure file looks reasonable.
        recognize_file = False
        if f == stem:
            recognize_file = True
        elif f.startswith("%s_K_A" % stem):
            recognize_file = True
        elif f.startswith("%s_K_G" % stem):
            recognize_file = True
        if not recognize_file:
            continue
        
        ext2file[e] = os.path.join(path, file)
    return ext2file
    
def read_data_set(file_or_stem, default=None):
    import arrayio
    from genomicode import parsefns
    from genomicode import Matrix
    from genomicode import clusterio

    files = find_data_files(file_or_stem)

    filename = file_or_stem
    if not os.path.exists(filename):
        # If this file does not exist, then look for a CDT, NRM, or
        # PCL file (in that order).
        DATA_EXTS = ["cdt", "nrm", "pcl"]
        DATA_EXTS = [x for x in DATA_EXTS if x in files]
        assert DATA_EXTS, "I could not find the expression data file."
        ext = DATA_EXTS[0]
        filename = files[ext]
    MATRIX = arrayio.read(filename)

    # If no gene IDs were provided, then just make some up.
    if not MATRIX.row_names():
        header = "GENE.ID"
        MATRIX._row_order.append(header)
        x = ["R%s" % x for x in parsefns.pretty_range(0, MATRIX.nrow())]
        MATRIX._row_names[header] = x
        synonyms[arrayio.ROW_ID] = header
        MATRIX = Matrix.add_synonyms(MATRIX, synonyms)

    # Read the clustering files.
    formats = [
        ("gtr", clusterio.read_gtr_file),
        ("atr", clusterio.read_atr_file),
        ("gtc", clusterio.read_gtc_file),
        ("atc", clusterio.read_atc_file),
        ("kgg", clusterio.read_kgg_file),
        ("kag", clusterio.read_kag_file),
        ]
    data = {}  # ext -> output
    for ext, read_fn in formats:
        if ext not in files:
            continue
        data[ext] = read_fn(files[ext])

    if default is None:
        default = ClusterData(None, None, None, None, None, None)

    cluster_data = ClusterData(
        data.get("gtr", default.gene_tree),
        data.get("atr", default.array_tree),
        data.get("gtc", default.gene_tree_cluster),
        data.get("atc", default.array_tree_cluster),
        data.get("kgg", default.gene_cluster),
        data.get("kag", default.array_cluster),
        )
    return MATRIX, cluster_data

def write_data_set(MATRIX, cluster_data, jobname):
    from arrayio import tab_delimited_format
    from genomicode import clusterio

    matrix_file = "%s.cdt" % jobname
    tab_delimited_format.write(MATRIX, open(matrix_file, 'w'))

    cd = cluster_data
    formats = [
        ("gtr", clusterio.write_gtr_file, cd.gene_tree),
        ("atr", clusterio.write_atr_file, cd.array_tree),
        ("gtc", clusterio.write_gtc_file, cd.gene_tree_cluster),
        ("atc", clusterio.write_atc_file, cd.array_tree_cluster),
        ("kgg", clusterio.write_kgg_file, cd.gene_cluster),
        ("kag", clusterio.write_kag_file, cd.array_cluster),
        ]

    for ext, write_fn, data in formats:
        if not data:
            continue
        outfile = "%s.%s" % (jobname, ext)
        write_fn(data, open(outfile, 'w'))

def plot(filename, MATRIX, cluster_data, plotlib, layout, coords):
    plot_width = coords.hm_x + layout.heatmap.width()
    plot_height = coords.hm_y + layout.heatmap.height()

    # Plot each element of the figure.
    image = plotlib.image(plot_width, plot_height)
    if layout.gene_dendrogram:
        plot_dendrogram(
            plotlib, image, MATRIX, coords.gd_x, coords.gd_y,
            layout.gene_dendrogram, "GENE", cluster_data.gene_tree)
    if layout.array_dendrogram:
        plot_dendrogram(
            plotlib, image, MATRIX, coords.ad_x, coords.ad_y,
            layout.array_dendrogram, "ARRAY", cluster_data.array_tree)
    if layout.gene_cluster:
        plot_gene_clusters(
            plotlib, image, MATRIX, coords.gc_x, coords.gc_y,
            layout.gene_cluster, cluster_data.gene_cluster)
    if layout.array_cluster:
        plot_array_clusters(
            plotlib, image, MATRIX, coords.ac_x, coords.ac_y,
            layout.array_cluster, cluster_data.array_cluster)
    if layout.gene_label:
        gene_labels = _get_gene_labels(MATRIX)
        plot_gene_labels(
            plotlib, image, MATRIX, coords.gl_x, coords.gl_y,
            layout.gene_label, gene_labels)
    if layout.array_label:
        array_labels = _get_array_labels(MATRIX)
        #print "HERE2", plot_width, plot_height, coords.al_x, coords.al_y
        plot_array_labels(
            plotlib, image, MATRIX, coords.al_x, coords.al_y,
            layout.array_label, array_labels)
    plot_matrix(
        plotlib, image, MATRIX, coords.hm_x, coords.hm_y, layout.heatmap)

    plotlib.write(image, open(filename, 'w'))
    
def plot_matrix(plotlib, image, MATRIX, xoff, yoff, layout):
    # Draw the actual matrix.
    X = MATRIX._X
    for i in range(MATRIX.nrow()):
        for j in range(MATRIX.ncol()):
            x = X[i][j]
            c = layout.color(x)

            # Find the coordinates and plot it.
            x, y, width, height = layout.coord(i, j)
            plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)

def plot_dendrogram(plotlib, image, MATRIX, xoff, yoff, layout, dim, tree):
    from genomicode import clusterio

    if dim == "GENE":
        ids = MATRIX.row_names("GID")
    elif dim == "ARRAY":
        ids = MATRIX.col_names("AID")
    else:
        raise AssertionError, "Unknown dim: %s" % dim

    # num is the row or column of the node.
    id2num = {}  # gene or node id -> num
    id2distance = {}  # gene or node id -> distance
    # Find id2num and id2distance for each of the leaves.
    for i, x in enumerate(ids):
        id = clusterio.parse_node(x)
        id2num[id] = i
        id2distance[id] = 1

    # Find id2num and id2distance for each of the internal nodes.
    for i, node in enumerate(tree):
        id = -(i+1)
        left, right, distance = node
        left_num = id2num[left]
        right_num = id2num[right]
        id2num[id] = (left_num + right_num)/2.0
        id2distance[id] = distance

    # Draw each of the nodes of the tree.
    for i, node in enumerate(tree):
        node_id = -(i+1)
        left_id, right_id, node_dist = node
        node_num = id2num[node_id]
        left_num = id2num[left_id]
        right_num = id2num[right_id]
        left_dist = id2distance[left_id]
        right_dist = id2distance[right_id]

        # Two left lines, then two right lines.
        lines = layout.lines(
            node_num, node_dist, left_num, left_dist, right_num, right_dist)
        ll1, ll2, lr1, lr2 = lines

        # Line 1 is the top joining line, and line 2 is the line that
        # leads to the node.
        cl1 = cl2 = layout.color(left_id)
        cr1 = cr2 = layout.color(right_id)
        # Color the node and everything beneath it.  Don't color the
        # lines on top of the node.  The exception is if there's only
        # a single leaf in the cluster, then color the bottom-most
        if cl1 != cr1:
            cl1 = cr1 = (0, 0, 0)
            if left_id < 0:
                cl2 = (0, 0, 0)
            if right_id < 0:
                cr2 = (0, 0, 0)

        data = [(ll1, cl1), (ll2, cl2), (lr1, cr1), (lr2, cr2)]
        for line, color in data:
            x, y, width, height = line
            plotlib.rectangle(image, x+xoff, y+yoff, width, height, color)

    c = layout.color(node_id)
    x, y, width, height = layout.root(node_num, node_dist)
    plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)

def plot_gene_clusters(plotlib, image, X, xoff, yoff, layout, clusters):
    import arrayio
    from genomicode import colorfns
    assert X.nrow() == len(clusters), "%d %d" % (X.nrow(), len(clusters))

    # Figure out what kind of IDs to use.
    ID_NAMES = ["GID", "NAME", arrayio.ROW_ID]
    ID_NAMES = [x for x in ID_NAMES if x in X.row_names()]
    ids = [x[0] for x in clusters]
    for ID_NAME in ID_NAMES:
        ID = X.row_names(ID_NAME)
        num_found = 0
        for id in ids:
            if id in ID:
                num_found += 1
        if num_found == len(ids):
            break
    else:
        raise AssertionError, "I could not find the cluster IDs."
    
    GID = X.row_names(ID_NAME)
    gid2i = {}
    for i, gid in enumerate(GID):
        gid2i[gid] = i
        
    max_cluster = max([x[1] for x in clusters])
    for gid, n in clusters:
        i = gid2i[gid]
        x, y, width, height = layout.coord(i)
        p = 0.5
        if max_cluster > 0:
            p = float(n) / max_cluster
        # Bug: This should be set in the layout.
        c = _get_color(p, colorfns.matlab_colors)
        plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)

def plot_array_clusters(plotlib, image, X, xoff, yoff, layout, clusters):
    import arrayio
    from genomicode import colorfns
    assert X.ncol() == len(clusters)

    # Figure out what kind of IDs to use.
    ID_NAMES = [
        "AID", arrayio.COL_ID, arrayio.tab_delimited_format.SAMPLE_NAME]
    ID_NAMES = [x for x in ID_NAMES if x in X.col_names()]
    ids = [x[0] for x in clusters]
    for ID_NAME in ID_NAMES:
        ID = X.col_names(ID_NAME)
        num_found = 0
        for id in ids:
            if id in ID:
                num_found += 1
        if num_found == len(ids):
            break
    else:
        raise AssertionError, "I could not find the array IDs."

    AID = X.col_names(ID_NAME)
    aid2i = {}
    for i, aid in enumerate(AID):
        aid2i[aid] = i
        
    max_cluster = max([x[1] for x in clusters])
    for aid, n in clusters:
        i = aid2i[aid]
        x, y, width, height = layout.coord(i)
        p = 0.5
        if max_cluster > 0:
            p = float(n) / max_cluster
        # Bug: This should be set in the layout.
        c = _get_color(p, colorfns.matlab_colors)
        plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)

def plot_gene_labels(plotlib, image, X, xoff, yoff, layout, labels):
    fontsize = plotlib.fit_fontsize_to_height(layout.item_height())
    if fontsize < MIN_FONTSIZE:
        return
    for i in range(X.nrow()):
        x, y, width, height = layout.coord(i)
        plotlib.text(image, xoff+x, yoff+y, labels[i], fontsize, (0, 0, 0))

def plot_array_labels(plotlib, image, X, xoff, yoff, layout, labels):
    fontsize = plotlib.fit_fontsize_to_height(layout.item_height())
    if fontsize < MIN_FONTSIZE:
        return
    for i in range(X.ncol()):
        x, y, width, height = layout.coord(i)
        #print "HERE1", xoff, yoff, x, y, width, height, fontsize
        plotlib.text90(image, xoff+x, yoff+y, labels[i], fontsize, (0, 0, 0))

def _cluster(MATRIX, *args, **params):
    import tempfile
    import subprocess
    import arrayio
    from genomicode import filefns

    cluster = params.get("cluster") or "cluster"

    path = "."
    x, filestem = tempfile.mkstemp(dir=path); os.close(x)
    filefns.safe_unlink(filestem)

    # Write the data set in PCL format.
    # This implementation requires a matrix in PCL format.
    MATRIX = convert_to_pcl(MATRIX)
    pcl_file = filestem + ".pcl"
    arrayio.pcl_format.write(MATRIX, open(pcl_file, 'w'))

    args = list(args)
    args.append('-f "%s"' % pcl_file)

    cmd = "%s %s" % (cluster, " ".join(args))
    #print cmd
    #w, r = os.popen4(cmd)
    p = subprocess.Popen(
        cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout
    w.close()
    output = r.read()
    if output.find("cluster: command not found") >= 0:
        raise AssertionError, "cluster: command not found"
    elif output.find("command not found") >= 0:
        raise AssertionError, "%s: command not found" % cluster

    return filestem

def _cleanup_cluster(filestem):
    # Just remove all the files with the filestem.
    from genomicode import filefns

    path, filestem = os.path.split(filestem)
    for file in os.listdir(path):
        if not file.startswith(filestem):
            continue
        filename = os.path.join(path, file)
        filefns.safe_unlink(filename)

def _get_gene_ids(MATRIX):
    names = MATRIX.row_names()

    KNOWN_NAMES = ["GID", "NAME", "GWEIGHT"]
    
    UNKNOWN = [x for x in names if x not in KNOWN_NAMES]
    
    # If only one column other than KNOWN_NAMES, then it must be the
    # one supplied by the user.  Use that as the ID.
    if len(UNKNOWN) == 1:
        return MATRIX.row_names(UNKNOWN[0])

    # If multiple columns other than KNOWN_NAMES, then look for
    # possible ones.  If none are available, then just arbitrarily
    # choose the first one, alphabetically.
    if len(UNKNOWN) > 1:
        POSSIBLE_IDS = ["Probe.Set.ID", "LocusLink", "Gene.Symbol"]
        for n in POSSIBLE_IDS:
            if n in names:
                return MATRIX.row_names(n)
        x = sorted(UNKNOWN)
        return MATRIX.row_names(x[0])
    
    # If no columns other than KNOWN_NAMES, then select the ones from
    # KNOWN_NAMES.
    for n in KNOWN_NAMES:
        if n in names:
            return MATRIX.row_names(n)

    raise AssertionError, "I could not find any possible gene IDs."

def _get_array_ids(MATRIX):
    import arrayio
    return MATRIX.col_names(arrayio.COL_ID)

def _get_gene_labels(MATRIX):
    name = _choose_gene_label(MATRIX)
    labels = MATRIX.row_names(name)[:]

    # Gene Symbols can contain "///".  If this exists, then truncate
    # it for readability.
    for i in range(len(labels)):
        j = labels[i].find("///")
        if j >= 0:
            labels[i] = labels[i][:j].strip()
            
    return labels

def _get_array_labels(MATRIX):
    import arrayio
    labels = MATRIX.col_names(arrayio.COL_ID)[:]

    # Array labels might be:
    # 122.CEL.gz
    # Cut off the .gz and .CEL for readability.
    for i in range(len(labels)):
        x = labels[i]
        x = x.replace(".gz", "")
        x = x.replace(".CEL", "")
        labels[i] = x
    return labels

_COLOR_CACHE = {}  # (fn, num) -> list
def _get_color(perc, color_fn, num_colors=256):
    # Convert a percentage into a (r, g, b) color.
    # r, g, b are numbers from 0 to 255.
    global _COLOR_CACHE
    import math

    x = color_fn, num_colors
    if x not in _COLOR_CACHE:
        _COLOR_CACHE[x] = color_fn(num_colors)
    colors = _COLOR_CACHE[x]
    i = min(int(math.floor(perc*num_colors)), num_colors-1)
    r, g, b = colors[i]
    r = min(int(math.floor(r*256)), 255)
    g = min(int(math.floor(g*256)), 255)
    b = min(int(math.floor(b*256)), 255)
    return r, g, b

def _exists_nz(filename):
    import stat
    if not os.path.exists(filename):
        return None
    if os.stat(filename)[stat.ST_SIZE] > 0:
        return filename
    return None
    
def main():
    from optparse import OptionParser, OptionGroup

    # Smallest boxes in which the labels can be clearly printed.
    DEF_WIDTH = 20
    DEF_HEIGHT = 20

    usage = "usage: %prog [options] filename\n" + \
            "  Values should range from -1 to 1."
    parser = OptionParser(usage=usage, version="%prog 02")

    parser.add_option(
        "--cluster", dest="autoclust", default=False, action="store_true",
        help="Will automatically use the options: --gc median --gn var "
        "-g -a --gl --al")
    parser.add_option(
        "", "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_option(
        "", "--cluster_app", dest="cluster", default=None,
        help="Path to cluster program.")

    group = OptionGroup(
        parser, "Input/Output",
        "If not specified, I will search the current directory for "
        "reasonable defaults.")
    group.add_option(
        "-i", "--gene_indexes", dest="gene_indexes", type="string",
        default=None,
        help="Indexes of genes to show, e.g. 1-50,75 (1-based, inclusive).")
    group.add_option(
        "-n", "--gene_names", dest="gene_names", type="string",
        action="append", default=[],
        help="Comma-separated list of genes to show.")
    group.add_option(
        "", "--gene_file", dest="gene_file", type="string", default=None,
        help="<file>[,<column num>] containing the names of genes.")
    group.add_option(
        "", "--array_indexes", dest="array_indexes", type="string",
        default=None,
        help="Indexes of arrays to show, e.g. 1-50,75 (1-based, inclusive).")
    group.add_option(
        "", "--array_file", dest="array_file", type="string", default=None,
        help="<file>[,<column num>] containing the names of arrays.")
    group.add_option(
        "-j", "--jobname", dest="jobname", type="string", default=None,
        help="Save the processed matrix to a file.")
    group.add_option(
        "-o", "--outfile", dest="outfile", type="string", default=None,
        help="Save the image to this file.")
    group.add_option(
        "--format", dest="image_format", type="choice",
        choices=["png", "svg"], default="png",
        help="Image format: png (default) or svg.")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Normalization")
    group.add_option(
        "-l", "--log_transform", dest="log_transform", default=False,
        action="store_true",
        help="Log transform the data first.")
    group.add_option(
        "--select_genes_var", dest="select_genes_var", type="int",
        default=None,
        help="Select this number of genes based on variance.")
    group.add_option(
        "--gc", "--gene_center", dest="gene_center", type="choice",
        choices=["mean", "median"], default=None, 
        help="Center each gene by: mean, median.")
    group.add_option(
        "--gn", "--gene_normalize", dest="gene_normalize", default=None,
        choices=["ss", "var"], 
        help="Normalize each gene by: ss (sum of squares), var (variance).")
    group.add_option(
        "--ac", "--array_center", dest="array_center", type="choice",
        choices=["mean", "median"], default=None, 
        help="Center each array by: mean, median.")
    group.add_option(
        "--an", "--array_normalize", dest="array_normalize", default=None,
        choices=["ss", "var"], 
        help="Normalize each array by: ss (sum of squares), var (variance).")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Clustering")
    group.add_option(
        "-g", "--cluster_genes", dest="cluster_genes", default=False,
        action="store_true", help="Cluster the genes.")
    group.add_option(
        "-a", "--cluster_arrays", dest="cluster_arrays", default=False,
        action="store_true", help="Cluster the arrays.")
    group.add_option(
        "--algorithm", dest="cluster_alg", type="choice",
        choices=["hierarchical", "kmeans"], default="hierarchical",
        help="Choose a clustering algorithm: hierarchical (default), "
        "kmeans.")
    group.add_option(
        "--distance", dest="distance", type="choice",
        choices=["uncent-cor", "pearson", "abs-uncent-cor", "abs-pearson",
                 "spearman", "kendall", "euclidean", "city-block"],
        default="uncent-cor", 
        help="Choose a distance metric: uncent-cor (default), pearson, "
        "abs-uncent-cor, abs-pearson, spearman, kendall, euclidean, "
        "or city-block.")
    group.add_option(
        "--method", dest="method", type="choice",
        choices=["complete", "single", "centroid", "average"],
        default="complete",
        help="Choose clustering method: complete (default), single, "
        "centroid, or average linkage.")
    group.add_option(
        "--gk", "--gene_k", dest="gene_k", type="int", default=None,
        help="For hierarchical clustering, cut genes into K clusters "
        "(default 0).")
    group.add_option(
        "--ak", "--array_k", dest="array_k", type="int", default=None,
        help="For hierarchical clustering, cut arrays into K clusters "
        "(default 0).")
    group.add_option(
        "-k", "--kmeans_k", dest="kmeans_k", type="int", default=5,
        help="For K-means clustering, choose K (default 5).")
    parser.add_option_group(group)

    group = OptionGroup(parser, "Graphics")
    group.add_option(
        "--gl", "--label_genes", dest="label_genes", action="store_true",
        default=False, help="Label the genes on the plot.")
    group.add_option(
        "--al", "--label_arrays", dest="label_arrays", action="store_true",
        default=False, help="Label the arrays on the plot.")
    group.add_option(
        "-x", "--width", dest="width", type="int", default=DEF_WIDTH,
        help="Width of boxes (default=%d)." % DEF_WIDTH)
    group.add_option(
        "-y", "--height", dest="height", type="int", default=DEF_HEIGHT,
        help="Height of boxes (default=%d)." % DEF_HEIGHT)
    group.add_option(
        "-s", "--scale", dest="scale", default=0, type="float",
        help="Add this to each expression value before plotting (default 0)."
        "  Scale applied, then gain applied.")
    group.add_option(
        "-m", "--gain", dest="gain", default=1, type="float",
        help="Multiply each expression value by this value before plotting "
        "(default 1).")
    group.add_option(
        "", "--no_autoscale", dest="autoscale", action="store_false",
        default=True, help="Disable autoscaling.")
    group.add_option(
        "--color", dest="color_scheme", type="choice", default="matlab",
        choices=["red", "red-green", "blue-yellow", "matlab", "bild",
                 "genespring", "yahoo"],
        help="Choose the color scheme to use: red, red-green, blue-yellow, "
        "matlab (default), bild, genespring, or yahoo.")
    group.add_option(
        "--gene_tree_scale", dest="gene_tree_scale", type="float", default=1.0,
        help="Scale the width of the gene tree by this factor.  "
        "Set to 0 to disable dendrogram.")
    group.add_option(
        "--array_tree_scale", dest="array_tree_scale", type="float",
        default=1.0,
        help="Scale the height of the array tree by this factor.  "
        "Set to 0 to disable dendrogram.")
    group.add_option(
        "--gene_tree_thickness", dest="gene_tree_thickness", type="float",
        default=1.0,
        help="Scale the thickness of the lines in the gene tree by this "
        "factor.")
    group.add_option(
        "--array_tree_thickness", dest="array_tree_thickness", type="float",
        default=1.0,
        help="Scale the thickness of the lines in the array tree by this "
        "factor.")
    
    parser.add_option_group(group)

    # Parse the input arguments.
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("Please specify an infile.")
    infile, = args
    if not os.path.exists(infile):
        parser.error("I could not find file %s." % infile)
        
    if options.libpath:
        sys.path = options.libpath + sys.path

    if options.autoclust:
        options.gene_center = "median"
        options.gene_normalize = "var"
        options.cluster_genes = True
        options.cluster_arrays = True
        options.label_genes = True
        options.label_arrays = True

    if not options.jobname:
        options.jobname = _guess_filestem(infile)

    outfile = options.outfile
            
    # Choose a plotting library.
    if options.image_format == "svg":
        plotlib = __import__(
            "genomicode.svgplot", globals(), locals(), ["svgplot"])
        if not outfile:
            outfile = "%s.svg" % options.jobname
    else:
        plotlib = __import__(
            "genomicode.pilplot", globals(), locals(), ["pilplot"])
        if not outfile:
            outfile = "%s.png" % options.jobname

    MATRIX, cluster_data = read_data_set(infile)
    # Do the normalization, clustering, etc. before plotting the
    # results.
    MATRIX, cluster_data = process_data_set(
        MATRIX, options.cluster, cluster_data, options.jobname,
        options.gene_indexes, options.gene_names, options.gene_file,
        options.select_genes_var,
        options.array_indexes, options.array_file,
        options.log_transform, options.gene_center, options.gene_normalize,
        options.array_center, options.array_normalize,
        options.cluster_genes, options.cluster_arrays,
        options.cluster_alg, options.distance, options.method,
        options.gene_k, options.array_k, options.kmeans_k,
        options.scale, options.gain, options.autoscale)
    layout = make_layout(
        MATRIX, cluster_data, plotlib, 
        options.width, options.height, options.color_scheme,
        options.cluster_genes, options.gene_tree_scale,
        options.gene_tree_thickness,
        options.cluster_arrays, options.array_tree_scale,
        options.array_tree_thickness,
        options.cluster_alg, 
        options.label_genes, options.label_arrays)
    
    megapixels = layout.heatmap.width() * layout.heatmap.height() / 1024 / 1024
    assert megapixels <= MAX_MEGAPIXELS, "%dx%d plot too big." % (
        layout.heatmap.width(), layout.heatmap.height())
    
    coords = calc_coords_for_layout(layout)
    plot(outfile, MATRIX, cluster_data, plotlib, layout, coords)

if __name__ == '__main__':
    main()
    #import profile; profile.run("main()")
