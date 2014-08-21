#!/usr/bin/env python

# To do:
# - Allow choice of gene label.
# - Make transparent background.
# - What if there are nan's in the probability file?

"""

Classes:
ClusterData

PlotLayout              Holds layouts for each individual pieces of the plot.
PlotCoords

HeatmapLayout
ColorbarLayout
DendrogramLayout
GeneDendrogramLayout
ArrayDendrogramLayout
GeneClusterLayout
ArrayClusterLayout
GeneLabelLayout
ArrayLabelLayout


Global Variables:
MIN_FONTSIZE
MAX_MEGAPIXELS


Functions:
process_data_set
make_layout
calc_coords_for_layout

plot
plot_matrix
plot_colorbar
plot_dendrogram
plot_gene_clusters
plot_array_clusters
plot_gene_labels
plot_array_labels

read_data_set
read_filecol

_pretty_scale_matrix
_calc_colorbar_size
_calc_colorbar_ticks
_get_color

_choose_gene_id
_choose_gene_label
_get_gene_ids        Return unique IDs for the genes.
_get_array_ids       Return unique IDs for the arrays.
_get_gene_labels     Return pretty labels for the genes.
_get_array_labels    Return pretty labels for the arrays.

_parse_gene_names
_parse_color

_exists_nz

"""

import os, sys

MIN_FONTSIZE = 6
MAX_MEGAPIXELS = 256  # No more than 256 megapixel plot.


class ClusterData:
    def __init__(
        self, gene_tree, array_tree, gene_tree_cluster, array_tree_cluster,
        gene_cluster, array_cluster):
        # gene_tree           gtr
        # array_tree          atr
        # gene_tree_cluster   gtc
        # array_tree_cluster  atc
        # gene_cluster        kgg
        # array_cluster       kag
        self.gene_tree = gene_tree
        self.array_tree = array_tree
        self.gene_tree_cluster = gene_tree_cluster
        self.array_tree_cluster = array_tree_cluster
        self.gene_cluster = gene_cluster
        self.array_cluster = array_cluster



class PlotLayout:
    def __init__(self, heatmap, colorbar, gene_dendrogram, array_dendrogram,
                 gene_cluster, array_cluster, gene_label, array_label):
        self.heatmap = heatmap
        self.colorbar = colorbar
        self.gene_dendrogram = gene_dendrogram
        self.array_dendrogram = array_dendrogram
        self.gene_cluster = gene_cluster
        self.array_cluster = array_cluster
        self.gene_label = gene_label
        self.array_label = array_label
        

class PlotCoords:
    def __init__(self, hm_x, hm_y, cb_x, cb_y, gd_x, gd_y, ad_x, ad_y,
                 gc_x, gc_y, ac_x, ac_y, gl_x, gl_y, al_x, al_y):
        self.hm_x, self.hm_y = hm_x, hm_y
        self.cb_x, self.cb_y = cb_x, cb_y
        self.gd_x, self.gd_y = gd_x, gd_y
        self.ad_x, self.ad_y = ad_x, ad_y
        self.gc_x, self.gc_y = gc_x, gc_y
        self.ac_x, self.ac_y = ac_x, ac_y
        self.gl_x, self.gl_y = gl_x, gl_y
        self.al_x, self.al_y = al_x, al_y
        

class HeatmapLayout:
    def __init__(
        self, nrow, ncol, boxwidth, boxheight, scale_border, grid,
        inverse_colors, black0, color_fn):
        # Looks OK with even 1 pixel.
        #MIN_GRID = 1
        #if boxwidth < MIN_GRID or boxheight < MIN_GRID:
        #    grid = False
        self.nrow = nrow
        self.ncol = ncol
        self.boxwidth = boxwidth
        self.boxheight = boxheight
        self.inverse_colors = inverse_colors
        self.black0 = black0
        self.color_fn = color_fn
        self.BORDER = int(round(min(boxwidth, boxheight)*0.10) * scale_border)
        self.GRID_SIZE = int(round(min(boxwidth, boxheight)*0.05))
        if not grid:
            self.GRID_SIZE = 0
        assert self.GRID_SIZE <= self.BORDER
    def width(self):
        return self.size()[0]
    def height(self):
        return self.size()[1]
    def size(self):
        height = self.BORDER*2
        width = self.BORDER*2
        height += self.boxheight * self.nrow
        width += self.boxwidth * self.ncol
        height += (self.nrow-1) * self.GRID_SIZE
        width += (self.ncol-1) * self.GRID_SIZE
        return width, height
    def coord(self, row, col):
        x = self.BORDER
        y = self.BORDER
        x += col * (self.boxwidth + self.GRID_SIZE)
        y += row * (self.boxheight + self.GRID_SIZE)
        return x, y, self.boxwidth, self.boxheight
    def color(self, x):
        # x is from [0, 1].  find the nearest color.
        import math
        if x is None or math.isnan(x):
            # Missing value.  Return a white box.
            #return _get_color(0.5, self.color_fn)
            return (255, 255, 255)
        assert x >= 0 and x <= 1, "x out of range: %g" % x
        return _get_color(
            x, self.color_fn, flip_colors=self.inverse_colors,
            black0=self.black0)


class ColorbarLayout:
    def __init__(
        self, cb_width, cb_height, signal_0, signal_1,
        ticks, tick_labels, label_sizes, fontsize, inverse_colors, color_fn):
        TICK_SIZE = 0.15       # relative to BAR_SHORT
        TICK_BUFFER = 0.15     # relative to BAR_SHORT

        assert len(ticks) == len(tick_labels)
        assert len(ticks) == len(label_sizes)

        self.TICK_SIZE = TICK_SIZE
        self.TICK_BUFFER = TICK_BUFFER
        self._cb_width = cb_width  # Width of just the bar.
        self._cb_height = cb_height
        self._signal_0 = signal_0
        self._signal_1 = signal_1
        self._ticks = ticks
        self._tick_labels = tick_labels
        self._label_sizes = label_sizes
        self._fontsize = fontsize
        self._inverse_colors = inverse_colors
        #self._black0 = black0
        self._color_fn = color_fn
    def is_vertical(self):
        return self._cb_height > self._cb_width
    def width(self):
        width = self._cb_width
        if self.is_vertical():
            # Vertical skinny colorbar.
            # Tick mark.
            width += self._cb_width * self.TICK_SIZE
            # BUFFER between tick mark and label.
            width += self._cb_width * self.TICK_BUFFER
            # Text.
            text_width = max([x[1] for x in self._label_sizes])
            # PIL doesn't calculate text widths very accurately.
            # Compensate with a fudge factor.  2 is not big enough.
            text_width *= 2.5
            width += text_width
        width = int(width)
        return width
    def height(self):
        height = self._cb_height
        # Bug: For vertical colorbar, does not take into account
        # height of labels.  Might be cut off.
        if not self.is_vertical():
            # Horizontal colorbar.
            # Tick mark.
            height += self._cb_height * self.TICK_SIZE
            # BUFFER between tick mark and label.
            height += self._cb_height * self.TICK_BUFFER
            # Text.
            text_height = max([x[0] for x in self._label_sizes])
            height += text_height
        height = int(height)
        return height
    def size(self):
        # Size taken by the entire color bar, including labels.
        return self.width(), self.height()
    def bar_width(self):
        # Size of just the bar.
        return self._cb_width
    def bar_height(self):
        return self._cb_height
    def num_ticks(self):
        return len(self._tick_labels)
    def tick_coord(self, i):
        assert i >= 0 and i < len(self._ticks)
        tick = self._ticks[i]
        perc = float(tick-self._signal_0)/(self._signal_1-self._signal_0)
        if self.is_vertical():
            width = self._cb_width * self.TICK_SIZE
            height = 1
            x = self._cb_width
            #y = perc * self._cb_height          # high numbers on bottom
            y = (1.0-perc) * self._cb_height    # high numbers on top
            y = min(y, self._cb_height-height)
        else:
            width = 1
            height = self._cb_height * self.TICK_SIZE
            x = perc * self._cb_width
            y = self._cb_height
            x = min(x, self._cb_width-width)
        x, y, width, height = int(x), int(y), int(width), int(height)
        return x, y, width, height
    def tick_label(self, i):
        assert i >= 0 and i < len(self._tick_labels)
        return self._tick_labels[i]
    def label_coord(self, i):
        x = self.tick_coord(i)
        tick_x, tick_y, tick_width, tick_height = x
        label_width, label_height = self._label_sizes[i]
        if self.is_vertical():
            x = tick_x + tick_width + self._cb_width*self.TICK_BUFFER
            y = tick_y - label_height/2.0
        else:
            x = tick_x - label_width/2.0
            y = tick_y + tick_height + self._cb_height*self.TICK_BUFFER
        x, y = int(x), int(y)
        return x, y
    def label_size(self, i):
        return self._label_sizes[i]
    def fontsize(self):
        return self._fontsize
    def color(self, x):
        # x is from [0, 1].  find the nearest color.
        assert x >= 0 and x <= 1, "x out of range: %g" % x
        return _get_color(x, self._color_fn, flip_colors=self._inverse_colors)


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
        #x = max(x1, x2)
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
        #min_ppi = pixels_per_item
        x = int(math.ceil(min_ppi*0.20 * thickness_scale))
        #x = min(max(x, 1), min_ppi)
        x = min(max(x, 1), pixels_per_item)
        self.LINEWIDTH = x
        #print min_ppi, thickness_scale, min_ppi, self.LINEWIDTH
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
            # If requested, should I use the inverse of the color here?
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
        height = self.item_size() + self.LINEWIDTH
        width = self.dist_size() + self.LINEWIDTH
        return width, height
    def coord(self, row, distance):
        x = self.dist_coord(distance)
        y = self.item_coord(row)
        return x, y
    def lines(self, node_num, node_dist, left_num, left_dist,
              right_num, right_dist):
        node_x, node_y = self.coord(node_num, node_dist)
        left_x, left_y = self.coord(left_num, left_dist)
        right_x, right_y = self.coord(right_num, right_dist)
        #print "NODE", node_x, node_y
        #print "LEFT", left_x, left_y
        #print "RIGHT", right_x, right_y

        # The right node is on top of the left node.
        #  3-----4   right  (2 lines: vertical and horizontal)
        #  |
        #  *         node
        #  |
        #  1-----2   left   (2 lines: vertical and horizontal)
        #
        # Separate out lines 1 and 3 in case they are different
        # colors.
        # 
        # Add two left lines, then two right lines.
        line1 = self.vthicken(node_x, node_y, 1, left_y-node_y+1)
        line2 = self.hthicken(node_x, left_y, left_x-node_x+1, 1)
        line3 = self.vthicken(node_x, right_y, 1, node_y-right_y+1)
        line4 = self.hthicken(node_x, right_y, right_x-node_x+1, 1)
        
        # For some reason, left and right are sometimes switched.  If
        # that's the case, then adjust the coordinates so there are no
        # negative heights.
        if line1[3] < 0:
            line1 = self.vthicken(node_x, left_y, 1, node_y-left_y+1)
        if line3[3] < 0:
            line3 = self.vthicken(node_x, node_y, 1, right_y-node_y+1)
        assert line1[3] >= 0, "%d %d %d" % (node_x, left_x, right_x)
        assert line3[3] >= 0, "%d %d %d" % (node_x, left_x, right_x)

        # Make sure the x-coordinates of 2,4 are aligned with 1,3.
        # Also, make sure the width is at least the same as the line
        # width.
        if line1[0] < line2[0]:
            delta = line2[0] - line1[0]
            x, y, width, height = line2
            line2 = x-delta, y, max(width+delta, self.LINEWIDTH), height
            x, y, width, height = line4
            line4 = x-delta, y, max(width+delta, self.LINEWIDTH), height
            
        lines = [line1, line2, line3, line4]
        return lines
    def root(self, node_num, node_dist):
        root_x, root_y = self.coord(node_num, self.lowest)
        node_x, node_y = self.coord(node_num, node_dist)
        x = root_x, root_y, node_x-root_x+1, 1
        return self.hthicken(*x)


class ArrayDendrogramLayout(DendrogramLayout):
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
        height = self.dist_size() + self.LINEWIDTH
        width = self.item_size() + self.LINEWIDTH
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

        #  1--*--3
        #  |     |
        #  2     4
        line1 = self.hthicken(left_x, node_y, node_x-left_x+1, 1)
        line2 = self.vthicken(left_x, node_y, 1, left_y-node_y+1)
        line3 = self.hthicken(node_x, node_y, right_x-node_x+1, 1)
        line4 = self.vthicken(right_x, node_y, 1, right_y-node_y+1)

        # For some reason, left and right are sometimes switched.  If
        # that's the case, then adjust the coordinates so there are no
        # negative widths.
        if line1[2] < 0:
            line1 = self.hthicken(node_x, node_y, left_x-node_x+1, 1)
        if line3[2] < 0:
            line3 = self.hthicken(right_x, node_y, node_x-right_x+1, 1)
        assert line1[2] >= 0, "%d %d %d" % (node_x, left_x, right_x)
        assert line3[2] >= 0, "%d %d %d" % (node_x, left_x, right_x)

        # Make sure the y-coordinates of 2,4 are aligned with 1,3.
        # Also, make sure the height is at least the same as the line
        # width.
        if line1[1] < line2[1]:
            delta = line2[1] - line1[1]
            x, y, width, height = line2
            line2 = x, y-delta, width, max(height+delta, self.LINEWIDTH)
            x, y, width, height = line4
            line4 = x, y-delta, width, max(height+delta, self.LINEWIDTH)
            
        #print node_x, node_y
        lines = [line1, line2, line3, line4]
        return lines
    def root(self, node_num, node_dist):
        root_x, root_y = self.coord(node_num, self.lowest)
        node_x, node_y = self.coord(node_num, node_dist)
        x = root_x, root_y, 1, node_y-root_y+1
        return self.vthicken(*x)


class GeneClusterLayout:
    def __init__(self, num_items, item_width, item_height, grid):
        array_layout = ArrayClusterLayout(
            num_items, item_height, item_width, grid)
        self.array_layout = array_layout
    def width(self):
        return self.size()[0]
    def height(self):
        return self.size()[1]
    def size(self):
        height, width = self.array_layout.size()
        return width, height
    def coord(self, num):
        y, x, height, width = self.array_layout.coord(num)
        return x, y, width, height


class ArrayClusterLayout:
    def __init__(self, num_items, item_width, item_height, grid):
        self.num_items = num_items
        self.item_width = item_width
        self.item_height = item_height
        self.BORDER = 1
        self.GRID_SIZE = 1
        if not grid:
            self.GRID_SIZE = 0
        assert self.GRID_SIZE <= self.BORDER
    def width(self):
        return self.size()[0]
    def height(self):
        return self.size()[1]
    def size(self):
        height = self.BORDER*2
        width = self.BORDER*2
        height += self.item_height
        width = self.item_width * self.num_items
        width += (self.num_items-1) * self.GRID_SIZE
        return width, height
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self.num_items
        x = self.BORDER
        y = self.BORDER
        x += num * (self.item_width + self.GRID_SIZE)
        return x, y, self.item_width, self.item_height


class GeneLabelLayout:
    def __init__(self, item_height, item_widths, fontsize):
        self._item_height = item_height
        self._item_widths = item_widths
        num_items = len(item_widths)
        self._width = max(item_widths)
        self._height = item_height * num_items
        self._num_items = num_items
        self._fontsize = fontsize
    def item_height(self):
        return self._item_height
    def width(self):
        return self._width
    def height(self):
        return self._height
    def size(self):
        return self.width(), self.height()
    def fontsize(self):
        return self._fontsize
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self._num_items
        x = 0
        y = num * self._item_height
        return x, y, self._item_widths[num], self._item_height


class ArrayLabelLayout:
    def __init__(self, item_height, item_widths, fontsize):
        # item_height refers to the text, not rotated.
        self._item_height = item_height
        self._item_widths = item_widths
        num_items = len(item_widths)
        # _width and _height refer to the layout object.
        self._width = item_height * num_items
        self._height = max(item_widths)
        self._num_items = num_items
        self._fontsize = fontsize
    def item_height(self):
        return self._item_height
    def width(self):
        return self._width
    def height(self):
        return self._height
    def size(self):
        return self.width(), self.height()
    def fontsize(self):
        return self._fontsize
    def coord(self, num):
        # Return a box that bounds the region.
        assert num >= 0 and num < self._num_items
        x = num * self._item_height
        y = self._height-self._item_widths[num]
        return x, y, self._item_height, self._item_widths[num]


def process_data_set(MATRIX, scale, gain, autoscale):
    assert MATRIX.nrow() > 0, "Matrix has no genes."
    
    x = _pretty_scale_matrix(MATRIX, scale, gain, autoscale)
    MATRIX_scaled, orig_min, orig_max = x
    return MATRIX_scaled, orig_min, orig_max


def make_layout(
    MATRIX, cluster_data, plotlib,
    # Heatmap
    boxwidth, boxheight, scale_border, grid, color_scheme, flip_colors,
    signal_0, signal_1, black0,
    # Labels
    label_genes, label_arrays,
    # Dendrograms
    gene_tree_scale, gene_tree_thickness,
    array_tree_scale, array_tree_thickness,
    # Colorbar
    colorbar,
    ):
    from genomicode import colorlib

    # Choose the color scheme.
    scheme2fn = {
        "red" : colorlib.red_shade,
        "white" : colorlib.white_shade,
        "red-green" : colorlib.rg_array_colors,
        "blue-yellow" : colorlib.by_array_colors,
        "matlab" : colorlib.matlab_colors,
        "bild" : colorlib.bild_colors,
        "genespring" : colorlib.genespring_colors,
        "yahoo" : colorlib.yahoo_weather_colors,
        }
    assert color_scheme in scheme2fn, "Unknown color scheme: %s" % color_scheme
    color_fn = scheme2fn[color_scheme]

    # Make the layout for the heatmap.
    hm_layout = HeatmapLayout(
        MATRIX.nrow(), MATRIX.ncol(), boxwidth, boxheight,
        scale_border, grid, flip_colors, black0, color_fn)

    # Make the layout for the colorbar.
    cb_layout = None
    if colorbar:
        x = _calc_colorbar_size(
            hm_layout.width(), hm_layout.height(), hm_layout.GRID_SIZE,
            boxwidth, boxheight)
        width, height = x
        x = _calc_colorbar_ticks(
            width, height, signal_0, signal_1, plotlib)
        ticks, tick_labels, label_sizes, fontsize = x
        cb_layout = ColorbarLayout(
            width, height, signal_0, signal_1,
            ticks, tick_labels, label_sizes, fontsize, flip_colors, color_fn)

    # Make layouts for the dendrograms.
    gd_layout = ad_layout = None
    #if(cluster_genes and cluster_data.gene_tree and gene_tree_scale > 0 and 
    #   cluster_alg == "hierarchical" and MATRIX.ncol() > 1):
    if(cluster_data.gene_tree and gene_tree_scale > 0 and MATRIX.ncol() > 1):
        # Only add the dendrogram if hierarchical clustering was
        # requested.  If clustering not done, then the matrix file
        # will not have the GID annotations, and there will be no way
        # to match up the genes with the clusters.
        # Also should add dendrogram if the clusters were supplied by
        # the user in a gtr file.
        #print "Making gd_layout."
        assert gene_tree_scale > 0
        assert gene_tree_thickness > 0
        width, height = boxwidth, boxheight
        width += hm_layout.GRID_SIZE
        height += hm_layout.GRID_SIZE
        gd_layout = GeneDendrogramLayout(
            MATRIX.nrow(), MATRIX.ncol(), width, height,
            gene_tree_scale, gene_tree_thickness, 
            cluster_data.gene_tree, cluster_data.gene_tree_cluster,
            colorlib.matlab_colors)
    if(cluster_data.array_tree and array_tree_scale > 0 and MATRIX.nrow() > 1):
        #print "Making ad_layout."
        assert array_tree_scale > 0
        assert array_tree_thickness > 0
        width, height = boxwidth, boxheight
        width += hm_layout.GRID_SIZE
        height += hm_layout.GRID_SIZE
        #print "HERE", width, height
        ad_layout = ArrayDendrogramLayout(
            MATRIX.nrow(), MATRIX.ncol(), width, height,
            array_tree_scale, array_tree_thickness, 
            cluster_data.array_tree, cluster_data.array_tree_cluster, 
            colorlib.matlab_colors)

    # Make layouts for the clusters.
    # Can plot these (k-means) clusters if either kmeans or
    # hierarchical clustering was requested.  Unlike hierarchical
    # clustering, plotting this does not require any annotations in
    # the matrix file.
    gc_layout = ac_layout = None
    if cluster_data.gene_cluster:
        gc_layout = GeneClusterLayout(MATRIX.nrow(), boxwidth, boxheight, grid)
    if cluster_data.array_cluster:
        # Make these variables.
        arraycluster_width, arraycluster_height = boxwidth, 20
        ac_layout = ArrayClusterLayout(
            MATRIX.ncol(), arraycluster_width, arraycluster_height, grid)

    # Make the layout for the gene or array labels.
    gl_layout = al_layout = None
    gene_labels = array_labels = None
    gl_fontsize = al_fontsize = None
    # If plotting both gene and array labels, make sure they aren't
    # wildly different sizes.
    if label_genes:
        gl_fontsize = plotlib.fit_fontsize_to_height(boxheight)
        if gl_fontsize < MIN_FONTSIZE:
            gl_fontsize = None
    if label_arrays:
        al_fontsize = plotlib.fit_fontsize_to_height(boxwidth)
        if al_fontsize < MIN_FONTSIZE:
            al_fontsize = None
    if gl_fontsize and al_fontsize:
        FONT_RATIO = 1.5
        gl_fontsize = int(min(gl_fontsize, al_fontsize*FONT_RATIO))
        al_fontsize = int(min(al_fontsize, gl_fontsize*FONT_RATIO))
    
    if label_genes and gl_fontsize:
        gene_labels = _get_gene_labels(MATRIX)
        height = boxheight
        height += hm_layout.GRID_SIZE
        widths = [plotlib.get_text_size(x, gl_fontsize)[0]
                  for x in gene_labels]
        gl_layout = GeneLabelLayout(height, widths, gl_fontsize)
    if label_arrays and al_fontsize:
        array_labels = _get_array_labels(MATRIX)
        width = boxwidth
        width += hm_layout.GRID_SIZE
        widths = [plotlib.get_text_size(x, al_fontsize)[0]
                  for x in array_labels]
        al_layout = ArrayLabelLayout(width, widths, al_fontsize)

    x = PlotLayout(
        hm_layout, cb_layout, gd_layout, ad_layout, gc_layout, ac_layout,
        gl_layout, al_layout)
    return x


def calc_coords_for_layout(layout):
    # Layout:
    #
    #                           <array dendrogram>
    #                           <array buffer 1>
    #                           <array label>
    #                           <array buffer 2>
    #                           <array cluster>
    #                           <array buffer 3>
    # <dend> <label> <cluster>  <heatmap>           <colobar> (if on right)
    #                           <colorbar> (if on bottom)

    x = y = 0
    def _safe_size(layout):
        if layout is None:
            return 0, 0
        return layout.size()

    # cluster should be spaced a little further
    gb1, gb2, gb3 = 6, 6, 10
    ab1, ab2, ab3 = 6, 6, 10
    if not layout.gene_dendrogram:
        gb1 = 0
    if not layout.gene_label:
        gb2 = 0
    if not layout.gene_cluster:
        gb3 = 0
    if not layout.array_dendrogram:
        ab1 = 0
    if not layout.array_label:
        ab2 = 0
    if not layout.array_cluster:
        ab3 = 0

    hm_width, hm_height = _safe_size(layout.heatmap)
    cb_width, cb_height = _safe_size(layout.colorbar)
    gd_width, gd_height = _safe_size(layout.gene_dendrogram)
    ad_width, ad_height = _safe_size(layout.array_dendrogram)
    gc_width, gc_height = _safe_size(layout.gene_cluster)
    ac_width, ac_height = _safe_size(layout.array_cluster)
    gl_width, gl_height = _safe_size(layout.gene_label)
    al_width, al_height = _safe_size(layout.array_label)
    
    # Position the heatmap based on the dendrograms.
    hm_x = x + gd_width + gb1 + gl_width + gb2 + gc_width + gb3
    hm_y = y + ad_height + ab1 + al_height + ab2 + ac_height + ab3

    # On X-axis: gene dendrogram, label, cluster, then heatmap.
    gd_x, gd_y = x, hm_y+layout.heatmap.BORDER
    gl_x, gl_y = gd_x+gd_width+gb1, gd_y
    gc_x, gc_y = gl_x+gl_width+gb2, gd_y
    
    # On Y-axis: array dendrogram, label, cluster, then heatmap.
    ad_x, ad_y = hm_x+layout.heatmap.BORDER, y
    al_x, al_y = ad_x, ad_y+ad_height+ab1
    ac_x, ac_y = ad_x, al_y+al_height+ab2

    # Add the colorbar.
    cb_x = cb_y = None
    if layout.colorbar:
        CB_BUFFER = 0.75  # separation from heatmap, relative to BAR_SHORT
        bar_width = layout.colorbar.bar_width()
        bar_height = layout.colorbar.bar_height()
        if layout.colorbar.is_vertical():
            cb_x = hm_x + hm_width + CB_BUFFER*bar_width
            cb_y = hm_y
            # If there are no dendrograms or labels, then need to add
            # a buffer so that the labels aren't cut off.
            if not layout.array_dendrogram and not layout.array_label:
                cb_y += layout.colorbar.fontsize()
        else:
            cb_x = hm_x
            cb_y = hm_y + hm_height + CB_BUFFER*bar_height
            if not layout.gene_dendrogram and not layout.gene_label:
                cb_x += layout.colorbar.fontsize()
        cb_x, cb_y = int(cb_x), int(cb_y)

    x = PlotCoords(
        hm_x, hm_y, cb_x, cb_y, gd_x, gd_y, ad_x, ad_y,
        gc_x, gc_y, ac_x, ac_y, gl_x, gl_y, al_x, al_y)
    return x


def plot(
    filename, MATRIX, cluster_data, plotlib, layout, coords,
    border_color, grid_color):
    # Calculate the plot width and height.
    plot_width = coords.hm_x + layout.heatmap.width()
    plot_height = coords.hm_y + layout.heatmap.height()
    if layout.colorbar:
        x = coords.cb_x + layout.colorbar.width()
        plot_width = max(plot_width, x)
        x = coords.cb_y + layout.colorbar.height()
        plot_height = max(plot_height, x)

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
            layout.gene_cluster, cluster_data.gene_cluster,
            border_color, grid_color)
    if layout.array_cluster:
        plot_array_clusters(
            plotlib, image, MATRIX, coords.ac_x, coords.ac_y,
            layout.array_cluster, cluster_data.array_cluster,
            border_color, grid_color)

    if layout.gene_label:
        gene_labels = _get_gene_labels(MATRIX)
        plot_gene_labels(
            plotlib, image, MATRIX, coords.gl_x, coords.gl_y,
            layout.gene_label, gene_labels)
    if layout.array_label:
        array_labels = _get_array_labels(MATRIX)
        plot_array_labels(
            plotlib, image, MATRIX, coords.al_x, coords.al_y,
            layout.array_label, array_labels)
    plot_matrix(
        plotlib, image, MATRIX, coords.hm_x, coords.hm_y, layout.heatmap,
        border_color, grid_color)
    if layout.colorbar:
        plot_colorbar(
            plotlib, image, coords.cb_x, coords.cb_y, layout.colorbar)

    plotlib.write(image, open(filename, 'w'))

    
def plot_matrix(
    plotlib, image, MATRIX, xoff, yoff, layout, border_color, grid_color):
    # (0, 0, 0) is too dark for small box sizes.  100 looks too washed
    # out.  50-75 is about right.
    #GRID_COLOR = (0, 0, 0)
    #GRID_COLOR = (75, 75, 75)
    #BORDER_COLOR = (0, 0, 0)

    width, height = layout.size()

    # Draw the underlying grid.
    plotlib.rectangle(image, xoff, yoff, width, height, grid_color)
    
    # Draw a border around the heatmap.
    # Draw top, right, bottom, and left borders.
    #plotlib.rectangle(image, xoff, yoff, width, height, None, border_color)
    plotlib.rectangle(image, xoff, yoff, width, layout.BORDER, border_color)
    plotlib.rectangle(
        image, xoff+width-layout.BORDER, yoff, layout.BORDER, height,
        border_color)
    plotlib.rectangle(
        image, xoff, yoff+height-layout.BORDER, width, layout.BORDER,
        border_color)
    plotlib.rectangle(image, xoff, yoff, layout.BORDER, height, border_color)
    
    # Draw the actual matrix.
    X = MATRIX._X
    for i in range(MATRIX.nrow()):
        for j in range(MATRIX.ncol()):
            x = X[i][j]
            c = layout.color(x)

            # Find the coordinates and plot it.
            x, y, width, height = layout.coord(i, j)
            plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)


def plot_colorbar(plotlib, image, xoff, yoff, layout):
    #yoff += 100
    #xoff += 100
    BLACK = (0, 0, 0)
    OUTLINE_COLOR = (0, 0, 0)
    TICK_COLOR = (50, 50, 50)

    # Draw the colorbar.
    cb_width, cb_height = layout.bar_width(), layout.bar_height()
    if layout.is_vertical():
        for i in range(cb_height):
            #color = layout.color(float(i)/cb_height)        # big on bottom
            color = layout.color(1.0-(float(i)/cb_height))  # big on top
            plotlib.line(image, xoff, yoff+i, cb_width, 1, color)
    else:
        for i in range(cb_width):
            color = layout.color(float(i)/cb_width)
            plotlib.line(image, xoff+i, yoff, 1, cb_height, color)
            
    plotlib.rectangle(
        image, xoff, yoff, cb_width, cb_height, None, outline=OUTLINE_COLOR)

    # Draw tickmarks.
    for i in range(layout.num_ticks()):
        x = layout.tick_coord(i)
        x, y, width, height = x
        plotlib.line(image, xoff+x, yoff+y, width, height, TICK_COLOR)

    # Label the tickmarks.
    fontsize = layout.fontsize()
    if fontsize < MIN_FONTSIZE:
        return
    
    labels = [layout.tick_label(i) for i in range(layout.num_ticks())]
    label_sizes = [layout.label_size(i) for i in range(layout.num_ticks())]
    max_width = max([x[0] for x in label_sizes])
    #max_height = max([x[1] for x in label_sizes])
                   
    for i, label in enumerate(labels):
        x, y = layout.label_coord(i)
        # Right align the vertical colorbar.
        if cb_height > cb_width:
            width, height = label_sizes[i]
            x += max_width - width
        plotlib.text(image, xoff+x, yoff+y, label, fontsize, BLACK)
        

def plot_dendrogram(plotlib, image, MATRIX, xoff, yoff, layout, dim, tree):
    import arrayio
    from genomicode import clusterio

    if dim == "GENE":
        n = "GID"  # Use the gene ID if available.
        assert n in MATRIX.row_names(), "Gene dendrogram not available."
        ids = MATRIX.row_names(n)
    elif dim == "ARRAY":
        n = "AID"
        assert n in MATRIX.col_names(), "Array dendrogram not available."
        ids = MATRIX.col_names(n)
    else:
        raise AssertionError, "Unknown dim: %s" % dim

    # num is the row or column of the node.
    id2num = {}       # gene or node id -> num
    id2distance = {}  # gene or node id -> distance
    # Find id2num and id2distance for each of the leaves.
    for i, x in enumerate(ids):
        id_ = clusterio.parse_node(x)
        id2num[id_] = i
        id2distance[id_] = 1
    #print tree

    # Set id2num and id2distance the internal nodes.
    for i, node in enumerate(tree):
        id_ = -(i+1)
        left, right, distance = node
        left_num = id2num[left]
        right_num = id2num[right]
        id2num[id_] = (left_num + right_num)/2.0
        id2distance[id_] = distance
    #print id2num

    # Draw the nodes of the tree.
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
        #print node
        #print lines

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
        #data = [(ll1, cl1), (lr1, cr1)]
        #data = [(ll2, cl2), (lr2, cr2)]
        for line, color in data:
            x, y, width, height = line
            #print x, y, width, height
            assert width >= 0 and height >= 0, "%d %d" % (width, height)
            plotlib.rectangle(image, x+xoff, y+yoff, width, height, color)

    c = layout.color(node_id)
    x, y, width, height = layout.root(node_num, node_dist)
    plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)


def plot_gene_clusters(
    plotlib, image, X, xoff, yoff, layout, clusters, border_color, grid_color):
    import arrayio
    from genomicode import colorlib
    assert X.nrow() == len(clusters), "%d %d" % (X.nrow(), len(clusters))

    #GRID_COLOR = (75, 75, 75)
    #BORDER_COLOR = (0, 0, 0)

    # Figure out what kind of IDs to use.
    ID_NAMES = ["GID", "NAME", arrayio.ROW_ID]
    ID_NAMES = [x for x in ID_NAMES if x in X.row_names() or x in X._synonyms]
    ids = [x[0] for x in clusters]
    for ID_NAME in ID_NAMES:
        ID = X.row_names(ID_NAME)
        num_found = 0
        for id_ in ids:
            if id_ in ID:
                num_found += 1
        #print ID_NAME, num_found, len(ids), ids[:3]
        if num_found == len(ids):
            break
    else:
        raise AssertionError, "I could not find the cluster IDs: %s" % \
              str(X.row_names())
    
    GID = X.row_names(ID_NAME)
    gid2i = {}
    for i, gid in enumerate(GID):
        gid2i[gid] = i
        
    # Draw the underlying grid, and a border around the whole thing.
    width, height = layout.size()
    plotlib.rectangle(image, xoff, yoff, width, height, grid_color)
    plotlib.rectangle(image, xoff, yoff, width, height, None, border_color)

    max_cluster = max([x[1] for x in clusters])
    for gid, n in clusters:
        i = gid2i[gid]
        x, y, width, height = layout.coord(i)
        c = 255, 255, 255
        if n is not None:
            p = 0.5
            if max_cluster > 0:
                p = float(n) / max_cluster
            # Bug: This should be set in the layout.
            c = _get_color(p, colorlib.matlab_colors)
        plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)


def plot_array_clusters(
    plotlib, image, X, xoff, yoff, layout, clusters, border_color, grid_color):
    import arrayio
    from genomicode import colorlib
    
    #assert X.ncol() == len(clusters)

    #GRID_COLOR = (75, 75, 75)
    #BORDER_COLOR = (0, 0, 0)

    # Figure out what kind of IDs to use.
    ID_NAMES = [
        "AID", arrayio.COL_ID, arrayio.tab_delimited_format.SAMPLE_NAME]
    ID_NAMES = [x for x in ID_NAMES if x in X.col_names() or x in X._synonyms]
    #ID_NAMES = [x for x in ID_NAMES if x in X.col_names()]
    ids = [x[0] for x in clusters]
    for ID_NAME in ID_NAMES:
        ID = X.col_names(ID_NAME)
        x = [x for x in ID if x in ids]
        # Accept this column name if any IDs match.
        #if len(x) == len(ID):
        if len(x) > 0:
            break
    else:
        raise AssertionError, "I could not find the array IDs."

    # Draw the underlying grid, and a border around the whole thing.
    width, height = layout.size()
    plotlib.rectangle(image, xoff, yoff, width, height, grid_color)
    plotlib.rectangle(image, xoff, yoff, width, height, None, border_color)

    aid2cluster = {}
    for (aid, cluster) in clusters:
        aid2cluster[aid] = cluster
    max_cluster = max([x[1] for x in clusters])

    AID = X.col_names(ID_NAME)
    for i, aid in enumerate(AID):
        x, y, width, height = layout.coord(i)
        c = 255, 255, 255   # default color, if no cluster
        n = aid2cluster.get(aid)
        if n is not None:
            p = 0.5
            if max_cluster > 0:
                p = float(n) / max_cluster
            # Bug: This should be set in the layout.
            c = _get_color(p, colorlib.matlab_colors)
        plotlib.rectangle(image, x+xoff, y+yoff, width, height, c)


def plot_gene_labels(plotlib, image, X, xoff, yoff, layout, labels):
    fontsize = layout.fontsize()
    if fontsize < MIN_FONTSIZE:
        return
    #print layout.__class__.__name__
    for i in range(X.nrow()):
        x, y, width, height = layout.coord(i)
        #print x, y, width, height, fontsize
        w, h = plotlib.get_text_size(labels[i], fontsize)
        # Right-align the text.  Need layout width, not width of item.
        x += max(layout.width()-w, 0)
        # Vertical align the text.
        y += (height - h)/2
        plotlib.text(image, xoff+x, yoff+y, labels[i], fontsize, (0, 0, 0))


def plot_array_labels(plotlib, image, X, xoff, yoff, layout, labels):
    fontsize = layout.fontsize()
    if fontsize < MIN_FONTSIZE:
        return
    for i in range(X.ncol()):
        x, y, width, height = layout.coord(i)
        w, h = plotlib.get_text_size(labels[i], fontsize)
        #print "HERE1", xoff, yoff, x, y, width, height, fontsize
        # Center the text.
        x += (width-h)/2
        plotlib.text90(image, xoff+x, yoff+y, labels[i], fontsize, (0, 0, 0))


def read_data_set(
    matrix_file, gene_cluster_file, array_cluster_file,
    gene_tree_file, array_tree_file,
    gene_tree_cluster_file, array_tree_cluster_file):
    import arrayio
    from genomicode import parselib
    from genomicode import Matrix
    from genomicode import clusterio

    MATRIX = arrayio.read(matrix_file)

    # If no gene IDs were provided, then just make some up.
    if not MATRIX.row_names():
        header = "GENE.ID"
        MATRIX._row_order.append(header)
        x = ["R%s" % x for x in parselib.pretty_range(0, MATRIX.nrow())]
        MATRIX._row_names[header] = x
        synonyms = {}
        synonyms[arrayio.ROW_ID] = header
        #MATRIX = Matrix.add_synonyms(MATRIX, synonyms)
        MATRIX._synonyms.update(synonyms)
    if not MATRIX.col_names():
        header = arrayio.tdf.SAMPLE_NAME
        MATRIX._col_order.append(header)
        x = ["C%s" % x for x in parselib.pretty_range(0, MATRIX.ncol())]
        MATRIX._col_names[header] = x
        synonyms = {}
        synonyms[arrayio.COL_ID] = header
        #MATRIX = Matrix.add_synonyms(MATRIX, synonyms)
        MATRIX._synonyms.update(synonyms)


    # Read the clustering files.
    readers = [
        ("gene_tree", gene_tree_file, clusterio.read_gtr_file),
        ("array_tree", array_tree_file, clusterio.read_atr_file),
        ("gene_tree_cluster", gene_tree_cluster_file, clusterio.read_gtc_file),
        ("array_tree_cluster", array_tree_cluster_file,
         clusterio.read_atc_file),
        ("gene_cluster", gene_cluster_file, clusterio.read_kgg_file),
        ("array_cluster", array_cluster_file, clusterio.read_kag_file),
        ]
    data = {}  # name -> output
    for (name, filename, read_fn) in readers:
        if filename is None:
            data[name] = None
            continue
        assert os.path.exists(filename), "I could not find file: %s" % filename
        data[name] = read_fn(filename)

    cluster_data = ClusterData(**data)
    return MATRIX, cluster_data


def read_filecol(filecol):
    from genomicode import iolib

    # filecol is either <filename> or <filename>,<col>.  commas
    # are not allowed in the filenames.  <col> should be 1-based
    # index.
    filename, colnum = filecol, 1
    if filecol.find(",") >= 0:
        x = filecol.split(",")
        assert len(x) == 2, "File should be specified: <filename>,<col>"
        filename, colnum = x
        colnum = int(colnum)
        assert colnum >= 1
    assert os.path.exists(filename), "could not find file %s" % filename
    data = iolib.split_tdf(open(filename).read())
    # Make sure colnum is correct.
    for x in data:
        assert colnum <= len(x)
    names = [x[colnum-1].strip() for x in data]
    names = [x for x in names if x]
    return names


def _pretty_scale_matrix(MATRIX, scale, gain, autoscale):
    # Find a good default gain value.  After scaling, values should
    # range from [-1, 1].  Then, for convenience, I will re-scale that
    # matrix to [0, 1].
    # Will change the MATRIX variable.
    import math
    from genomicode import jmath

    MATRIX = MATRIX.matrix()
    nrow, ncol = MATRIX.dim()
    X = MATRIX._X

    # Choose a default scale so that the average expression level is
    # 0.
    defscale = 0.0
    if autoscale:
        x_all = []
        for x in X:
            x_all.extend(x)
        # Use safe_mean to handle missing values.
        defscale = -jmath.safe_mean(x_all)

    # Apply the scale specified by the user.
    for i in range(nrow):
        for j in range(ncol):
            # Ignore missing values.
            if X[i][j] is None:
                continue
            X[i][j] = X[i][j] + defscale + scale

    # Choose a default gain so that the maximum expression level is 1.
    defgain = 1.0
    if autoscale:
        x_max = None
        for i in range(nrow):
            for j in range(ncol):
                # Ignore missing values.
                if X[i][j] is None or math.isnan(X[i][j]):
                    continue
                if x_max is None or abs(X[i][j]) > x_max:
                    x_max = abs(X[i][j])
        if x_max is not None:
            defgain = 1.0/x_max
        # By default, automatically multiply by 2.0, or else
        # everything is too dark (empirically).
        defgain = defgain * 2.0

    # Apply the gain specified by the user.
    for i in range(nrow):
        for j in range(ncol):
            if X[i][j] is None:
                continue
            # The gain is scaled from the default gain.
            X[i][j] = X[i][j] * defgain * gain

    #x_min, x_max = 1E9, -1E9
    #for i in range(nrow):
    #    for j in range(ncol):
    #        x_min = min(x_min, X[i][j])
    #        x_max = max(x_max, X[i][j])
    #print x_min, x_max, defgain, gain, defscale, scale
    #for x in X:
    #   print "\t".join(map(str, x))

    # Finally, rescale to [0, 1].
    for i in range(nrow):
        for j in range(ncol):
            if X[i][j] is None:
                continue
            x = X[i][j]
            x = (x + 1) * 0.5
            x = max(min(x, 1), 0)
            X[i][j] = x

    #print defgain, gain, defscale, scale
    assert not math.isnan(defgain) and not math.isnan(defscale)
    ORIG_min = (0.0*2.0 - 1.0)/(defgain*gain) - (defscale+scale)
    ORIG_max = (1.0*2.0 - 1.0)/(defgain*gain) - (defscale+scale)
    
    return MATRIX, ORIG_min, ORIG_max


def _calc_colorbar_size(hm_width, hm_height, grid_size, box_width, box_height):
    # Calculate the dimensions of the colorbar.  The size of the
    # bar should be calculated based on the size of the heatmap,
    # and also the size of the boxes in the heatmap.
    #
    # MAX_BOXES of 100 is too big for signature heatmap from pybinreg.
    BAR_LONG = 0.50     # the long dimension, relative to heatmap
    BAR_SHORT = 0.075   # short dimension, relative to long_ratio
    MAX_BOXES = 50      # Maximum boxes in the long dimension.
    MIN_BOXES = 1       # Minimum boxes in the short dimension.

    vertical = hm_height > hm_width

    if vertical:
        # x1 is the upper limit.  Do not make bigger than this.
        x1 = hm_height * BAR_LONG
        # These are both lower and upper limits.  Make the bigger one
        # of these.
        x2 = (box_height+grid_size) * MAX_BOXES
        x3 = ((box_width+grid_size)*MIN_BOXES)/BAR_SHORT
        x4 = max(x2, x3)
        height = max(min(x1, x4), 1)
        width = max(height * BAR_SHORT, 1)
    else:
        x1 = hm_width * BAR_LONG
        x2 = (box_width+grid_size) * MAX_BOXES
        x3 = ((box_height+grid_size)*MIN_BOXES)/BAR_SHORT
        x4 = max(x2, x3)
        width = max(min(x1, x4), 1)
        height = max(width * BAR_SHORT, 1)
        
    width, height = int(width), int(height)
    return width, height


def _calc_colorbar_ticks(
    cb_width, cb_height, signal_0, signal_1, plotlib):
    import math
    from genomicode import graphlib

    TEXT_SIZE = 0.75
    MAX_TICKS = 20

    vertical = cb_height > cb_width

    # Calculate the minimum and maximum number to label.
    assert not math.isnan(signal_0) and not math.isnan(signal_1)
    assert signal_0 < signal_1, "%g %g" % (signal_0, signal_1)
    delta = signal_1 - signal_0
    num_decimals = max(-int(math.floor(math.log(delta, 10))), 0)+1
    # Where to start labels.  Give a larger range than might fit
    # and shrink it later based on where the tick marks are.
    label_min = math.floor(signal_0 * 10**num_decimals)/10**num_decimals
    label_max = math.ceil(signal_1 * 10**num_decimals)/10**num_decimals
    #print signal_0, signal_1, label_min, label_max, num_decimals
    assert label_min < label_max

    # Calculate the size of the font for the labels.
    text_height = int(min(cb_width, cb_height) * TEXT_SIZE)
    fontsize = plotlib.fit_fontsize_to_height(text_height)

    # Try different number of tick marks until the labels fit in the
    # colorbar.
    # Can't have more ticks than the size of the colorbar.
    num_ticks = min(MAX_TICKS, max(cb_width, cb_height)/2)
    while num_ticks > 0:
        # Calculate the ticks and remove any that are off the scale.
        ticks = graphlib.place_ticks(label_min, label_max, num_ticks=num_ticks)
        ticks = [x for x in ticks if x >= signal_0 and x <= signal_1]
        assert ticks, "I couldn't place any tick marks."

        # Format the tick labels.
        x = [max(len(str(abs(x)%1))-2, 0) for x in ticks]
        digits = max(x)
        tick_labels = ["%.*f" % (digits, x) for x in ticks]

        # Calculate the sizes of the tick labels.
        label_sizes = [plotlib.get_text_size(x, fontsize) for x in tick_labels]

        # See if this fits.
        if vertical:
            total = sum([x[1] for x in label_sizes])
        else:
            total = sum([x[0] for x in label_sizes])
        if total < max(cb_width, cb_height)/2:
            break
        num_ticks = min(num_ticks, len(ticks))-1
    assert num_ticks, "I couldn't place any tick marks."

    return ticks, tick_labels, label_sizes, fontsize


_COLOR_CACHE = {}  # (fn, num) -> list
def _get_color(perc, color_fn, num_colors=256, black0=False,
               flip_colors=False):
    # Convert a percentage into a (r, g, b) color.
    # r, g, b are numbers from 0 to 255.
    global _COLOR_CACHE
    import math

    assert perc >= 0.0 and perc <= 1.0
    if black0 and perc < 1.0/num_colors:
        return 0, 0, 0
    if flip_colors:
        perc = 1.0 - perc
    x = color_fn, num_colors
    if x not in _COLOR_CACHE:
        _COLOR_CACHE[x] = color_fn(num_colors)
    colors = _COLOR_CACHE[x]
    i = min(int(math.floor(perc*num_colors)), num_colors-1)
    r, g, b = colors[i]
    r = min(int(math.floor(r*256)), 255)
    g = min(int(math.floor(g*256)), 255)
    b = min(int(math.floor(b*256)), 255)
    #print perc, r, g, b
    return r, g, b


def _choose_gene_id(MATRIX):
    # Given a user-specified matrix, try to pick a good unique ID for
    # the genes.
    import arrayio
    
    headers = MATRIX.row_names()

    # Prioritize some potential ones.  Don't use the standard headers,
    # e.g. arrayio.ROW_ID, so that we can preserve the user's header.
    IDS = ["Probe.Set.ID"]
    for id_ in IDS:
        if id_ in headers:
            return id_

    # If no known headers are found, then choose a standard one.
    IDS = [arrayio.AFFY_PROBESET_ID, arrayio.GENE_ID, arrayio.ROW_ID]
    for id_ in IDS:
        if id_ in headers:
            return id_

    # If no standard ones are found, then just arbitrarily use the
    # first column that is not missing any values.
    for header in headers:
        names = MATRIX.row_names(header)
        missing = [x for x in names if not x.strip()]
        if not missing:
            return header
    
    raise AssertionError, "I could not find an ID for the matrix."


def _choose_gene_label(MATRIX):
    import arrayio

    names = MATRIX.row_names()

    # Prioritize some potential ones.
    IDS = [
        arrayio.GENE_SYMBOL, "Gene.Symbol", "Gene Symbol", "Symbol",
        #arrayio.GENE_DESCRIPTION, "Description",
        "DESCRIPTION",       # For GCT files.  Use the pretty name.
        "NAME",
        arrayio.GENE_ID, "LocusLink",
        arrayio.AFFY_PROBESET_ID, "Probe.Set.ID",
        arrayio.ROW_ID
        ]
    # Exception: If the GCT files have generic descriptions,
    # e.g. DESC0001, then use the name field instead.
    if "DESCRIPTION" in names:
        desc = MATRIX.row_names("DESCRIPTION")
        if desc[0].startswith("DESC"):
            i = IDS.index("DESCRIPTION")
            IDS.pop(i)
    for id_ in IDS:
        if id_ in names:
            return id_
    if names:
        return names[0]
    raise AssertionError, "I could not find an ID for the matrix."


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


def _parse_gene_names(gene_name_list):
    # This can a list of comma separated genes, e.g.
    # ["E2F1", "E2F2,E2F3"]
    # Need to separate them out.
    gene_names = []
    for x in gene_name_list:
        x = x.split(",")
        gene_names.extend(x)
    return gene_names


def _parse_color(color_str):
    # color_str is <R>,<G>,<B> where each number is an integer from
    # 0-255.  Return tuple of (<R>, <G>, <B>).
    x = color_str.split(",")
    assert len(x) == 3, "color should be <R>,<G>,<B>"
    x = [int(x) for x in x]
    for i in range(len(x)):
        assert x[i] >= 0 and x[i] < 256, "color should be 0-255"
    return tuple(x)
    

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

    usage = "usage: %prog [options] matrix_file png_file\n" + \
            "  Values should range from -1 to 1."
    parser = OptionParser(usage=usage, version="%prog 02")

    parser.add_option(
        "--libpath", action="append", default=[],
        help="Add to the Python library search path.")
    # XXX need way to save processed matrix to a file

    # For discrete clusters, like in k-means.
    parser.add_option("--gene_cluster_file", help="kgg file")
    parser.add_option("--array_cluster_file", help="kag file")
    # For trees.
    parser.add_option("--gene_tree_file", help="gtr file")
    parser.add_option("--array_tree_file", help="atr file")
    # A way to associate trees with specific cluster for the branches.
    # (Not required for trees).
    parser.add_option("--gene_tree_cluster_file", help="gtc file")
    parser.add_option("--array_tree_cluster_file", help="atc file")
    

    group = OptionGroup(parser, "Heatmap")
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
        "--no_autoscale", dest="autoscale", action="store_false",
        default=True, help="Disable autoscaling.")
    group.add_option(
        "--color", dest="color_scheme", type="choice", default="bild",
        choices=["red", "white", "red-green", "blue-yellow", "matlab", "bild",
                 "genespring", "yahoo"],
        help="Choose the color scheme to use: red, white, red-green, "
        "blue-yellow, matlab, bild (default), genespring, or yahoo.")
    group.add_option(
        "--black0", dest="black0", action="store_true", default=False,
        help="Color 0 values black (no matter the color scheme).")
    group.add_option(
        "--inverse", dest="inverse", action="store_true", default=False,
        help="Flip the colors for the heatmap.")
    group.add_option(
        "--scale_border", default=1.0, type="float",
        help="Scale the border thicker or thinner.")
    group.add_option(
        "--border_color",
        help="Specify the color of the border.  "
        "Format: <R>,<G>,<B>  (e.g. 128,128,128)")
    group.add_option(
        "--grid", action="store_true", default=False,
        help="Add a grid around the boxes in the heatmap.")
    group.add_option(
        "--grid_color",
        help="Specify the color of the grid.  "
        "Format: <R>,<G>,<B>  (e.g. 128,128,128)")
    parser.add_option_group(group)
    
    group = OptionGroup(parser, "Labels")
    group.add_option(
        "--gl", "--label_genes", dest="label_genes", action="store_true",
        default=False, help="Label the genes on the plot.")
    group.add_option(
        "--al", "--label_arrays", dest="label_arrays", action="store_true",
        default=False, help="Label the arrays on the plot.")
    parser.add_option_group(group)
    
    group = OptionGroup(parser, "Dendrogram")
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

    group = OptionGroup(parser, "Colorbar")
    group.add_option(
        "--colorbar", dest="colorbar", default=False, action="store_true",
        help="Add a colorbar to the plot.")
    parser.add_option_group(group)

    # Parse the input arguments.
    options, args = parser.parse_args()
    if not args or len(args) != 2:
        parser.error("Please specify an infile and an outfile.")
    infile, outfile = args
    if not os.path.exists(infile):
        parser.error("I could not find file %s." % infile)
        
    if options.libpath:
        sys.path = options.libpath + sys.path

    # Check the options.
    # Not completely implemented yet.
    assert options.scale_border > 0 and options.scale_border < 5.0
    border_color = 0, 0, 0
    # (0, 0, 0) is too dark for small box sizes.  100 looks too washed
    # out.  50-75 is about right.
    grid_color = 75, 75, 75
    if options.border_color:
        border_color = _parse_color(options.border_color)
    if options.grid_color:
        grid_color = _parse_color(options.grid_color)

    # Choose a plotting library.
    plotlib = __import__(
        "genomicode.pilplot", globals(), locals(), ["pilplot"])

    MATRIX, cluster_data = read_data_set(
        infile, options.gene_cluster_file, options.array_cluster_file,
        options.gene_tree_file, options.array_tree_file,
        options.gene_tree_cluster_file, options.array_tree_cluster_file)
    x = process_data_set(
        MATRIX, options.scale, options.gain, options.autoscale)
    MATRIX, signal_0, signal_1 = x
    layout = make_layout(
        MATRIX, cluster_data, plotlib,
        # Heatmap
        options.width, options.height, options.scale_border, options.grid,
        options.color_scheme, options.inverse,
        signal_0, signal_1,  options.black0,
        # Labels
        options.label_genes, options.label_arrays,
        # Dendrograms
        options.gene_tree_scale, options.gene_tree_thickness,
        options.array_tree_scale, options.array_tree_thickness,
        # Colorbar
        options.colorbar,
        )

    megapixels = layout.heatmap.width() * layout.heatmap.height() / 1024 / 1024
    assert megapixels <= MAX_MEGAPIXELS, "%dx%d plot too big [%d:%d]." % (
        layout.heatmap.width(), layout.heatmap.height(),
        options.width, options.height)
    
    coords = calc_coords_for_layout(layout)
    plot(
        outfile, MATRIX, cluster_data, plotlib, layout, coords,
        border_color, grid_color)


if __name__ == '__main__':
    main()
    #import profile; profile.run("main()")
