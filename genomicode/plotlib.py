"""

Functions:
plot_heatmap

find_tall_heatmap_size
find_wide_heatmap_size

"""

def plot_heatmap(
    infile, outfile, xpix, ypix, color=None, 
    gene_label=False, cluster_genes=False,
    gene_center=None, gene_normalize=None,
    array_label=False, cluster_arrays=False,
    scale=None, gain=None, no_autoscale=False, 
    python=None, arrayplot=None, cluster=None, libpath=None):
    import os
    import sys
    import subprocess

    # If arrayplot is not supplied, then use the default arrayplot.py.
    # This may not be in the current path, so be sure not to include
    # python.
    if not arrayplot:
        python = None
        arrayplot = "arrayplot.py"
    color = color or "bild"

    cmd = [
        python,
        arrayplot,
        "-x %d" % xpix,
        "-y %d" % ypix,
        "-o %s" % outfile,
        "--color=%s" % color,
        ]
    if gene_label:
        cmd.append("--gl")
    if cluster_genes:
        cmd.append("-g")
    if gene_center:
        assert gene_center in ["mean", "median"]
        cmd.append("--gc=%s" % gene_center)
    if gene_normalize:
        assert gene_normalize in ["ss", "var"]
        cmd.append("--gn=%s" % gene_normalize)
    if array_label:
        cmd.append("--al")
    if cluster_arrays:
        cmd.append("-a")
    if scale:
        cmd.append("--scale=%s" % scale)
    if gain:
        cmd.append("--gain=%s" % gain)
    if no_autoscale:
        cmd.append("--no_autoscale")
    if cluster is not None:
        cmd.append("--cluster_app=%s" % cluster)
    if libpath:
        for path in libpath:
            cmd.append("--libpath=%s" % path)
    cmd.append("'%s'" % infile)

    # If python is None, then ignore it.
    cmd = [x for x in cmd if x is not None]

    cmd = " ".join(cmd)
    #print cmd
    #w, r = os.popen4(cmd)
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout
    w.close()
    output = r.read()

    # Make sure the signature was generated correctly.  An error could
    # mean that arrayplot.py or cluster is missing.
    if not os.path.exists(outfile):
        print >>sys.stderr, output
        raise AssertionError, "Failed to make dataset."

    return output

def find_tall_heatmap_size(
    nrow, ncol, min_box_height=None, min_box_width=None,
    max_box_height=None, max_box_width=None,
    max_total_height=None, max_total_width=None, height_width_ratio=None):
    # Return tuple of the pixels for each box as (xpix, ypix).
    # Minimum sizes take precedence over maximum sizes.

    # Set some defaults.
    min_box_height = min_box_height or 1.0
    min_box_width = min_box_width or 1.0
    max_total_height = max_total_height or 600
    max_total_width = max_total_width or 800
    max_box_height = max_box_height or max(min_box_height,
                                           float(max_total_height)/nrow)
    max_box_width = max_box_width or max(min_box_width,
                                         float(max_total_width)/ncol)
    height_width_ratio = height_width_ratio or 1.618  # HEIGHT/WIDTH

    min_box_height = float(min_box_height)
    min_box_width = float(min_box_width)
    max_box_height = float(max_box_height)
    max_box_width = float(max_box_width)
    max_total_height = float(max_total_height)
    max_total_width = float(max_total_width)
    height_width_ratio = float(height_width_ratio)
    
    # Start with the minimum size.
    ypix = min_box_height
    # Set the width to match the height.
    total_y = ypix * nrow
    total_x = total_y / height_width_ratio
    #print total_x, total_y
    xpix = float(total_x) / ncol
    #print "HERE1", xpix, ypix
    # If the width is too small, then use this to set the minimum.
    if xpix < min_box_width:
        xpix = min_box_width
        total_x = xpix * ncol
        total_y = total_x * height_width_ratio
        ypix = float(total_y) / nrow
        assert ypix >= min_box_height
        #print "HERE3", xpix, ypix
        
    # Increase size up to the maximum allowed.
    max_xpix = min(max_box_width, float(max_total_width) / ncol)
    max_ypix = min(max_box_height, float(max_total_height) / nrow)
    #print "HERE2", max_xpix, max_ypix, xpix, ypix
    x_ratio = max_xpix / xpix
    y_ratio = max_ypix / ypix
    ratio = min(y_ratio, x_ratio)
    if ratio > 1.0:
        xpix, ypix = xpix*ratio, ypix*ratio
    #print xpix, ypix

    xpix, ypix = int(xpix), int(ypix)
    height = nrow * ypix
    width = ncol * xpix
    #print "DIM %dx%d (%d*%d, %d*%d)" % (height, width, nrow, ypix, ncol, xpix)
    return xpix, ypix

def find_wide_heatmap_size(
    nrow, ncol, min_box_height=None, min_box_width=None,
    max_box_height=None, max_box_width=None,
    max_total_height=None, max_total_width=None, height_width_ratio=None):
    inv_height_width_ratio = height_width_ratio
    if inv_height_width_ratio is not None:
        inv_height_width_ratio = 1.0 / inv_height_width_ratio
    x = find_tall_heatmap_size(
        ncol, nrow,
        min_box_height=min_box_width, min_box_width=min_box_height,
        max_box_height=max_box_width, max_box_width=max_box_height,
        max_total_height=max_total_width, max_total_width=max_total_height,
        height_width_ratio=inv_height_width_ratio)
    xpix, ypix = x
    xpix, ypix = ypix, xpix
    height = nrow * ypix
    width = ncol * xpix
    #print "DIM %dx%d (%d*%d, %d*%d)" % (height, width, nrow, ypix, ncol, xpix)
    return xpix, ypix

def test_find_heatmap_size():
    print find_tall_heatmap_size(10, 10)  # 12, 20
    print find_tall_heatmap_size(5, 10)   # 6, 20
    print find_tall_heatmap_size(10, 10, min_box_width=100)  # 100, 161
    print find_tall_heatmap_size(5, 10, min_box_width=100)   # 100, 323

    print find_wide_heatmap_size(10, 10)  # 20, 12
    print find_wide_heatmap_size(5, 10)   # 16, 20

if __name__ == '__main__':
    test_find_heatmap_size()
