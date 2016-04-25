#!/usr/bin/env python

def _parse_breaks_seq(s):
    # Format: <start>,<stop>,<skip>
    from genomicode import jmath
    x = s.split(",")
    assert len(x) == 3, "Format: <start>,<stop>,<skip>"

    start, stop, skip = x
    start = start.replace("n", "-")
    stop = stop.replace("n", "-")
    start, stop, skip = float(start), float(stop), float(skip)
    assert stop > start
    assert skip > 0

    tol = min(1E-5, skip/100.0)
    breaks = []
    i = 0
    while True:
        s = start + i*skip
        if s >= stop-tol:
            break
        breaks.append(s)
        i += 1
    return breaks


def _parse_cluster(options_cluster, indexes_include_headers, MATRIX):
    # Return a vector of clusters, where each cluster is an integer
    # from 0 to K-1.  K is the total number of clusters.  The length
    # of the vector should be the same as the number of annotations in
    # the matrix.
    from genomicode import parselib

    index2cluster = {}
    for clust_i, s in enumerate(options_cluster):
        ranges = parselib.parse_ranges(s)
        for s, e in ranges:
            for i in range(s-1, e):
                # Convert from 1-based to 0-based indexes.
                i -= 1
                if indexes_include_headers:
                    i -= 1  # assume 1 header
                    #i -= len(MATRIX._row_names)
                assert i < MATRIX.num_annots(), \
                       "Index %d out of range" % i
                assert i not in index2cluster, \
                       "Index %d in multiple clusters" % i
                index2cluster[i] = clust_i
    #cluster = [len(options_cluster)] * MATRIX.ncol()
    cluster = [None] * MATRIX.num_annots()
    for i, g in index2cluster.iteritems():
        cluster[i] = g
    return cluster


def write_prism_file(filename, hist):
    # hist is R list from hist function.
    from genomicode import jmath

    # XY plot in Prism.

    # Get "breaks" out of histogram return value.
    breaks = [x for x in hist.rx2("breaks")]
    breaks = breaks[:-1]
    counts = [x for x in hist.rx2("counts")]
    density = [x for x in hist.rx2("density")]
    mids = [x for x in hist.rx2("mids")]
    assert len(breaks) == len(counts)
    assert len(breaks) == len(density)
    assert len(breaks) == len(mids)

    header = ["Mids", "Left", "Counts", "Density"]
    x = [mids, breaks, counts, density]
    x = jmath.transpose(x)
    x = [header] + x
    handle = open(filename, 'w')
    for x in x:
        print >>handle, "\t".join(map(str, x))


def main():
    import os
    import argparse
    
    from genomicode import jmath
    from genomicode import AnnotationMatrix
    from genomicode import colorlib
    from genomicode import pcalib

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("datafile", help="Tab-delimited data file.")
    parser.add_argument("x_header", help="Which column for X values.")
    parser.add_argument("y_header", help="Which column for Y values.")
    parser.add_argument(
        "plot_file", help="Name of image file, e.g. outfile.png.  "
        "Will generate PNG format by default.  If this file name ends with "
        ".pdf, will generate a PDF file instead.")
    

    group = parser.add_argument_group(title="General Appearance")
    group.add_argument(
        "--no_box", action="store_true",
        help="Turn off the box around the plot.")
    group.add_argument(
        "--height", type=int, help="Height (in pixels) of the plot.")
    group.add_argument(
        "--width", type=int, help="Width (in pixels) of the plot.")
    group.add_argument(
        "--mar_left", default=1.0, type=float,
        help="Scale margin at left of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--mar_bottom", default=1.0, type=float,
        help="Scale margin at bottom of plot.  Default 1.0.")
    group.add_argument(
        "--scale_points", default=1.0, type=float,
        help="Scale the size of the points.  Default 1.0")
    #group.add_argument(
    #    "--xlabel_size", default=1.0, type=float,
    #    help="Scale the size of the labels on X-axis.  Default 1.0.")
    group.add_argument(
        "--log_x", action="store_true",
        help="Plot the X-axis on a log scale.")
    group.add_argument(
        "--log_y", action="store_true",
        help="Plot the Y-axis on a log scale.")

    group = parser.add_argument_group(title="Plot Labels")
    group.add_argument("--title", help="Put a title on the plot.")
    group.add_argument("--xlab", help="Label the X-axis.")
    group.add_argument("--ylab", help="Label the Y-axis.")
    group.add_argument(
        "--add_regression", action="store_true",
        help="Put a regression line on the plot.")

    group = parser.add_argument_group(title="Point Labels")
    group.add_argument(
        "--label_header",
        help="Label each point with the values in this column.")
    group.add_argument(
        "--label_size", type=float, 
        help="Scale the size of the labels by this value.")
    group.add_argument(
        "--label_pos", default="top",
        choices=["top", "bottom", "left", "right"],
        help="Where to label the points.")
    
    group = parser.add_argument_group(title="Colors")
    group.add_argument(
        "-c", "--cluster", action="append",
        help="Group samples into a cluster (e.g. -c 1-5); 1-based.")
    group.add_argument(
        "--indexes_include_headers", "--iih", action="store_true",
        help="If not given (default), then index 1 is the first row "
        "with data.  If given, then index 1 is the very first row "
        "in the file, including the headers.")
    group.add_argument(
        "--default_color",
        help="Default color of points.  Format: #000000.")
    

    # Parse the input arguments.
    args = parser.parse_args()
    if not os.path.exists(args.datafile):
        parser.error("File not found: %s" % args.datafile)
    if args.width is not None:
        assert args.width > 10, "too small"
        assert args.width < 4096*16, "width too big"
    if args.height is not None:
        assert args.height > 10, "too small"
        assert args.height < 4096*16, "height too big"
    assert args.mar_bottom > 0 and args.mar_bottom < 10
    assert args.mar_left > 0 and args.mar_left < 10
    #assert args.xlabel_size > 0 and args.xlabel_size < 10

    MATRIX = AnnotationMatrix.read(args.datafile, False)
    assert MATRIX.num_headers() and MATRIX.num_annots(), "Empty matrix."
    assert args.x_header in MATRIX.headers, \
           "header not found: %s" % args.x_header
    assert args.y_header in MATRIX.headers, \
           "header not found: %s" % args.y_header
    if args.label_header:
        assert args.label_header in MATRIX.headers, \
               "header not found: %s" % args.label_header
    if args.label_size is not None:
        assert args.label_size > 0 and args.label_size <= 20

    cluster = None
    if args.cluster:
        cluster = _parse_cluster(
            args.cluster, args.indexes_include_headers, MATRIX)
    
    height = args.height or 2400
    width = args.width or 3200

    # Pull out the values for the plot.
    x1 = MATRIX[args.x_header]
    x2 = MATRIX[args.y_header]
    x_values = map(float, x1)
    y_values = map(float, x2)
    assert len(x_values) == len(y_values)

    col = "#000000"
    if args.default_color:
        col = args.default_color
    assert len(col) == 7
    assert col[0] == "#"
    if cluster is not None:
        col_rgb = pcalib.choose_colors(cluster)
        col = [col] * len(col_rgb)
        for i in range(len(col_rgb)):
            if col_rgb[i] is None:
                continue
            col[i] = colorlib.rgb2hex(col_rgb[i], prefix="#")
        assert len(col) == len(x_values)

        
    # Start R and set up the environment.
    R = jmath.start_R()

    main = jmath.R_var("NA")
    if args.title:
        main = args.title
    sub = ""
    xlab = args.x_header
    if args.xlab:
        xlab = args.xlab
    ylab = args.y_header
    if args.xlab:
        ylab = args.ylab
    
    
    lwd = 2
    cex = 1 * args.scale_points
    cex_lab = 1.5
    cex_main = 2.0
    cex_sub = 1.0

    plot_log = ""
    if args.log_x:
        plot_log += "x"
    if args.log_y:
        plot_log += "y"

    assert x_values
    assert y_values
    jmath.R_equals(x_values, "X")
    jmath.R_equals(y_values, "Y")

    bm_type = "png16m"
    if args.plot_file.lower().endswith(".pdf"):
        bm_type = "pdfwrite"
    jmath.R_fn(
        "bitmap", args.plot_file, type=bm_type, 
        height=height, width=width, units="px", res=300)
    
    # Set the margins.
    x = 5*1.2*args.mar_bottom, 4*1.2*args.mar_left, 4, 2
    mar = [x+0.1 for x in x]
    jmath.R_fn("par", mar=mar, RETVAL="op")

    jmath.R_fn(
        "plot", jmath.R_var("X"), jmath.R_var("Y"),
        main="", xlab="", ylab="", pch=19, cex=cex,
        log=plot_log, col=col,
        axes=jmath.R_var("FALSE"), RETVAL="x")
    # Make plot area solid white.
    #jmath.R('usr <- par("usr")')
    #jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    #jmath.R_fn(
    #    "hist", jmath.R_var("X"), plot=jmath.R_var("FALSE"),
    #    main=main, xlab="", ylab="", axes=jmath.R_var("FALSE"),
    #    add=jmath.R_var("TRUE"))

    if args.label_header:
        cex = 1
        if args.label_size is not None:
            cex = args.label_size
        pos2specifier = {
            "top" : 3,
            "bottom" : 1,
            "left" : 2,
            "right" : 4,
            }
        pos = pos2specifier[args.label_pos]
        point_labels = MATRIX[args.label_header]
        jmath.R_fn(
            "text", jmath.R_var("X"), jmath.R_var("Y"),
            labels=point_labels, cex=cex, pos=pos)
    

    # Calculate correlation, and other statistics.
    r = jmath.R("cor(X, Y)")
    p_value = jmath.R("cor.test(X, Y)$p.value")
    r = r[0]
    p_value = p_value[0]

    # Add a regression line.
    if args.add_regression:
        jmath.R("fit <- lm(Y ~ X)")
        coef = jmath.R("fit$coefficients")
        assert len(coef) == 2
        b, m = coef
        x1 = min(x_values)
        y1 = x1*m + b
        x2 = max(x_values)
        y2 = x2*m + b
        jmath.R_fn("lines", [x1, x2], [y1, y2], lwd=3, lty=2, col="#C63F31")
        sub = "R=%.2f (p=%.2g)" % (r, p_value)
        header = "X", "Y", "R", "p"
        print "\t".join(header)
        x = xlab, ylab, r, p_value
        print "\t".join(map(str, x))

    if not args.no_box:
        jmath.R_fn("box", lwd=lwd)
    jmath.R_fn("axis", 1, lwd=lwd, **{ "cex.axis" : 1.5 })
    jmath.R_fn("axis", 2, lwd=lwd, **{ "cex.axis" : 1.5 })
    jmath.R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
        **{ "cex.lab" : cex_lab, "cex.main" : cex_main, "cex.sub" : cex_sub })
    R("par(op)")
    jmath.R_fn("dev.off")

    #if args.prism_file:
    #    write_prism_file(args.prism_file, R["x"])
        
if __name__ == '__main__':
    main()
