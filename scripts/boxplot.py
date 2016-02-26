#!/usr/bin/env python

import os


def plot_boxplot(
    filename, group_names, group2values, height=None, width=None, cluster=None,
    title="", subtitle="", sub="", xlab="", ylab="",
    subtitle_size=1.0, subtitle_line=0.5, subtitle_col="#000000",
    xlabel_size=1.0, xlabel_off=False,
    mar_bottom=1.0, mar_left=1.0, mar_top=1.0):
    # group_names is a list of the names for each group.
    # group2values is a dictionary of group_name -> list of values.
    # Also, can be matrix (values x groups).
    # subtitle goes under title.  sub goes under plot.
    from genomicode import config
    from genomicode import jmath
    from genomicode import colorlib
    from genomicode import pcalib

    # Start R and set up the environment.
    R = jmath.start_R()
    path = config.changlab_Rlib
    plotlib = os.path.join(path, "plotlib.R")
    assert os.path.exists(plotlib), "I cannot find: %s" % plotlib
    jmath.R_fn("source", plotlib)

    main = jmath.R_var("NA")
    if title:
        main = title
    sub = sub
    xlab = xlab
    ylab = ylab
    xlabel = group_names
    if xlabel_off:
        xlabel = jmath.R_var("FALSE")

    col = jmath.R_var("NULL")
    if cluster is not None:
        x = pcalib.choose_colors(cluster)
        x = [colorlib.rgb2hex(x) for x in x]
        x = [x.replace("0x", "#") for x in x]
        col = x

    lwd = 2
    las = 3   # vertical labels
    at = jmath.R_var("NULL")
    if xlabel != jmath.R_var("FALSE"):
        at = range(1, len(xlabel)+1)
    cex_labels = 1.25*xlabel_size
    #cex_legend = 1
    cex_xlab = 1.5
    cex_ylab = 2.0
    cex_sub = 1.5

    if type(group2values) is type([]):
        # Is matrix.  Should do more checking here.
        jmath.R_equals(group2values, "X")
    else:
        R("X <- list()")
        for i, n in enumerate(group_names):
            x = group2values.get(n, [])
            x = [x for x in x if x is not None]
            jmath.R_equals(x, "s")
            R("X[[%d]] <- s" % (i+1))

    #try:
    #    #jmath.R_equals(MATRIX._X, "X")
    #    jmath.R_equals(X, "X")
    #except ValueError, x:
    #    # Not needed anymore.  Missing values are now implemented in jmath.
    #    ## Look for missing values.
    #    #for i in range(len(MATRIX._X)):
    #    #    assert None not in MATRIX._X[i], \
    #    #           "Missing values in row %d (0-based)." % i
    #    ## Cannot diagnose error.  Raise the original exception.
    #    raise
    
    jmath.R_equals(xlabel, "labels")
    jmath.R_equals(at, "at")

    bm_type = "png16m"
    if filename.lower().endswith(".pdf"):
        bm_type = "pdfwrite"
    jmath.R_fn(
        "bitmap", filename, type=bm_type,
        height=height, width=width, units="px", res=300)
    
    # Set the margins.
    # default is 5.1, 4.1, 4.1, 2.1
    label_adjust = 1.0
    if xlabel == jmath.R_var("FALSE"):
        label_adjust = 0.2
    x = 5*2.0*mar_bottom*label_adjust, 4*1.2*mar_left, 4*mar_top, 2
    mar = [x+0.1 for x in x]
    jmath.R_fn("par", mar=mar, RETVAL="op")
        
    jmath.R_fn(
        "boxplot", jmath.R_var("X"), col=col, main="", xlab="", ylab="",
        axes=jmath.R_var("FALSE"), pch=19, cex=1, ylim=jmath.R_var("NULL"))
    # Make plot area solid white.
    jmath.R('usr <- par("usr")')
    jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    jmath.R_fn(
        "boxplot", jmath.R_var("X"), col=col, main="", xlab="", ylab="",
        axes=jmath.R_var("FALSE"), pch=19, cex=1, ylim=jmath.R_var("NULL"),
        add=jmath.R_var("TRUE"))
    
    jmath.R_fn("box", lwd=lwd)
    jmath.R_fn(
        "axis", 1, lwd=lwd, labels=jmath.R_var("labels"),
        at=jmath.R_var("at"), las=las, **{ "cex.axis" : cex_labels })
    jmath.R_fn(
        "axis", 2, lwd=lwd, **{ "cex.axis" : 1.5 })
    jmath.R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
        **{ "cex.lab" : cex_xlab, "cex.main" : 2.0, "cex.sub" : cex_sub,
            "col.sub" : "#A60400" })
    if subtitle:
        jmath.R_fn(
            "mtext", subtitle, cex=1.0*subtitle_size,
            line=subtitle_line, col=subtitle_col)
    R("par(op)")
    jmath.R_fn("dev.off")


def main():
    import argparse
    
    import arrayio

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("expression_file", help="Gene expression file.")
    parser.add_argument(
        "plot_file", help="Name of image file, e.g. outfile.png.  "
        "Will generate PNG format by default.  If this file name ends with "
        ".pdf, will generate a PDF file instead.")

    parser.add_argument(
        "--num_header_cols", type=int,
        help="This number of columns are headers.  If not given, will guess.")
    parser.add_argument(
        "-v", "--verbose", default=False, action="store_true",
        help="")

    # Plot labels.
    group = parser.add_argument_group(title="Plot labels")
    group.add_argument(
        "--title", default=None, help="Put a title on the plot.")

    # Margins and sizes.
    group = parser.add_argument_group(title="Margins and sizes")
    group.add_argument(
        "--height", default=None, type=int,
        help="Height (in pixels) of the plot.")
    group.add_argument(
        "--width", default=None, type=int,
        help="Width (in pixels) of the plot.")
    group.add_argument(
        "--mar_left", default=1.0, type=float,
        help="Scale margin at left of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--mar_bottom", default=1.0, type=float,
        help="Scale margin at bottom of plot.  Default 1.0.")
    group.add_argument(
        "--xlabel_size", default=1.0, type=float,
        help="Scale the size of the labels on X-axis.  Default 1.0.")
    group.add_argument(
        "--xlabel_off", default=False, action="store_true",
        help="Turn off the X labels.")

    group = parser.add_argument_group(title="Coloring")
    group.add_argument(
        "-c", "--cluster", default=[], action="append",
        help="Group samples into a cluster (e.g. -c 1-5); 1-based.")
    group.add_argument(
        "--indexes_include_headers", "--iih", action="store_true",
        help="If not given (default), then index 1 is the first column "
        "with data.  If given, then index 1 is the very first column "
        "in the file, including the headers.")

    # Parse the input arguments.
    args = parser.parse_args()
    if not os.path.exists(args.expression_file):
        parser.error("I could not find file %s." % args.expression_file)
    if args.num_header_cols is not None:
        assert args.num_header_cols > 0 and args.num_header_cols < 100
    if args.width is not None:
        assert args.width > 10, "too small"
        assert args.width < 4096*16, "width too big"
    if args.height is not None:
        assert args.height > 10, "too small"
        assert args.height < 4096*16, "height too big"
    assert args.mar_bottom > 0 and args.mar_bottom < 10
    assert args.mar_left > 0 and args.mar_left < 10
    assert args.xlabel_size > 0 and args.xlabel_size < 10
    height = args.height or 1600
    width = args.width or 1600

    MATRIX = arrayio.read(args.expression_file, hcols=args.num_header_cols)
    assert MATRIX.nrow() and MATRIX.ncol(), "Empty matrix."

    cluster = None
    if args.cluster:
        import pcaplot
        cluster = pcaplot._parse_cluster(
            args.cluster, args.indexes_include_headers, MATRIX)

    group_names = MATRIX.col_names(arrayio.COL_ID)
    #group2values = {}
    #for i in range(len(group_names)):
    #    values = [x[i] for x in MATRIX._X]
    #    group2values[group_names[i]] = values
    plot_boxplot(
        args.plot_file, group_names, MATRIX._X,
        height=height, width=width, cluster=cluster,
        title=args.title, xlab="", ylab="Gene Expression", 
        xlabel_size=args.xlabel_size, xlabel_off=args.xlabel_off,
        mar_bottom=args.mar_bottom, mar_left=args.mar_left)


if __name__ == '__main__':
    main()
