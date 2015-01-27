#!/usr/bin/env python

import os

def main():
    import argparse
    
    import arrayio
    from genomicode import config
    from genomicode import jmath

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("expression_file", help="Gene expression file.")
    parser.add_argument(
        "plot_file", help="Name of image file, e.g. outfile.png.  "
        "Will generate PNG format by default.  If this file name ends with "
        ".pdf, will generate a PDF file instead.")

    parser.add_argument(
        "-v", "--verbose", default=False, action="store_true",
        help="")

    # Plot labels.
    parser.add_argument(
        "--title", default=None, help="Put a title on the plot.")

    # Margins and sizes.
    parser.add_argument(
        "--height", default=None, type=int,
        help="Height (in pixels) of the plot.")
    parser.add_argument(
        "--width", default=None, type=int,
        help="Width (in pixels) of the plot.")
    parser.add_argument(
        "--mar_left", default=1.0, type=float,
        help="Scale margin at left of plot.  Default 1.0 (no scaling).")
    parser.add_argument(
        "--mar_bottom", default=1.0, type=float,
        help="Scale margin at bottom of plot.  Default 1.0.")
    parser.add_argument(
        "--xlabel_size", default=1.0, type=float,
        help="Scale the size of the labels on X-axis.  Default 1.0.")
    parser.add_argument(
        "--xlabel_off", default=False, action="store_true",
        help="Turn off the X labels.")
    

    # Parse the input arguments.
    args = parser.parse_args()
    if not os.path.exists(args.expression_file):
        parser.error("I could not find file %s." % args.expression_file)
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

    MATRIX = arrayio.read(args.expression_file)
    assert MATRIX.nrow() and MATRIX.ncol(), "Empty matrix."

    # Start R and set up the environment.
    R = jmath.start_R()
    path = config.changlab_Rlib
    plotlib = os.path.join(path, "plotlib.R")
    assert os.path.exists(plotlib), "I cannot find: %s" % plotlib
    jmath.R_fn("source", plotlib)

    main = jmath.R_var("NA")
    if args.title:
        main = args.title
    sub = ""
    xlab = ""
    ylab = "Gene Expression"
    labels = jmath.R_var("FALSE")
    col = jmath.R_var("NULL")

    if not args.xlabel_off:
        labels = MATRIX.col_names(arrayio.COL_ID)

    lwd = 2
    las = 3   # vertical labels
    at = jmath.R_var("NULL")
    if labels != jmath.R_var("FALSE"):
        at = range(1, len(labels)+1)
    cex_labels = 1.25*args.xlabel_size
    cex_legend = 1
    cex_lab = 1.5
    cex_sub = 1.5

    try:
        jmath.R_equals(MATRIX._X, "X")
    except ValueError, x:
        # Look for missing values.
        for i in range(len(MATRIX._X)):
            assert None not in MATRIX._X[i], \
                   "Missing values in row %d (0-based)." % i
        # Cannot diagnose error.  Raise the original exception.
        raise
    
    jmath.R_equals(labels, "labels")
    jmath.R_equals(at, "at")

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
        "boxplot", jmath.R_var("X"), col=col, main=main, xlab="", ylab="",
        axes=jmath.R_var("FALSE"), pch=19, cex=1, ylim=jmath.R_var("NULL"))
    # Make plot area solid white.
    jmath.R('usr <- par("usr")')
    jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    jmath.R_fn(
        "boxplot", jmath.R_var("X"), col=col, main=main, xlab="", ylab="",
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
        **{ "cex.lab" : cex_lab, "cex.main" : 2.0, "cex.sub" : cex_sub })
    R("par(op)")
    jmath.R_fn("dev.off")


if __name__ == '__main__':
    main()
