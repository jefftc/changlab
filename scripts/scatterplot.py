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

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("datafile", help="Tab-delimited data file.")
    parser.add_argument("x_header", help="Which column for X values.")
    parser.add_argument("y_header", help="Which column for Y values.")
    parser.add_argument(
        "plot_file", help="Name of image file, e.g. outfile.png.  "
        "Will generate PNG format by default.  If this file name ends with "
        ".pdf, will generate a PDF file instead.")

    group = parser.add_argument_group(title="Plot Labels")
    group.add_argument("--title", help="Put a title on the plot.")
    group.add_argument("--xlab", help="Label the X-axis.")

    group = parser.add_argument_group(title="Margins and Sizes")
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
    #group.add_argument(
    #    "--xlabel_size", default=1.0, type=float,
    #    help="Scale the size of the labels on X-axis.  Default 1.0.")
    

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
        
    height = args.height or 2400
    width = args.width or 3200

    MATRIX = AnnotationMatrix.read(args.datafile, False)
    assert MATRIX.num_headers() and MATRIX.num_annots(), "Empty matrix."
    assert args.x_header in MATRIX.headers, \
           "header not found: %s" % args.x_header
    assert args.y_header in MATRIX.headers, \
           "header not found: %s" % args.y_header

    # Pull out the values for the plot.
    x1 = MATRIX.get_annots(args.x_header)
    x2 = MATRIX.get_annots(args.y_header)
    x_values = map(float, x1)
    y_values = map(float, x2)

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
    
    lwd = 2
    cex_lab = 1.5
    cex_main = 2.0
    cex_sub = 1.5

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
        main=main, xlab="", ylab="",
        axes=jmath.R_var("FALSE"), RETVAL="x")
    # Make plot area solid white.
    #jmath.R('usr <- par("usr")')
    #jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    #jmath.R_fn(
    #    "hist", jmath.R_var("X"), plot=jmath.R_var("FALSE"),
    #    main=main, xlab="", ylab="", axes=jmath.R_var("FALSE"),
    #    add=jmath.R_var("TRUE"))
    
    #jmath.R_fn("box", lwd=lwd)
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
