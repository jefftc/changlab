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


def main():
    import os
    import argparse
    
    from genomicode import jmath
    from genomicode import AnnotationMatrix

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("datafile", help="Tab-delimited data file.")
    parser.add_argument("header", help="Which column contains data to plot.")
    parser.add_argument(
        "plot_file", help="Name of image file, e.g. outfile.png.  "
        "Will generate PNG format by default.  If this file name ends with "
        ".pdf, will generate a PDF file instead.")

    group = parser.add_argument_group(title="Calculations")
    group.add_argument(
        "--breaks_seq",
        help="Set the breakpoints.  Format: <start>,<stop>,<skip>.")
    
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
    group.add_argument(
        "--xlabel_size", default=1.0, type=float,
        help="Scale the size of the labels on X-axis.  Default 1.0.")
    

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
    assert args.xlabel_size > 0 and args.xlabel_size < 10
        
    height = args.height or 2400
    width = args.width or 3200

    MATRIX = AnnotationMatrix.read(args.datafile, False)
    assert MATRIX.num_headers() and MATRIX.num_annots(), "Empty matrix."
    assert args.header in MATRIX.headers, "header not found: %s" % args.header

    # Pull out the values for the histogram.
    x = MATRIX.get_annots(args.header)
    values = map(float, x)

    value_min = value_max = None

    # Start R and set up the environment.
    R = jmath.start_R()

    main = jmath.R_var("NA")
    if args.title:
        main = args.title
    sub = ""
    xlab = ""
    if args.xlab:
        xlab = args.xlab
    ylab = "Frequency"
    
    breaks = "Sturges"
    if args.breaks_seq:
        breaks = _parse_breaks_seq(args.breaks_seq)
        value_min, value_max = min(breaks), max(breaks)
        jmath.R_equals(breaks, "breaks")
        breaks = jmath.R_var("breaks")

    if value_min is not None:
        values = [x for x in values if x >= value_min]
    if value_max is not None:
        values = [x for x in values if x < value_max]

    lwd = 2
    cex_lab = 1.5
    cex_main = 2.0
    cex_sub = 1.5

    assert values
    jmath.R_equals(values, "X")

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
        "hist", jmath.R_var("X"), breaks=breaks, main=main, xlab="", ylab="",
        axes=jmath.R_var("FALSE"))
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


if __name__ == '__main__':
    main()
