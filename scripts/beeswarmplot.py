#!/usr/bin/env python

def read_prism_file(filename):
    # <header1>  <header2>
    # <values>   <values>
    # ...
    from genomicode import genesetlib

    MATRIX = []
    for x in genesetlib.read_gmx(filename, allow_duplicates=True):
        name, description, genes = x
        x = [description] + genes
        x = [float(x) for x in x]
        values = x
        x = name, values
        MATRIX.append(x)
    return MATRIX
        
        
def main():
    import os
    import argparse
    
    from genomicode import jmath
    #from genomicode import AnnotationMatrix
    #from genomicode import colorlib
    #from genomicode import pcalib
    from genomicode import hashlib

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "datafile", help="Tab-delimited text file in Prism format.  "
        "Each column is a series.  First row is header.")
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

    group = parser.add_argument_group(title="Plot Labels")
    group.add_argument("--title", help="Put a title on the plot.")
    group.add_argument("--xlab", help="Label the X-axis.")
    group.add_argument("--ylab", help="Label the Y-axis.")

    group = parser.add_argument_group(title="Legend")
    group.add_argument(
        "--add_legend", action="store_true", help="Add a legend to the plot.")
    group.add_argument(
        "--legend_inset", type=float, default=0.05,
        help="")
    LEGEND_LOCATIONS = [
        "bottomright", "bottom", "bottomleft", "left", "topleft", "top",
        "topright", "right", "center",]
    group.add_argument(
        "--legend_loc", choices=LEGEND_LOCATIONS, 
        help="Where to draw the legend.")
    
    group = parser.add_argument_group(title="Point Appearance")
    group.add_argument(
        "--scale_points", default=1.0, type=float,
        help="Scale the size of the points.  Default 1.0")
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

    assert args.legend_inset >= 0 and args.legend_inset < 10
    if args.legend_loc is None:
        args.legend_loc = "bottomright"
    assert args.scale_points > 0 and args.scale_points < 20

    if args.default_color:
        assert len(args.default_color) == 7
        assert args.default_color[0] == "#"

    # Read the data file.
    # List of (name, values).
    MATRIX = read_prism_file(args.datafile)
    
    height = args.height or 2400
    width = args.width or 3200

    # Pull out the values and colors for the plot.
    default_color = "#000000"
    if args.default_color:
        default_color = args.default_color
        

    # Start R and set up the environment.
    R = jmath.start_R()
    R("library(beeswarm)")

    main = jmath.R_var("NA")
    if args.title:
        main = args.title
    sub = ""
    xlab = ""
    if args.xlab:
        xlab = args.xlab
    ylab = ""
    if args.xlab:
        ylab = args.ylab
    
    lwd_box = 2
    lwd_axis = 2
    #lwd_regr = 3
    cex = 1.0 * args.scale_points
    cex_lab = 1.5
    cex_main = 2.0
    cex_sub = 1.0

    
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

    R("X <- list()")
    for title, values in MATRIX:
        title_h = hashlib.hash_var(title)
        jmath.R_equals(values, "x")
        R('X[["%s"]] <- x' % title_h)

    keywds = {
        "cex.axis" : cex_lab,   # Y-axis
        "cex.names" : cex_lab,  # X-axis
        }
    jmath.R_fn(
        "beeswarm", jmath.R_var("X"),
        main="", xlab="", ylab="", pch=19, cex=cex,
        #axes=jmath.R_var("FALSE"),
        RETVAL="x", **keywds)
    # Make plot area solid white.
    #jmath.R('usr <- par("usr")')
    #jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    #jmath.R_fn(
    #    "hist", jmath.R_var("X"), plot=jmath.R_var("FALSE"),
    #    main=main, xlab="", ylab="", axes=jmath.R_var("FALSE"),
    #    add=jmath.R_var("TRUE"))

    # Calculate correlation, and other statistics.
    # TODO: Should calculate this for each series.
    #r = jmath.R("cor(X, Y)")
    #p_value = jmath.R("cor.test(X, Y)$p.value")
    #r = r[0]
    #p_value = p_value[0]
    #print "R = %.2f" % r
    #print "p = %.2g" % p_value

    if not args.no_box:
        jmath.R_fn("box", lwd=lwd_box)
    jmath.R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
        **{ "cex.lab" : cex_lab, "cex.main" : cex_main, "cex.sub" : cex_sub })
    R("par(op)")
    jmath.R_fn("dev.off")

    # Print out some statistics.
    #means = []
    #for (title, values) in MATRIX:
    #    m = jmath.mean(values)
    #    means.append(m)
    #header = "Classes", "Means", "p-value"
    

        
if __name__ == '__main__':
    main()
