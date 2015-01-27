#!/usr/bin/env python

import os

def find_gene_names(MATRIX, gene_names):
    # gene_names is a list of names.  Each element can be a
    # comma-separated list of names.  To recover a single list of
    # names, join each member of the list with a comma, and then split
    # by commas.
    x = [x.strip() for x in gene_names]
    x = ",".join(x)
    names = x.split(",")
    I_row, I_col = MATRIX._index(row=names)
    return I_row

def get_pretty_gene_name(MATRIX, gene_i):
    from genomicode import arrayplatformlib as apl

    # Optimization.  Look for known header names and assume they're
    # right.
    headers = ["Gene Symbol", "Gene.Symbol"]
    for h in headers:
        if h in MATRIX.row_names():
            return MATRIX.row_names(h)[gene_i]

    # Try finding a gene symbol.
    header = apl.find_header(MATRIX, apl.GENE_SYMBOL)
    if header is not None:
        gene_name = MATRIX.row_names(header)[gene_i]
        if gene_name:
            return gene_name

    # Try finding the gene ID.
    header = apl.find_header(MATRIX, apl.GENE_ID)
    if header is not None:
        gene_name = MATRIX.row_names(header)[gene_i]
        if gene_name:
            return "GeneID %s" % gene_name

    # Just use a probe ID.
    header = apl.find_header(MATRIX, apl.PROBE_ID)
    if header is not None:
        gene_name = MATRIX.row_names(header)[gene_i]
        if gene_name:
            return "ProbeID %s" % gene_name

    # Use whatever's in the first column.
    if MATRIX.row_names():
        gene_name = MATRIX.row_names(MATRIX.row_names()[0])[gene_i]
        if gene_name:
            return gene_name

    # Just return the index of the gene
    return "Gene %04d" % gene_i

def main():
    import argparse
    
    import arrayio
    from genomicode import config
    from genomicode import colorlib
    from genomicode import jmath
    from genomicode.jmath import R_fn, R_var, R_equals

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("expression_file", help="Gene expression file.")
    parser.add_argument(
        "plot_file", help="Name of image file, e.g. outfile.png.  "
        "Will generate PNG format by default.  If this file name ends with "
        ".pdf, will generate a PDF file instead.")

    parser.add_argument(
        "--title", default=None, help="Put a title on the plot.")
    parser.add_argument(
        "--gene_names", default=[], action="append",
        help="Comma-separated list of IDs (e.g. probes, gene names) "
        "to include.")
    parser.add_argument(
        "--all_genes", default=False, action="store_true",
        help="Plot all genes in the file.")
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
        help="Scale margin at bottom of plot.  Default 1.0 (no scaling).")
    parser.add_argument(
        "--xlabel_size", default=1.0, type=float,
        help="Scale the size of the labels on X-axis.  Default 1.0.")
    parser.add_argument(
        "--xlabel_off", default=False, action="store_true",
        help="Turn off the X labels.")
    parser.add_argument(
        "-v", "--verbose", default=False, action="store_true",
        help="")

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
    assert args.gene_names or args.all_genes, \
           "Please specify some genes to plot."
    assert args.mar_bottom > 0 and args.mar_bottom < 10
    assert args.mar_left > 0 and args.mar_left < 10
    assert args.xlabel_size > 0 and args.xlabel_size < 10

    height = args.height or 1600
    width = args.width or 1600

    MATRIX = arrayio.read(args.expression_file)
    assert MATRIX.nrow() and MATRIX.ncol(), "Empty matrix."

    I = None
    if args.gene_names:
        I = find_gene_names(MATRIX, args.gene_names)
    elif args.all_genes:
        I = range(MATRIX.nrow())
    assert I, "No genes found."
    assert len(I) < 50, "Too many genes."
    MATRIX = MATRIX.matrix(I, None)

    gene_names = [get_pretty_gene_name(MATRIX, i)
                  for i in range(MATRIX.nrow())]

    # Start R and set up the environment.
    R = jmath.start_R()
    path = config.changlab_Rlib
    plotlib = os.path.join(path, "plotlib.R")
    assert os.path.exists(plotlib), "I cannot find: %s" % plotlib
    R_fn("source", plotlib)

    main = R_var("NA")
    if args.title:
        main = args.title
    sub = ""
    xlab = ""
    ylab = "Gene Expression"
    labels = jmath.R_var("FALSE")
    #labels = MATRIX.col_names(arrayio.COL_ID)
    col = R_var("NULL")
    xlim = [1, MATRIX.ncol()+1]
    y_max = jmath.max(jmath.max(MATRIX._X))
    y_min = jmath.min(jmath.min(MATRIX._X))
    ylim = [y_min-1, y_max+1]

    if not args.xlabel_off:
        labels = MATRIX.col_names(arrayio.COL_ID)

    lwd = 2
    las = 3   # vertical labels
    at = R_var("NULL")
    if labels != jmath.R_var("FALSE"):
        at = range(1, len(labels)+1)
    cex_labels = 1*args.xlabel_size
    cex_legend = 1
    cex_lab = 1.5
    cex_sub = 1.5
    x = colorlib.bild_colors(len(gene_names))
    x = [colorlib.rgb2hex(x) for x in x]
    x = [x.replace("0x", "#") for x in x]
    col = x

    R_equals(MATRIX._X, "X")
    R_equals(labels, "labels")
    R_equals(at, "at")

    bm_type = "png16m"
    if args.plot_file.lower().endswith(".pdf"):
        bm_type = "pdfwrite"
    R_fn(
        "bitmap", args.plot_file, type=bm_type,
        height=height, width=width, units="px", res=300)

    # Set the margins.
    x = 5*1.2*args.mar_bottom, 4*1.2*args.mar_left, 4, 2
    mar = [x+0.1 for x in x]
    R_fn("par", mar=mar, RETVAL="op")
    
    R_fn(
        "plot", R_var("NA"), type="n", axes=R_var("FALSE"), xlab="", ylab="",
        xlim=xlim, ylim=ylim)
    jmath.R('usr <- par("usr")')
    jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    jmath.R_fn("box", lwd=lwd)
    jmath.R_fn(
        "axis", 1, lwd=lwd, labels=R_var("labels"),
        at=R_var("at"), las=las, **{ "cex.axis" : cex_labels })
    jmath.R_fn(
        "axis", 2, lwd=lwd, **{ "cex.axis" : 1.5 })
    jmath.R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
        **{ "cex.lab" : cex_lab, "cex.main" : 2.0, "cex.sub" : cex_sub })
    
    for i in range(MATRIX.nrow()):
        y = MATRIX._X[i]
        x = range(1, len(y)+1)
        R_fn("lines", x, y, lwd=lwd, col=col[i])
        R_fn("points", x, y, pch=19, cex=1, col=col[i])

    R_fn(
        "legend", "bottomleft", legend=gene_names, fill=col, cex=1,
        inset=0.05, **{ "box.lwd" : 1.5 })
        
    R_fn("par", R_var("op"))
    R_fn("dev.off")
       

if __name__ == '__main__':
    main()
