#!/usr/bin/env python


def _parse_breaks_seq(s):
    # Format: <start>,<stop>,<skip>
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

    # density from R doesn't sum up to 1.  (e.g. sum to 2).
    # Recalculate so that it sums to 1.
    total = sum(counts)
    for i in range(len(density)):
        density[i] = counts[i] / float(total)
    header = ["Mids", "Left", "Counts", "Density"]
    x = [mids, breaks, counts, density]
    x = jmath.transpose(x)
    x = [header] + x
    handle = open(filename, 'w')
    for x in x:
        print >>handle, "\t".join(map(str, x))


PALETTE2FN = [
    ("red", "red_shade"),
    ("white", "white_shade"),
    ("red-green", "rg_array_colors"),
    ("blue-yellow", "by_array_colors"),
    ("red-green-soft", "red_green_soft"),
    ("red-blue-soft", "red_blue_soft"),
    ("matlab", "matlab_colors"),
    ("bild", "bild_colors"),
    ("genepattern", "broad_colors"),
    ("genespring", "genespring_colors"),
    ("yahoo", "yahoo_weather_colors"),
    ("brewer-prgn-div", "brewer_prgn_div"),
    ("brewer-rdbu-div", "brewer_rdbu_div"),
    ("brewer-rdylbu-div", "brewer_rdylbu_div"),
    ("brewer-rdylgn-div", "brewer_rdylgn_div"),
    ("brewer-spectral-div", "brewer_spectral_div"),
    ]

def _fmt_palettes():
    palettes = [x[0] for x in PALETTE2FN]
    if len(palettes) == 1:
        return palettes[0]

    # Make it comma separated.
    x = ", ".join(palettes[:-1])
    if len(palettes) > 2:
        x = "%s," % x
    x = "%s or %s" % (x, palettes[-1])
    return x


def get_palette_fn(name):
    # Choose the palette.
    from genomicode import colorlib

    for palette, fn_name in PALETTE2FN:
        if palette == name:
            break
    else:
        raise AssertionError, "Unknown color scheme: %s" % name
    assert hasattr(colorlib, fn_name)
    x = getattr(colorlib, fn_name)
    return x


def _make_col_palette(palette_name, num_bars, symmetric):
    # Format the col variable for the R hist function.
    from genomicode import colorlib

    color_fn = get_palette_fn(palette_name)
    if not symmetric:
        x = color_fn(num_bars)
        print x
    else:
        # Color order goes from left to right.  The colors that was on
        # the left should now be in the middle.
        if num_bars % 2:  # Odd.
            x = color_fn(num_bars/2+1)
            x = list(reversed(x)) + x[:-1]
        else:
            x = color_fn(num_bars/2)
            x = list(reversed(x)) + x
    x = [colorlib.rgb2hex(x) for x in x]
    x = [x.replace("0x", "#") for x in x]
    assert len(x) == num_bars
    return x

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
    parser.add_argument(
        "--prism_file", help="Write Prism-formatted results to this file.")

    group = parser.add_argument_group(title="Calculations")
    group.add_argument(
        "--breaks_seq",
        help="Set the breakpoints.  Format: <start>,<stop>,<skip>.")
    group.add_argument(
        "--num_breaks", type=int, help="Number of breakpoints.")
    group.add_argument(
        "--ymax", type=int,
        help="Set the maximum value for the Y axis.")
    
    group = parser.add_argument_group(title="Plot Labels")
    group.add_argument("--title", help="Put a title on the plot.")
    group.add_argument("--xlab", help="Label the X-axis.")
    group.add_argument(
        "--xlabel_size", default=1.0, type=float,
        help="Scale the size of the labels on X-axis.  Default 1.0.")
    group.add_argument(
        "--xlabel_off", action="store_true", help="Do not label the X axis.")
    group.add_argument(
        "--ylabel_off", action="store_true", help="Do not label the Y axis.")
    group.add_argument(
        "--xtick_label_off", action="store_true",
        help="Do not draw the tick labels on the X axis.")

    group = parser.add_argument_group(title="Colors")
    group.add_argument(
        "--bar_color",  help="Set the color of the bars.  Default #FFFFFF")
    x = _fmt_palettes()
    group.add_argument(
        "--bar_palette", help="Color the bars according to a palette: %s." % x)
    group.add_argument(
        "--symmetric_palette", action="store_true",
        help="Make the color symmetric.")

    group = parser.add_argument_group(title="Appearance")
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
        "--xaxis_off", action="store_true", help="Do not show the X axis.")
    group.add_argument(
        "--yaxis_off", action="store_true", help="Do not show the Y axis.")


    # Parse the input arguments.
    args = parser.parse_args()
    if not os.path.exists(args.datafile):
        parser.error("File not found: %s" % args.datafile)
    assert not (args.breaks_seq and args.num_breaks)
    if args.num_breaks:
        assert args.num_breaks >= 2 and args.num_breaks <= 1000
    if args.width is not None:
        assert args.width > 10, "too small"
        assert args.width < 4096*16, "width too big"
    if args.height is not None:
        assert args.height > 10, "too small"
        assert args.height < 4096*16, "height too big"
    assert args.mar_bottom > 0 and args.mar_bottom < 10
    assert args.mar_left > 0 and args.mar_left < 10
    assert args.xlabel_size > 0 and args.xlabel_size < 10
    assert not (args.bar_color and args.bar_palette)
    assert not args.symmetric_palette or args.bar_palette
    assert args.ymax is None or args.ymax > 0


    height = args.height or 2400
    width = args.width or 3200

    MATRIX = AnnotationMatrix.read(args.datafile, False)
    assert MATRIX.num_headers() and MATRIX.num_annots(), "Empty matrix."
    assert args.header in MATRIX.headers, "header not found: %s" % args.header

    # Pull out the values for the histogram.
    x = MATRIX[args.header]
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
    xtick_labels = jmath.R_var("TRUE")
    ytick_labels = jmath.R_var("TRUE")

    if args.xlabel_off:
        xlab = ""
    if args.ylabel_off:
        ylab = ""
    if args.xtick_label_off:
        xtick_labels = jmath.R_var("FALSE")

    breaks = "Sturges"
    if args.breaks_seq:
        breaks = _parse_breaks_seq(args.breaks_seq)
        value_min, value_max = min(breaks), max(breaks)
        jmath.R_equals(breaks, "breaks")
        breaks = jmath.R_var("breaks")
    if args.num_breaks:
        breaks = args.num_breaks

    if value_min is not None:
        values = [x for x in values if x >= value_min]
    if value_max is not None:
        values = [x for x in values if x < value_max]

    lwd = 2
    cex_lab = 1.5
    cex_main = 2.0
    cex_sub = 1.5
    ylim = jmath.R_var("NULL")
    if args.ymax is not None:
        ylim = [0, args.ymax]

    assert values
    jmath.R_equals(values, "X")

    # Figure out the colors.  Do it after X is assigned.
    col = jmath.R_var("NULL")
    if args.bar_color:
        assert args.bar_color.startswith("#")
        col = args.bar_color
    elif args.bar_palette:
        # Figure out how many breaks there are.  Number of bars is num
        # breaks + 1.
        jmath.R_fn(
            "hist", jmath.R_var("X"), breaks=breaks, plot=jmath.R_var("FALSE"),
            RETVAL="x")
        breaks = [x for x in R["x"].rx2("breaks")]
        num_bars = len(breaks) + 1
        col = _make_col_palette(
            args.bar_palette, num_bars, args.symmetric_palette)

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
        "hist", jmath.R_var("X"), breaks=breaks, main="", xlab="", ylab="",
        ylim=ylim, axes=jmath.R_var("FALSE"), col=col, RETVAL="x")
    # Make plot area solid white.
    #jmath.R('usr <- par("usr")')
    #jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    #jmath.R_fn(
    #    "hist", jmath.R_var("X"), plot=jmath.R_var("FALSE"),
    #    main=main, xlab="", ylab="", axes=jmath.R_var("FALSE"),
    #    add=jmath.R_var("TRUE"))
    
    #jmath.R_fn("box", lwd=lwd)
    # x-axis
    if not args.xaxis_off:
        jmath.R_fn(
            "axis", 1, lwd=lwd, labels=xtick_labels, **{ "cex.axis" : 1.5 })
    # y-axis
    if not args.yaxis_off:
        jmath.R_fn(
            "axis", 2, lwd=lwd, labels=ytick_labels, **{ "cex.axis" : 1.5 })
    jmath.R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
        **{ "cex.lab" : cex_lab, "cex.main" : cex_main, "cex.sub" : cex_sub })
    R("par(op)")
    jmath.R_fn("dev.off")

    if args.prism_file:
        write_prism_file(args.prism_file, R["x"])
        
if __name__ == '__main__':
    main()
