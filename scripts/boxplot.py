#!/usr/bin/env python

import os

def main():
    from optparse import OptionParser, OptionGroup
    
    import arrayio
    from genomicode import config
    from genomicode import jmath

    usage = "usage: %prog [options] filename outfile.pdf"
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "", "--title", dest="title", default=None, type="str",
        help="Put a title on the plot.")
    parser.add_option(
        "", "--width", dest="width", default=None, type="int",
        help="Width (in pixels) of the plot.")
    parser.add_option(
        "-v", "--verbose", dest="verbose", default=False, action="store_true",
        help="")
    
    # Parse the input arguments.
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Please specify an infile and an outfile.")
    elif len(args) > 2:
        parser.error("Too many input parameters (%d)." % len(args))
    filename, outfile = args
    if not os.path.exists(filename):
        parser.error("I could not find file %s." % filename)
    if options.width is not None:
        assert options.width > 10, "too small"
        assert options.width < 4096*16, "width too big"

    MATRIX = arrayio.read(filename)
    assert MATRIX.nrow() and MATRIX.ncol(), "Empty matrix."

    # Start R and set up the environment.
    R = jmath.start_R()
    path = config.changlab_Rlib
    plotlib = os.path.join(path, "plotlib.R")
    assert os.path.exists(plotlib), "I cannot find: %s" % plotlib

    jmath.R_equals(MATRIX._X, "X")
    jmath.R_equals(MATRIX.col_names(arrayio.COL_ID), "sample.name")

    jmath.R_fn("source", plotlib)
    jmath.R_fn(
        "bitmap", outfile, type="pdfwrite", height=1600, width=1600,
        units="px", res=300)
    main = jmath.R_var("NA")
    if options.title:
        main = options.title
    jmath.R_fn(
        "my.boxplot", jmath.R_var("X"), ylab="Gene Expression", main=main,
        labels=jmath.R_var("sample.name"))
    jmath.R_fn("dev.off")


if __name__ == '__main__':
    main()
