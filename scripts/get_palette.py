#!/usr/bin/env python


def main():
    import argparse
    from genomicode import parallel
    import arrayplot2

    parser = argparse.ArgumentParser()

    DEFAULT_COLOR_SCHEME = "brewer-rdylbu-div"
    COLOR_SCHEMES = [x[0] for x in arrayplot2.SCHEME2FN]
    assert DEFAULT_COLOR_SCHEME in COLOR_SCHEMES

    #x = arrayplot2._fmt_color_schemes(COLOR_SCHEMES, DEFAULT_COLOR_SCHEME)
    parser.add_argument(
        "--color", dest="color_scheme", 
        default=DEFAULT_COLOR_SCHEME, choices=COLOR_SCHEMES,
        help="Choose the color scheme to use.  "
        "brewer_qual_set1 works best when n is 8.")
    parser.add_argument(
        "-n", "--num_colors", default=8, type=int,
        help="Number of colors to get.")

    # Parse the input arguments.
    args = parser.parse_args()
    assert args.num_colors >= 2 and args.num_colors < 1000

    print "--color=%s" % args.color_scheme

    color_fn = arrayplot2.get_color_scheme_fn(args.color_scheme)
    for i in range(args.num_colors):
        perc = 1.0 - 1.00 / (args.num_colors-1) * i
        # Return as r, g, b (0-255)
        x = arrayplot2._get_color(perc, color_fn, num_colors=args.num_colors)
        r, g, b = x
        print "%d, %d, %d" % (r, g, b)

    # brewer-reds-seq goes from dark to light.

    # Make a heatmap with the colors.
    MATRIX_FILE = "palette.txt"
    HEATMAP_FILE = "heatmap.png"
    handle = open(MATRIX_FILE, 'w')
    header = "Color", "Value"
    print >>handle, "\t".join(header)
    for i in range(args.num_colors):
        perc = 1.0 - 1.00 / (args.num_colors-1) * i
        x = "COLOR_%d" % i, perc
        print >>handle, "\t".join(map(str, x))
    handle.close()

    cmd = [
        "arrayplot2.py",
        "--no_autoscale",
        "-m", "2.0",
        "-s", "-0.5",
        "--color", args.color_scheme,
        MATRIX_FILE,
        HEATMAP_FILE,
        ]
    cmd = " ".join(cmd)
    x = parallel.sshell(cmd)
    if x:
        print x,



if __name__ == '__main__':
    main()
