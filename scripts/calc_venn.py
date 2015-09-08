#!/usr/bin/env python


def match_gene_sets(genesets):
    # Return a tuple.  First element is list of genesets where the _UP
    # and _DN (_DOWN) are combined.  Second element are the pretty
    # names.  Preserve the order of the gene sets.
    
    pretty_names = genesets[:]
    i = 0
    while i < len(genesets)-1:
        gs1 = genesets[i]
        ugs1 = gs1.upper()
        if ugs1.endswith("_DOWN"):
            ugs1 = ugs1[:-5] + "_DN"
        
        # If there is already a positive and negative, skip it.
        if gs1.find(",") >= 0:
            i += 1
            continue

        match = False
        for j in range(i+1, len(genesets)):
            gs2 = genesets[j]
            ugs2 = gs2.upper()
            if ugs2.endswith("_DOWN"):
                ugs2 = ugs2[:-5] + "_DN"

            # Sort them for easier comparison.
            s1, s2 = ugs1, ugs2
            if s1 > s2:
                s1, s2 = s2, s1
            if s1.endswith("_DN") and s2.endswith("_UP") and \
                   s1[:-3] == s2[:-3]:
                match = True
                break
        if not match:
            i += 1
            continue
        x = "%s,%s" % (gs1, gs2)
        genesets[i] = x
        del genesets[j]
        x = gs1
        x = x.replace("_UP", "")
        x = x.replace("_DN", "")
        x = x.replace("_DOWN", "")
        pretty_names[i] = x
        del pretty_names[j]
    return genesets, pretty_names


def draw_venn(
    filename, all_names, name2genes, pair2common, args_margin,
    args_label_size, args_count_size):
    import sys
    import StringIO
    from genomicode import jmath
    R_fn = jmath.R_fn
    R_var = jmath.R_var
    R_equals = jmath.R_equals
    
    R = jmath.start_R()

    # Prevent R from writing junk to the screen.
    handle = StringIO.StringIO()
    old_stdout = sys.stdout
    sys.stdout = handle
    R_fn('library', R_var('VennDiagram'))
    sys.stdout = old_stdout

    if len(all_names) == 2:
        raise NotImplementedError
    elif len(all_names) == 3 or len(all_names) == 5:
        varnames = ["A", "B", "C", "D", "E"]
        for i in range(len(all_names)):
            n = all_names[i]
            R_equals(name2genes[n], varnames[i])
        #n1, n2, n3 = all_names
        #R_equals(name2genes[n1], "A")
        #R_equals(name2genes[n2], "B")
        #R_equals(name2genes[n3], "C")
        if len(all_names) == 3:
            R('x <- list(A=A, B=B, C=C)')
        elif len(all_names) == 5:
            R('x <- list(A=A, B=B, C=C, D=D, E=E)')
        else:
            raise NotImplementedError
        for i in range(len(all_names)):
            n = all_names[i]
            R('names(x)[%d] <- "%s"' % (i+1, n))
        #R('names(x)[1] <- "%s"' % n1)
        #R('names(x)[2] <- "%s"' % n2)
        #R('names(x)[3] <- "%s"' % n3)

        cex = 1*args_count_size         # Size of number in each circle.
        cat_cex = 1.5*args_label_size   # Size of category labels.
        margin = 0.05*args_margin   # Amount of space around plot.
        # Bigger margin is smaller figure.
        
        if len(all_names) == 3:
            fill = ["cornflowerblue", "green", "yellow"]
            cat_col = ["darkblue", "darkgreen", "orange"]
        elif len(all_names) == 5:
            fill = [
                "dodgerblue", "goldenrod1", "darkorange1", "seagreen3",
                "orchid3"]
            cat_col = [
                "dodgerblue", "goldenrod1", "darkorange1", "seagreen3",
                "orchid3"]
            margin = 0.25*args_margin
            cat_cex = 0.75*args_label_size
            cex = 0.65*args_count_size
        else:
            raise NotImplementedError
        
        params = {
            "col" : "transparent",   # color of outer lines
            
            "fill" : fill,
            "alpha" : 0.50,

            # Number of items.
            #"lty" : "blank",
            "cex" : cex,
            
            # Labels
            "cat.cex" : cat_cex,
            "cat.col" : cat_col,
            "cat.default.pos" : "text",
            "margin" : margin,
            }
        R_fn(
            "venn.diagram", R_var("x"), filename=filename,
            **params)
        
        ## area1 = len(name2genes[n1])
        ## area2 = len(name2genes[n2])
        ## area3 = len(name2genes[n3])
        ## n12 = len(pair2common[(n1, n2)])
        ## n23 = len(pair2common[(n2, n3)])
        ## n13 = len(pair2common[(n1, n3)])
        ## x1 = pair2common[(n1, n2)]
        ## x2 = name2genes[n3]
        ## n123 = len(set(x1).intersection(x2))
        ## category = all_names
        ## fontfamily = ["Helvetica"]*7
        ## cat_fontfamily = ["Helvetica"]*3
        ## #fill = ["blue", "red", "green"]
        ## fill = ["cornflowerblue", "green", "yellow"]
        ## cat_col = ["cornflowerblue", "green", "yellow"]
        
        ## R_fn(
        ##     "bitmap", filename, type="png256",
        ##     height=1600, width=1600, units="px", res=300)
        ## params = {
        ##     "fontfamily" : fontfamily,
        ##     #"col" : "transparent",
        ##     "fill" : fill,
            
        ##     #"lty" : "blank",
        ##     #"cex" : 2,
            
        ##     # Labels
        ##     "category" : category,
        ##     "cat.fontfamily" : cat_fontfamily,
        ##     #"cat.cex" : 2,
        ##     "cat.col" : cat_col,
        ##     }
        ## R_fn(
        ##     "draw.triple.venn", area1, area2, area3, n12, n23, n13, n123,
        ##     **params)
        ## jmath.R_fn("dev.off")
    elif len(all_names) > 5:
        raise AssertionError, "Can't draw venn diagram with %d circles." % \
              len(all_names)
    else:
        raise NotImplementedError, len(all_names)

def main():
    import os
    import argparse
    import itertools
    from genomicode import genesetlib

    # Option for case insensitive comparison?
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "geneset_file", help="File in GMX or GMT format.")
    parser.add_argument(
        "-o", dest="outfile", help="Save the intersection to this "
        "(GMX or GMT) file.")
    #parser.add_argument(
    #    "--upper_diagonal", action="store_true",
    #    help="Calculate upper diagonal only.")
    parser.add_argument(
        "--automatch", action="store_true", 
        help="Will match _UP with _DN (or _DOWN).")

    group = parser.add_argument_group(title="Plot")
    group.add_argument(
        "--plotfile", help="Save a TIFF plot to this file.")
    group.add_argument(
        "--margin", default=1.0, type=float,
        help="Increase or decrease the margin.")
    group.add_argument(
        "--label_size", default=1.0, type=float,
        help="Increase or decrease the font for the labels.")
    group.add_argument(
        "--count_size", default=1.0, type=float,
        help="Increase or decrease the font for the counts.")

    args = parser.parse_args()
    assert os.path.exists(args.geneset_file), \
           "File not found: %s" % args.geneset_file

    assert args.margin > 0 and args.margin < 100
    assert args.label_size > 0 and args.label_size < 10
    assert args.count_size > 0 and args.count_size < 10

    name2genes = {}
    all_names = []  # preserve order of gene sets
    for x in genesetlib.read_genesets(args.geneset_file):
        name, desc, genes = x
        name2genes[name] = genes
        all_names.append(name)

    if args.automatch:
        all_names, pretty_names = match_gene_sets(all_names)
        # Merge the gene sets.
        merged = {}
        for name, pretty_name in zip(all_names, pretty_names):
            x = name.split(",")
            assert len(x) == 2
            n1, n2 = x
            g1, g2 = name2genes[n1], name2genes[n2]
            x = g1 + g2
            genes = sorted({}.fromkeys(x))
            merged[pretty_name] = genes
        name2genes = merged
        all_names = pretty_names

    # Count all pairwise intersections.
    pair2common = {}  # (name1, name2) -> common genes
    for x in itertools.product(all_names, all_names):
        name1, name2 = x
        genes1 = name2genes[name1]
        genes2 = name2genes[name2]
        common = set(genes1).intersection(genes2)
        pair2common[(name1, name2)] = sorted(common)

    # Write out the matrix.
    header = ["Count"] + all_names
    print "\t".join(header)
    for name1 in all_names:
        x = [len(pair2common[(name1, name2)]) for name2 in all_names]
        x = [name1] + x
        assert len(x) == len(header)
        print "\t".join(map(str, x))

    # Write to an output file.
    if args.outfile:
        genesets = []  # list of GeneSet objects
        for i in range(len(all_names) - 1):
            for j in range(i, len(all_names)):
                name1 = all_names[i]
                name2 = all_names[j]
                genes = pair2common[(name1, name2)]
                name = "%s___%s" % (name1, name2)
                x = genesetlib.GeneSet(name, "na", genes)
                genesets.append(x)

        write_fn = genesetlib.write_gmx
        if args.outfile.lower().endswith(".gmt"):
            write_fn = genesetlib.write_gmt
        write_fn(args.outfile, genesets)

    if args.plotfile:
        draw_venn(
            args.plotfile, all_names, name2genes, pair2common,
            args.margin, args.label_size, args.count_size)


if __name__ == '__main__':
    main()
