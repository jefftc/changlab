#!/usr/bin/env python

def _calc_p(combo2common, name1, name2, total_items):
    from genomicode import jmath
    
    #             not set 2  in set 2
    # not set 1      x11       x12       x1d
    # in set 1       x21       x22       x2d
    #                xd1       xd2       xdd
    xdd = total_items
    x2d = len(combo2common[(name1, name1)])
    xd2 = len(combo2common[(name2, name2)])
    x22 = len(combo2common[(name1, name2)])
    x21 = x2d - x22
    x12 = xd2 - x22
    x11 = xdd - x12 - x21 - x22

    p = jmath.fisher_test(x11, x12, x21, x22)
    return p


def _name_replace(name_replace, all_names, name2genes):
    if not name_replace:
        return all_names, name2genes
    replace = []  # list of (from, no)
    for x in name_replace:
        x = x.split(",", 1)
        assert len(x) == 2
        from_, to = x
        replace.append((from_, to))
    new_all_names = []
    new_name2genes = {}
    for name in all_names:
        new_name = name
        for (from_, to) in replace:
            new_name = new_name.replace(from_, to)
        assert new_name not in new_all_names, "dup"
        new_all_names.append(new_name)
        new_name2genes[new_name] = name2genes[name]
    return new_all_names, new_name2genes


def product_genesets(all_names, num_to_compare):
    import itertools
    
    x = [all_names] * num_to_compare

    for combo in itertools.product(*x):
        assert len(combo) == num_to_compare
        # Need to preserve order.
        #uniq_combo = {}.fromkeys(combo).keys()
        uniq_combo = []
        for x in combo:
            if x not in uniq_combo:
                uniq_combo.append(x)
        yield tuple(uniq_combo)


def match_gene_sets(genesets, delimiter):
    # Return a tuple.  First element is list of genesets where the _UP
    # and _DN (_DOWN) are combined.  Second element are the pretty
    # names.  Preserve the order of the gene sets.

    # Make sure delimiter is not already used.
    for x in genesets:
        assert delimiter not in x, "Error: bad delimiter"
    
    pretty_names = genesets[:]
    i = 0
    while i < len(genesets)-1:
        gs1 = genesets[i]
        ugs1 = gs1.upper()
        if ugs1.endswith("_DOWN"):
            ugs1 = ugs1[:-5] + "_DN"
        
        # If there is already a positive and negative, skip it.
        if gs1.find(delimiter) >= 0:
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
        x = "%s%s%s" % (gs1, delimiter, gs2)
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
    filename, all_names, name2genes, 
    args_margin, args_label_size, args_count_size):
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

    if len(all_names) in [2, 3, 5]:
        varnames = ["A", "B", "C", "D", "E"]
        for i in range(len(all_names)):
            n = all_names[i]
            R_equals(name2genes[n], varnames[i])
        #n1, n2, n3 = all_names
        #R_equals(name2genes[n1], "A")
        #R_equals(name2genes[n2], "B")
        #R_equals(name2genes[n3], "C")
        if len(all_names) == 2:
            R('x <- list(A=A, B=B)')
        elif len(all_names) == 3:
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

        if len(all_names) == 2:
            fill = ["cornflowerblue", "darkorchid1"]
            cat_col = ["cornflowerblue", "darkorchid1"]
            margin = 0.10*args_margin
            cat_cex = 0.75*args_label_size
            cex = 0.65*args_count_size
        elif len(all_names) == 3:
            fill = ["cornflowerblue", "green", "yellow"]
            cat_col = ["darkblue", "darkgreen", "orange"]
            margin = 0.10*args_margin
            cat_cex = 0.75*args_label_size
            cex = 0.65*args_count_size
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
        "--num_to_compare", type=int, default=2, 
        help="Num gene sets to compare at once.  "
        "Default is 2 (pairwise comparisons).")
    parser.add_argument(
        "--geneset", default=[], action="append",
        help="Which gene sets to include in the VENN diagram.  "
        "If automatch, this is the name without the _UP or _DN suffix.")
    parser.add_argument(
        "--all_genesets", action="store_true",
        help="Calculate the intersection of all gene sets.")
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
    group.add_argument(
        "--name_replace", default=[], action="append",
        help="For the plot, replace a string in the geneset name with "
        "another string.  Format: <from>,<to>  (MULTI)")

    group = parser.add_argument_group(title="Statistics")
    #group.add_argument(
    #    "--calc_p", action="store_true",
    #    help="Calculate a p-value for significance of overlap.")
    group.add_argument(
        "--p_num_items", default=0, type=int,
        help="Total number of items in the set.")

    args = parser.parse_args()
    assert os.path.exists(args.geneset_file), \
           "File not found: %s" % args.geneset_file
    assert args.num_to_compare >= 2 and args.num_to_compare <= 5

    assert args.margin > 0 and args.margin < 100
    assert args.label_size > 0 and args.label_size < 10
    assert args.count_size > 0 and args.count_size < 10

    assert not (args.geneset and args.all_genesets)
    if not args.all_genesets:
        assert len(args.geneset) > 1, "Must compare multiple gene sets."
        assert len(args.geneset) >= args.num_to_compare

    assert args.p_num_items >= 0

    name2genes = {}
    all_names = []  # preserve order of gene sets
    for x in genesetlib.read_genesets(args.geneset_file):
        name, desc, genes = x
        name2genes[name] = genes
        all_names.append(name)

    if args.automatch:
        DELIMITERS = [",", ";", "-"]
        delimiter = None
        for delim in DELIMITERS:
            for x in all_names:
                if delim in x:
                    break
            else:
                delimiter = delim
                break
        assert delimiter, "I could not find a delimiter"
        all_names, pretty_names = match_gene_sets(all_names, delimiter)
        # Merge the gene sets.
        merged = {}
        for name, pretty_name in zip(all_names, pretty_names):
            if delimiter not in name:  # not automatched
                genes = name2genes[name]
            else:
                x = name.split(delimiter)
                assert len(x) == 2, "Bad split: %s %s" % (name, delimiter)
                n1, n2 = x
                g1, g2 = name2genes[n1], name2genes[n2]
                x = g1 + g2
                genes = sorted({}.fromkeys(x))
            merged[pretty_name] = genes
        name2genes = merged
        all_names = pretty_names

    # Select only the gene sets of interest.
    if args.geneset:
        for x in args.geneset:
            assert x in name2genes, "Missing geneset: %s" % x
        all_names = args.geneset

    # Count all pairwise intersections.
    combo2common = {}  # (gs1, gs2[, ...]) -> list of common genes
    for combo in product_genesets(all_names, args.num_to_compare):
        assert len(combo) >= 1
        common = set(name2genes[combo[0]])
        for i in range(1, len(combo)):
            x = name2genes[combo[i]]
            common = common.intersection(x)
        combo2common[combo] = sorted(common)
        # Hack to make writing out matrix below easier.
        if args.num_to_compare == 2:
            if len(combo) == 1:
                x = (combo[0], combo[0])
                combo2common[x] = sorted(common)

    # Write out the matrix.
    if args.num_to_compare == 2:
        header = ["Count"] + all_names
        print "\t".join(header)
        for name1 in all_names:
            x = [len(combo2common[(name1, name2)]) for name2 in all_names]
            x = [name1] + x
            assert len(x) == len(header)
            print "\t".join(map(str, x))
    else:
        print "No matrix for comparing %d gene sets." % args.num_to_compare
    print

    # Write out the statistics.
    if args.p_num_items > 0:
        if args.num_to_compare != 2:
            raise NotImplementedError, "Can't calc p-values with > 2 groups."
        for name1 in all_names:
            for name2 in all_names:
                if name1 == name2:
                    continue
                if name1 > name2:
                    continue
                p = _calc_p(combo2common, name1, name2, args.p_num_items)
                print "%s %s p=%s" % (name1, name2, p)
        

    # Write to an output file.
    if args.outfile:
        genesets = []  # list of GeneSet objects
        seen = {}
        for combo in product_genesets(all_names, args.num_to_compare):
            combo_s = tuple(sorted(combo))
            if combo_s in seen:
                continue
            seen[combo_s] = 1
            
            genes = combo2common[combo]
            name = "___".join(combo)
            x = genesetlib.GeneSet(name, "na", genes)
            genesets.append(x)

        write_fn = genesetlib.write_gmx
        if args.outfile.lower().endswith(".gmt"):
            write_fn = genesetlib.write_gmt
        write_fn(args.outfile, genesets)

    if args.plotfile:
        x = _name_replace(args.name_replace, all_names, name2genes)
        all_names, name2genes = x
        draw_venn(
            args.plotfile, all_names, name2genes, 
            args.margin, args.label_size, args.count_size)


if __name__ == '__main__':
    main()
