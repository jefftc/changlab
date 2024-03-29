#!/usr/bin/env python


# Functions:
# parse_phenotypes
# parse_groups
# parse_ignore_samples
#
# ignore_samples
# calc_association
# center_scores
#
# write_prism_file
# plot_boxplot
# plot_waterfall

def parse_phenotypes(phenotypes):
    # list of phenotypes.
    # e.g. ["STEM", "EMT"]
    # Return (potentially empty) list of phenotypes.
    return phenotypes


def parse_groups(center_by_groups):
    # Return tuple of batch_header, list of group 1, list of group 2.
    # Format: <BATCH_HEADER>;
    #         <GROUP 1 VALUE>[,<GROUP 1 VALUE>,...];
    #         <GROUP 2 VALUE>[,<GROUP 2 VALUE>,...]
    # If not given, return a tuple of None's.
    if not center_by_groups:
        return None, None, None
    
    x = center_by_groups.split(";")
    assert len(x) == 3
    batch_header, x1, x2 = x
    group1 = x1.split(",")
    group2 = x2.split(",")
    return batch_header, group1, group2


def parse_ignore_samples(ignore_samples):
    # Return a tuple of <annot>, <value>
    # Format: <annot>,<value>
    x = ignore_samples.split(",")
    assert len(x) == 2
    return x


def ignore_samples(M, clinical_annots, ignore):
    x = parse_ignore_samples(ignore)
    annot, value = x

    assert annot in clinical_annots, "Missing annot: %s" % annot
    values = clinical_annots[annot]
    
    I = []  # indexes to keep
    for i in range(len(values)):
        if value != values[i]:
            I.append(i)
    assert len(I) < len(values), "I could not find any %s=%s" % (annot, value)

    M_f = M.matrix(None, I)
    annots_f = {}
    for name, values in clinical_annots.iteritems():
        values_f = [values[i] for i in I]
        annots_f[name] = values_f
        assert len(values_f) == M_f.ncol()
    return M_f, annots_f


def calc_association(phenotypes, scores, ignore_insufficient_groups):
    # Return a dictionary with keys:
    # n                    Number of samples.
    # m                    Number of groups.
    # scores               n-list of <float>
    # delta                None or <float>
    # phenotypes           n-list of <string>
    # groups               n-list of <int>  [0, length(group_names)-1]
    # group_names          m-list of <string>  (unique list of pheno)
    # num_samples          dict of <group (int)> : <int>
    # mean_score           dict of <group (int)> : <float>
    # p_value              <float>
    # relationship         <string>
    #
    # May return None if there is only 1 group, and
    # ignore_insufficient_groups is a true value.
    from genomicode import jmath
    from genomicode import sortlib
    
    # Select only the samples with phenotype and score information.
    I1 = [i for (i, x) in enumerate(phenotypes) if x]
    I2 = [i for (i, x) in enumerate(scores) if x != ""]
    I = sorted(set.intersection(set(I1), set(I2)))
    assert I, "No valid samples."

    phenotypes = [phenotypes[i] for i in I]
    scores = [float(scores[i]) for i in I]

    # Figure out the groupings.
    #group_names = sorted({}.fromkeys(phenotypes))
    group_names = sortlib.sort_natural({}.fromkeys(phenotypes))
    if len(group_names) < 2 and ignore_insufficient_groups:
        return None
    assert len(group_names) >= 2, "Need at least 2 groups (%s)." % \
           str(group_names)
    groups = [None] * len(phenotypes)
    for i in range(len(phenotypes)):
        x = group_names.index(phenotypes[i])
        groups[i] = x

    # Calculate the association.
    group2scores = {}  # group -> list of scores
    for i in range(len(scores)):
        n = groups[i]
        if n not in group2scores:
            group2scores[n] = []
        group2scores[n].append(scores[i])

    y = scores
    x = [[0]*len(group_names) for i in range(len(y))]
    for i in range(len(groups)):
        x[i][groups[i]] = 1
    jmath.start_R()
    jmath.R_equals(x, "x")
    jmath.R_equals(y, "y")
    jmath.R("m <- aov(y~x)")
    p_value = jmath.R('summary(m)[[1]][["Pr(>F)"]][1]')[0]

    # Count other things.
    num_samples = {}
    for n in group2scores:
        num_samples[n] = len(group2scores[n])
    mean_score = {}
    for n in group2scores:
        mean_score[n] = jmath.mean(group2scores[n])
    # If there are exactly 2 groups, then find the difference between
    # the two groups.
    delta = None   # list of deltas
    if len(group_names) == 2:
        delta = mean_score[1] - mean_score[0]

    # Figure out the relationship.
    relationship = ""
    assert len(group_names) >= 2
    high_score = None
    for n, score in mean_score.iteritems():
        if high_score is not None and score <= high_score:
            continue
        high_score = score
        x1 = "Higher"
        if len(group_names) > 2:
            x1 = "Highest"
        relationship = "%s in %s" % (x1, group_names[n])
    
    SCORE = {}
    SCORE["n"] = len(scores)
    SCORE["m"] = len(group_names)
    SCORE["scores"] = scores
    SCORE["phenotypes"] = phenotypes
    SCORE["groups"] = groups
    SCORE["group_names"] = group_names
    SCORE["num_samples"] = num_samples
    SCORE["mean_score"] = mean_score
    SCORE["delta"] = delta
    SCORE["p_value"] = p_value
    SCORE["relationship"] = relationship
    return SCORE


def center_scores(scores, batches, phenotypes, group1, group2):
    from genomicode import jmath
    
    assert len(scores) == len(phenotypes)
    assert len(batches) == len(phenotypes)
    batches_all = sorted({}.fromkeys(batches))
    scores_c = [None] * len(scores)
    for batch in batches_all:
        I = [i for i in range(len(batches)) if batches[i] == batch]
        
        scores1, scores2 = [], []
        for i in I:
            pheno = phenotypes[i]
            if pheno in group1:
                scores1.append(scores[i])
            elif pheno in group2:
                scores2.append(scores[i])
            else:
                raise AssertionError, "%s not in groups" % pheno
        assert scores1, "No samples from group1 in batch %s" % batch
        assert scores2, "No samples from group2 in batch %s" % batch

        mean1 = jmath.mean(scores1)
        mean2 = jmath.mean(scores2)
        n = (mean1 + mean2)/2.0
        for i in I:
            scores_c[i] = scores[i] - n
    assert None not in scores_c
    return scores_c


def write_prism_file(filename, scores, phenotypes, group_names):
    for x in phenotypes:
        assert x in group_names

    pheno2scores = {}
    for pheno, score in zip(phenotypes, scores):
        if pheno not in pheno2scores:
            pheno2scores[pheno] = []
        pheno2scores[pheno].append(score)
        
    matrix = []
    matrix.append(group_names)
    x = [[""]*len(group_names) for i in range(len(scores))]
    matrix.extend(x)
    for j in range(len(group_names)):
        scores = pheno2scores.get(group_names[j], [])
        for i in range(len(scores)):
            matrix[i+1][j] = scores[i]
            
    # Delete all the empty rows in the bottom.
    while matrix:
        x = matrix[-1]
        if x == [""]*len(x):
            del matrix[-1]
        else:
            break
        
    handle = open(filename, 'w')
    for x in matrix:
        print >>handle, "\t".join(map(str, x))


## def plot_boxplot(
##     filename, scores, phenotypes, group_names, p_value, gene_id,
##     mar_bottom, mar_left, mar_top):
##     import os
##     from genomicode import jmath
##     from genomicode.jmath import R_fn, R_var, R_equals
##     from genomicode import config

##     xlabel_size = 1.0
##     height = 1600
##     width = 1600

##     pheno2scores = {}
##     for pheno, score in zip(phenotypes, scores):
##         if pheno not in pheno2scores:
##             pheno2scores[pheno] = []
##         pheno2scores[pheno].append(score)

##     R = jmath.start_R()
##     path = config.changlab_Rlib
##     plotlib = os.path.join(path, "plotlib.R")
##     assert os.path.exists(plotlib), "I cannot find: %s" % plotlib
##     R_fn("source", plotlib)

##     #main = R_var("NA")
##     main = gene_id
##     #sub = ""
##     sub = "%.2g" % p_value
##     xlab = ""
##     ylab = "Gene Expression"
##     labels = group_names
##     col = R_var("NULL")

##     lwd = 2
##     las = 3   # vertical labels
##     at = R_var("NULL")
##     if labels:
##         at = range(1, len(labels)+1)
##     cex_labels = 1.25*xlabel_size
##     #cex_legend = 1
##     #cex_lab = 1.5
##     cex_xlab = 2.0
##     cex_ylab = 2.0
##     cex_sub = 1.5

##     R_equals(labels, "labels")
##     R_equals(at, "at")
##     R("X <- list()")
##     for i, n in enumerate(group_names):
##         s = pheno2scores.get(n, [])
##         R_equals(s, "s")
##         R("X[[%d]] <- s" % (i+1))

##     bm_type = "png16m"
##     if filename.lower().endswith(".pdf"):
##         bm_type = "pdfwrite"
##     R_fn(
##         "bitmap", filename, type=bm_type,
##         height=height, width=width, units="px", res=300)
    
##     # Set the margins.
##     # default is 5.1, 4.1, 4.1, 2.1
##     x = 10*mar_bottom, 5*mar_left, 4*mar_top, 2
##     mar = [x+0.1 for x in x]
##     R_fn("par", mar=mar, RETVAL="op")
        
##     R_fn(
##         "boxplot", R_var("X"), col=col, main="", xlab="", ylab="",
##         axes=R_var("FALSE"), pch=19, cex=1, ylim=R_var("NULL"))
##     # Make plot area solid white.
##     jmath.R('usr <- par("usr")')
##     jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
##     R_fn(
##         "boxplot", R_var("X"), col=col, main="", xlab="", ylab="",
##         axes=R_var("FALSE"), pch=19, cex=1, ylim=R_var("NULL"),
##         add=R_var("TRUE"))
    
##     R_fn("box", lwd=lwd)
##     R_fn(
##         "axis", 1, lwd=lwd, labels=R_var("labels"),
##         at=R_var("at"), las=las, **{ "cex.axis" : cex_labels })
##     R_fn(
##         "axis", 2, lwd=lwd, **{ "cex.axis" : 1.5 })
##     R_fn(
##         "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
##         **{ "cex.lab" : cex_xlab, "cex.main" : 2.0, "cex.sub" : cex_sub,
##             "col.sub" : "#A60400" })
##     R("par(op)")
##     R_fn("dev.off")


def plot_waterfall(
    filename, scores, phenotypes, group_names, sample_names, p_value, gene_id,
    mar_bottom, mar_left, mar_top, xlabel_off):
    import os
    from genomicode import jmath
    from genomicode.jmath import R_fn, R_var, R_equals
    from genomicode import config
    from genomicode import colorlib
    import analyze_clinical_outcome as aco

    # Sort by increasing score.
    O = jmath.order(scores)
    scores = [scores[i] for i in O]
    phenotypes = [phenotypes[i] for i in O]
    sample_names = [sample_names[i] for i in O]

    # Plot the colors.
    assert len(group_names) >= 2
    colors = ['#1533AD', '#FFB300']
    if len(group_names) > 2:
        x = colorlib.bild_colors(len(group_names))
        x = [aco.colortuple2hex(*x) for x in x]
        colors = x

    xlabel_size = 1.0
    height = 1600
    width = 1600

    R = jmath.start_R()
    path = config.changlab_Rlib
    plotlib = os.path.join(path, "plotlib.R")
    assert os.path.exists(plotlib), "I cannot find: %s" % plotlib
    R_fn("source", plotlib)

    #main = R_var("NA")
    main = gene_id
    sub = ""
    #sub = "%.2g" % p_value
    xlab = ""
    ylab = "Gene Expression"
    labels = sample_names
    col = [colors[group_names.index(x)] for x in phenotypes]
    x = range(1, len(scores)+1)
    y = scores

    r = (max(y)-min(y))*0.10
    mn = min(y)-r
    mx = max(y)+r
    ylim = (mn, mx)

    lwd = 2
    las = 3   # vertical labels
    cex_labels = 1.25*xlabel_size
    cex_ytick = 1.5
    #cex_legend = 1
    cex_xlab = 2.0
    cex_ylab = 2.0
    cex_sub = 2.0
    legend_x = "topleft"

    R_equals(labels, "labels")
    R_equals(y, "y")

    bm_type = "png16m"
    if filename.lower().endswith(".pdf"):
        bm_type = "pdfwrite"
    R_fn(
        "bitmap", filename, type=bm_type,
        height=height, width=width, units="px", res=300)
    
    # Set the margins.
    xlabel_bottom = 2.0
    if xlabel_off:
        R_equals(R_var("FALSE"), "labels")
        xlabel_bottom = 0.5
    x = 5*mar_bottom*xlabel_bottom, 5*mar_left, 4*mar_top, 2
    mar = [x+0.1 for x in x]
    R_fn("par", mar=mar, RETVAL="op")
    
    R_fn(
        "barplot", R_var("y"), xlab="", ylab="",
        axes=R_var("FALSE"), ylim=ylim, xpd=R_var("FALSE"),
        RETVAL="mp")
    # Make plot area solid white.
    jmath.R('usr <- par("usr")')
    jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    R_fn("box", lwd=lwd)
    mgp = 3, 1.5, 0
    R_fn("par", mgp=mgp, RETVAL="op2")
    R_fn(
        "axis", 1, lwd=lwd, labels=R_var("labels"),
        at=R_var("mp"), las=las, **{ "cex.axis" : cex_labels })
    R("par(op2)")
    R_fn("axis", 2, lwd=lwd, **{ "cex.axis" : cex_ytick })
    R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab="",
        **{ "cex.lab" : cex_xlab, "cex.main" : 2.0, "cex.sub" : cex_sub,
            "col.sub" : "#A60400" })
    R_fn("title", ylab=ylab, **{ "cex.lab" : cex_ylab } )
    R_fn(
        "barplot", R_var("y"), col=col, xlab="", ylab="",
        axes=R_var("FALSE"), ylim=ylim, add=R_var("TRUE"), xpd=R_var("FALSE"))
    R_fn(
        "legend", legend_x, legend=group_names, fill=colors, inset=0.05,
        bg="#FFFFFF")
    R("par(op)")
    R_fn("dev.off")
        

def main():
    import os
    import sys
    import itertools
    import argparse

    import arrayio
    import analyze_clinical_outcome as aco
    import boxplot
    from genomicode import parallel

    parser = argparse.ArgumentParser(
        description="Associate gene expression patterns with a "
        "categorical phenotype.")
    
    parser.add_argument(
        'expression_file',
        help='Either a gene expression file (GCT,CDT,PCL format) or gene set '
        'scores from score_geneset.py.')
    parser.add_argument(
        'phenotype_file', help="Table of phenotypes (tab-delimited text "
        "file).")
    parser.add_argument(
        "--ignore_samples", help="Ignore the samples where an annotation "
        "(a column in the phenotype file) matches a specific value.  "
        "Format:<header>,<value>")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of processors to use.")
    
    group = parser.add_argument_group(title='Analysis')
    group.add_argument(
        '--phenotype', default=[], action='append',
        help='Header in the phenotype file (MULTI).  Format: <header>')
    group.add_argument(
        '--all_phenotypes', action="store_true",
        help="Analyze all phenotypes in the file.")
    parser.add_argument(
        "--ignore_phenotype", default=[], action="append",
        help="Ignore this column in the phenotype file.  "
        "Helpful to get rid of the sample column when using "
        "--all_phenotypes.  Format: <header>  (MULTI)")
    group.add_argument(
        '--ignore_insufficient_groups', action="store_true",
        help="If a phenotype only has one group, then ignore it rather "
        "than raising an error.")
    group.add_argument(
        '--gene', default=[], action='append',
        help='Comma separated name or ID of genes to analyze.  '
        'I will search for this gene in the annotations of the '
        'expression_file.  '
        'You can use this parameter multiple times to search more genes.')
    group.add_argument(
        "--empty_vs_filled", action="store_true",
        help="Instead of categorizing by the contents of the cells, "
        "compare the ones that are empty against the ones that are filled.")
    group.add_argument(
        "--all_genes", action="store_true",
        help="Run analysis on all genes in this file.")
    group.add_argument(
        '--geneset', default=[], action='append',
        help='Name of the geneset to analyze. To specify multiple gene sets, '
        'use this parameter multiple times.')
    group.add_argument(
        "--center_by_phenotype",
        help="Center the scores or gene expression values seen for a "
        "phenotype to 0.  Only one --phenotype can be analyzed in this way "
        "at a time.  This phenotype should have two possible values.  "
        "If there are more values, they need to be merged into two groups.  "
        "Each phenotype must be seen in each BATCH.  "
        "Format: <BATCH_HEADER>;<PHENO 1 VALUE>[,<PHENO 1 VALUE>,...];"
        "<PHENO 2 VALUE>[,<PHENO 2 VALUE>,...]")
    group = parser.add_argument_group(title='Output')
    group.add_argument(
        '-o', dest='filestem', default=None,
        help='Prefix used to name files.  e.g. "myanalysis".')
    group.add_argument(
        "--gene_header", action="append", default=[],
        help="When naming the output file, use the gene name(s) under this "
        "Header (MULTI).  If not given, will try to use a combination of the "
        "probe ID and gene symbol.")

    group = parser.add_argument_group(title='Formatting the boxplot')
    group.add_argument(
        "--box_mar_left", default=1.0, type=float,
        help="Scale margin at left of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--box_mar_bottom", default=1.0, type=float,
        help="Scale margin at bottom of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--box_mar_top", default=1.0, type=float,
        help="Scale margin at top of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--water_mar_left", default=1.0, type=float,
        help="Scale margin at left of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--water_mar_bottom", default=1.0, type=float,
        help="Scale margin at bottom of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--water_mar_top", default=1.0, type=float,
        help="Scale margin at top of plot.  Default 1.0 (no scaling).")
    group.add_argument(
        "--water_xlabel_off", action="store_true",
        help="Do not label the X axis on the waterfall plot.")

    ## group.add_argument(
    ##     '--km_title', default=None, help='Title for the Kaplan-Meier plot.')
    ## group.add_argument(
    ##     '--km_title_size', default=1.0, type=float,
    ##     help='Scale the size of the title.  Default 1.0 (no scaling).')
    ## group.add_argument(
    ##     '--km_mar_title', default=1.0, type=float, 
    ##     help="Scale margin for the title.  Default 1.0 (no scaling).")
    ## group.add_argument(
    ##     '--km_subtitle_size', default=1.0, type=float,
    ##     help='Scale the size of the subtitle.  Default 1.0 (no scaling).')
    ## group.add_argument(
    ##     '--km_mar_subtitle', default=1.0, type=float, 
    ##     help="Scale margin for the subtitle.  Default 1.0 (no scaling).")
    ## group.add_argument(
    ##     '--km_xlab', default=None, 
    ##     help='x-axis label for the Kaplan-Meier plot.')
    ## group.add_argument(
    ##     '--km_ylab', default=None, 
    ##     help='y-axis label for the Kaplan-Meier plot.')
    ## group.add_argument(
    ##     '--km_legend_size', default=1.0, type=float,
    ##     help='Scale the size of the legend.  Default 1.0 (no scaling).')
    
    args = parser.parse_args()

    # Check inputs.
    assert args.expression_file, (
        'Please specify a gene expression or gene set score file.')
    assert os.path.exists(args.expression_file), "File not found: %s" % \
           args.expression_file
    assert args.phenotype_file, 'Please specify a phenotype file.'
    assert os.path.exists(args.phenotype_file), "File not found: %s" % \
           args.phenotype_file
    assert args.num_procs >= 1 and args.num_procs < 100

    assert args.phenotype or args.all_phenotypes, \
           'Please specify the phenotype to analyze.'
    assert not (args.phenotype and args.all_phenotypes)
    assert args.gene or args.geneset or args.all_genes, \
           'Please specify a gene or gene set.'
    assert not (args.gene and args.all_genes)
    has_gene = args.gene or args.all_genes
    assert not (has_gene and args.geneset), \
        'Please specify either a gene or a gene set, not both.'

    assert args.box_mar_bottom > 0 and args.box_mar_bottom < 10
    assert args.box_mar_left > 0 and args.box_mar_left < 10
    assert args.box_mar_top > 0 and args.box_mar_top < 10
    assert args.water_mar_bottom > 0 and args.water_mar_bottom < 10
    assert args.water_mar_left > 0 and args.water_mar_left < 10
    assert args.water_mar_top > 0 and args.water_mar_top < 10
    ## assert args.km_title_size > 0 and args.km_title_size < 10
    ## assert args.km_mar_title > 0 and args.km_mar_title < 10
    ## assert args.km_subtitle_size > 0 and args.km_subtitle_size < 10
    ## assert args.km_mar_subtitle > 0 and args.km_mar_subtitle < 10
    ## assert args.km_legend_size > 0 and args.km_legend_size < 10
    
    # Clean up the input.
    phenotypes = parse_phenotypes(args.phenotype)
    genes = aco.parse_genes(args.gene)
    gene_sets = aco.parse_gene_sets(args.geneset)
    x = parse_groups(args.center_by_phenotype)
    center_batch, center_group1, center_group2 = x
    filestem = aco.parse_filestem(args.filestem)

    if center_batch:
        assert len(phenotypes) == 1, \
               "Only 1 phenotype can be centered by groups."

    # Read the input files.
    M = aco.read_expression_or_geneset_scores(
        genes, args.all_genes, gene_sets, args.expression_file)
    x = aco.read_clinical_annotations(M, args.phenotype_file)
    M, clinical_annots = x

    # Filter the phenotype files.
    if args.ignore_samples:
        x = ignore_samples(M, clinical_annots, args.ignore_samples)
        M, clinical_annots = x

    if args.all_phenotypes:
        phenotypes = sorted(clinical_annots)
    phenotypes = [x for x in phenotypes if x not in args.ignore_phenotype]

    # Make sure at least one of the phenotypes are in the clinical
    # annotations.
    x = [x for x in phenotypes if x in clinical_annots]
    assert x, "Could not find phenotypes: %s" % ", ".join(phenotypes)
    phenotypes = x

    # Select the genes or gene sets of interest.
    if not args.all_genes:
        x = genes or gene_sets
        M = M.matrix(row=x)
    assert M.nrow(), "I could not find any of the genes or gene sets."

    # Make sure the batch information is valid.
    if center_batch:
        assert center_batch in clinical_annots, "Missing annotation: %s" % \
               center_batch
        assert len(phenotypes) == 1
        pheno = phenotypes[0]
        values = clinical_annots[pheno]
        for x in values:
            assert x in center_group1 or x in center_group2, \
                   "Unknown phenotype: %s" % x

    # Calculate the association of each gene and each phenotype.
    #expression_or_score = "Expression"
    #if gene_sets:
    #    expression_or_score = "Score"

    jobs = []  # list of (function, args, keywds)
    keys = []
    for x in itertools.product(phenotypes, range(M.nrow())):
        pheno_header, i = x
        phenotype = clinical_annots[pheno_header]
        if args.empty_vs_filled:
            x = ["0"] * len(phenotype)
            for j in range(len(phenotype)):
                if phenotype[j].strip():
                    x[j] = "1"
            phenotype = x
        
        scores = M.value(i, None)
        if center_batch:
            batch = clinical_annots[center_batch]
            scores = center_scores(
                scores, batch, phenotype, center_group1, center_group2)

        x = phenotype, scores, args.ignore_insufficient_groups
        x = calc_association, x, {}
        jobs.append(x)
        keys.append((pheno_header, i))
    retvals = parallel.pyfun(jobs, num_procs=args.num_procs)
    assert len(retvals) == len(keys)

    # (header, gene_index) -> returned from calc_association
    gene_phenotype_scores = {}
    for (pheno_header, i), x in zip(keys, retvals):
        if x is None:
            continue
        gene_phenotype_scores[(pheno_header, i)] = x
        

    # Files generated:
    # <filestem>.stats.txt      Or to STDOUT if no <filestem> given.
    # <filestem>.<outcome>.<gene_id>.waterfall.png
    # <filestem>.<outcome>.<gene_id>.boxplot.png
    # <filestem>.<outcome>.<gene_id>.prism.txt        Prism format.

    # Write the output in a table with headers:
    # <headers>            # From the expression or gene set file.
    # Phenotype
    # Groups               # one for each group
    # Num Samples          # one for each group, separated by semicolon
    # Average Expression   # one for each group, separated by semicolon
    # Relationship
    # p-value

    outhandle = sys.stdout
    if filestem:
        outhandle = open("%s.stats.txt" % filestem, 'w')

    # Figure out the header for the table.
    header = M.row_names() + [
        "Phenotype", "Groups", "Num Samples", "Average Expression",
        "Delta", "Relationship", "p-value"]
    print >>outhandle, "\t".join(header)

    # Write out each row of the table.
    for x in itertools.product(phenotypes, range(M.nrow())):
        pheno_header, gene_i = x
        SCORE = gene_phenotype_scores.get((pheno_header, gene_i))
        if not SCORE:   # couldn't calculate.
            continue

        gene_names = [M.row_names(x)[gene_i] for x in M.row_names()]
        phenotype = pheno_header
        group_names = SCORE["group_names"]
        I = range(len(group_names))
        num_samples = [SCORE["num_samples"][x] for x in I]
        mean_score = [SCORE["mean_score"][x] for x in I]
        delta = ""
        if len(group_names) == 2:
            delta = SCORE["delta"]
        relationship = SCORE["relationship"]
        p_value = SCORE["p_value"]

        _fmt = aco._format_list
        x = gene_names + [
            phenotype, _fmt(group_names), _fmt(num_samples), _fmt(mean_score),
            delta, relationship, p_value]
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))
    if filestem:
        outhandle.close()

    # Write out other files.
    if not filestem:
        return

    jobs = []  # list of (fn, args, keywds)
    for x in itertools.product(phenotypes, range(M.nrow())):
        pheno_header, gene_i = x
        SCORE = gene_phenotype_scores.get((pheno_header, gene_i))
        if not SCORE:
            continue
        
        # Write the PRISM file.
        gene_id = aco.format_gene_name(M, None, gene_i)
        sample_names = M.col_names(arrayio.COL_ID)
        filename = aco._make_filename(
            M, gene_i, filestem, pheno_header, args.gene_header,
            "prism", "txt")
        x1 = (filename,
                SCORE["scores"], SCORE["phenotypes"], SCORE["group_names"])
        x = write_prism_file, x1, {}
        jobs.append(x)

        
        # Make a boxplot.
        filename = aco._make_filename(
            M, gene_i, filestem, pheno_header, args.gene_header,
            "boxplot", "png")
        pretty_gene = aco.pretty_gene_name(M, args.gene_header, gene_i)

        group_names = SCORE["group_names"]
        pheno2scores = {}
        for pheno, score in zip(SCORE["phenotypes"], SCORE["scores"]):
            if pheno not in pheno2scores:
                pheno2scores[pheno] = []
            pheno2scores[pheno].append(score)
        p_value = "p=%.2g" % SCORE["p_value"]
        x1 = (filename, group_names, pheno2scores)
        x2 = {
            "height" : 1600,
            "width" : 1600,
            "title" : pretty_gene,
            "subtitle" : p_value,
            "subtitle_col" : "#A60400",
            "subtitle_size" : 1.2,
            "subtitle_line" : 0.5,
            "ylab" : "Gene Expression",
            "mar_bottom" : args.box_mar_bottom,
            "mar_left" : args.box_mar_left,
            "mar_top" : 1.25,
            }
        x = boxplot.plot_boxplot, x1, x2
        jobs.append(x)
            
        # Make a waterfall plot.
        #filename = "%s%s.%s.waterfall.png" % (
        #    filestem, pheno_header, gene_id_h)
        filename = aco._make_filename(
            M, gene_i, filestem, pheno_header, args.gene_header,
            "waterfall", "png")
        pretty = aco.pretty_gene_name(M, args.gene_header, gene_i)
        x1 = (
            filename, SCORE["scores"], SCORE["phenotypes"],
            SCORE["group_names"], sample_names, SCORE["p_value"], pretty,
            args.water_mar_bottom, args.water_mar_left, args.water_mar_top,
            args.water_xlabel_off)
        x = plot_waterfall, x1, {}
        jobs.append(x)
    parallel.pyfun(jobs, num_procs=args.num_procs)

            
if __name__ == '__main__':
    #import profile
    #profile.run("main()")
    main()
