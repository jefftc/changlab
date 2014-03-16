#!/usr/bin/env python


# Functions:
# parse_phenotypes
# calc_association


def parse_phenotypes(phenotypes):
    # list of phenotypes.
    # e.g. ["STEM", "EMT"]
    # Return (potentially empty) list of phenotypes.
    return phenotypes


def calc_association(phenotypes, scores, expression_or_score):
    # Return a dictionary with keys:
    # scores               n-list of <float>
    # phenotypes           n-list of <string>
    # groups               n-list of <int>  [0, length(group_names)-1]
    # group_names          m-list of <string>  (unique list of pheno)
    # num_samples          dict of <group> : <int>
    # mean_score           dict of <group> : <float>
    # p_value              <float>
    # relationship         <string>
    from genomicode import jmath
    
    # Select only the samples with phenotype and score information.
    I1 = [i for (i, x) in enumerate(phenotypes) if x]
    I2 = [i for (i, x) in enumerate(scores) if x]
    I = sorted(set.intersection(set(I1), set(I2)))
    assert I, "No valid samples."

    phenotypes = [phenotypes[i] for i in I]
    scores = [float(scores[i]) for i in I]

    # Figure out the groupings.
    group_names = sorted({}.fromkeys(phenotypes))
    assert len(group_names) >= 2, "Need at least 2 groups."
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
    R = jmath.start_R()
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

    # Figure out the relationship.
    relationship = ""

    SCORE = {}
    SCORE["scores"] = scores
    SCORE["phenotypes"] = phenotypes
    SCORE["group_names"] = group_names
    SCORE["num_samples"] = num_samples
    SCORE["mean_score"] = mean_score
    SCORE["p_value"] = p_value
    SCORE["relationship"] = relationship
    return SCORE


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


def plot_boxplot(
    filename, scores, phenotypes, group_names, p_value, gene_id,
    mar_bottom, mar_left, mar_top):
    import os
    from genomicode import jmath
    from genomicode.jmath import R_fn, R_var, R_equals
    from genomicode import config

    xlabel_size = 1.0
    height = 1600
    width = 1600

    pheno2scores = {}
    for pheno, score in zip(phenotypes, scores):
        if pheno not in pheno2scores:
            pheno2scores[pheno] = []
        pheno2scores[pheno].append(score)

    R = jmath.start_R()
    path = config.changlab_Rlib
    plotlib = os.path.join(path, "plotlib.R")
    assert os.path.exists(plotlib), "I cannot find: %s" % plotlib
    R_fn("source", plotlib)

    #main = R_var("NA")
    main = gene_id
    sub = ""
    xlab = ""
    ylab = "Gene Expression"
    labels = group_names
    col = R_var("NULL")

    lwd = 2
    las = 3   # vertical labels
    at = R_var("NULL")
    if labels:
        at = range(1, len(labels)+1)
    cex_labels = 1.25*xlabel_size
    cex_legend = 1
    cex_lab = 1.5
    cex_sub = 1.5

    R_equals(labels, "labels")
    R_equals(at, "at")
    R("X <- list()")
    for i, n in enumerate(group_names):
        s = pheno2scores.get(n, [])
        R_equals(s, "s")
        R("X[[%d]] <- s" % (i+1))

    bm_type = "png16m"
    if filename.lower().endswith(".pdf"):
        bm_type = "pdfwrite"
    R_fn(
        "bitmap", filename, type=bm_type,
        height=height, width=width, units="px", res=300)
    
    # Set the margins.
    x = 5*mar_bottom, 5*mar_left, 4*mar_top, 2
    mar = [x+0.1 for x in x]
    R_fn("par", mar=mar, RETVAL="op")
        
    R_fn(
        "boxplot", R_var("X"), col=col, main="", xlab="", ylab="",
        axes=R_var("FALSE"), pch=19, cex=1, ylim=R_var("NULL"))
    # Make plot area solid white.
    jmath.R('usr <- par("usr")')
    jmath.R('rect(usr[1], usr[3], usr[2], usr[4], col="#FFFFFF")')
    R_fn(
        "boxplot", R_var("X"), col=col, main="", xlab="", ylab="",
        axes=R_var("FALSE"), pch=19, cex=1, ylim=R_var("NULL"),
        add=R_var("TRUE"))
    
    R_fn("box", lwd=lwd)
    R_fn(
        "axis", 1, lwd=lwd, labels=R_var("labels"),
        at=R_var("at"), las=las, **{ "cex.axis" : cex_labels })
    R_fn(
        "axis", 2, lwd=lwd, **{ "cex.axis" : 1.5 })
    R_fn(
        "title", main=main, sub=sub, xlab=xlab, ylab=ylab,
        **{ "cex.lab" : cex_lab, "cex.main" : 2.0, "cex.sub" : cex_sub })
    R("par(op)")
    R_fn("dev.off")



def plot_waterfall(
    filename, scores, phenotypes, group_names, sample_names, p_value, gene_id):
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
    mar_bottom = 1.0
    mar_left = 1.0


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
    cex_legend = 1
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
    x = 16*mar_bottom, 6*mar_left, 4, 2
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
    from genomicode import hashlib

    parser = argparse.ArgumentParser(
        description="Associate gene expression patterns with a "
        "categorical phenotype.")
    
    parser.add_argument(
        'expression_file',
        help='Either a gene expression file (GCT,CDT,PCL format) or gene set '
        'scores from score_geneset.py.')
    parser.add_argument('phenotype_file', help='Table of phenotypes.')
    
    group = parser.add_argument_group(title='Analysis')
    group.add_argument(
        '--phenotype', default=[], action='append',
        help='Header in the phenotype file.  Format: <header>')
    group.add_argument(
        '--gene', default=[], action='append',
        help='Comma separated name or ID of genes to analyze.  '
        'I will search for this gene in the annotations of the '
        'expression_file.  '
        'You can use this parameter multiple times to search more genes.')
    group.add_argument(
        '--geneset', default=[], action='append',
        help='Name of the geneset to analyze. To specify multiple gene sets, '
        'use this parameter multiple times.')
    
    group = parser.add_argument_group(title='Output')
    group.add_argument(
        '-o', dest='filestem', default=None,
        help='Prefix used to name files.  e.g. "myanalysis".')

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
    assert args.phenotype_file, (
        'Please specify a phenotype file.')
    assert os.path.exists(args.phenotype_file), "File not found: %s" % \
           args.phenotype_file
    
    assert args.phenotype, 'Please specify the phenotype to analyze.'
    assert args.gene or args.geneset, 'Please specify a gene or gene set.'
    assert not (args.gene and args.geneset), (
        'Please specify either a gene or a gene set, not both.')

    assert args.box_mar_bottom > 0 and args.box_mar_bottom < 10
    assert args.box_mar_left > 0 and args.box_mar_left < 10
    assert args.box_mar_top > 0 and args.box_mar_top < 10
    ## assert args.km_title_size > 0 and args.km_title_size < 10
    ## assert args.km_mar_title > 0 and args.km_mar_title < 10
    ## assert args.km_subtitle_size > 0 and args.km_subtitle_size < 10
    ## assert args.km_mar_subtitle > 0 and args.km_mar_subtitle < 10
    ## assert args.km_legend_size > 0 and args.km_legend_size < 10
    
    # Clean up the input.
    genes = aco.parse_genes(args.gene)
    gene_sets = aco.parse_gene_sets(args.geneset)
    phenotypes = parse_phenotypes(args.phenotype)
    filestem = aco.parse_filestem(args.filestem)

    # Read the input files.
    M = aco.read_expression_or_geneset_scores(
        genes, gene_sets, args.expression_file)
    x = aco.read_clinical_annotations(M, args.phenotype_file)
    M, clinical_annots = x

    # Make sure at least one of the phenotypes are in the clinical
    # annotations.
    phenotypes = [x for x in phenotypes if x in clinical_annots]
    assert phenotypes, "No phenotypes found."

    # Select the genes or gene sets of interest.
    x = genes or gene_sets
    M = M.matrix(row=x)
    assert M.nrow(), "I could not find any of the genes or gene sets."

    # Calculate the association of each gene and each phenotype.
    expression_or_score = "Expression"
    if gene_sets:
        expression_or_score = "Score"

    # (header, gene_index) -> returned from calc_association
    gene_phenotype_scores = {}  
    for x in itertools.product(phenotypes, range(M.nrow())):
        pheno_header, i = x
        phenotype = clinical_annots[pheno_header]
        scores = M.value(i, None)

        x = calc_association(phenotype, scores, expression_or_score)
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
        outhandle = open("%sstats.txt" % filestem, 'w')

    # Figure out the header for the table.
    header = M.row_names() + [
        "Phenotype", "Groups", "Num Samples", "Average Expression",
        "Relationship", "p-value"]
    print >>outhandle, "\t".join(header)

    # Write out each row of the table.
    for x in itertools.product(phenotypes, range(M.nrow())):
        pheno_header, gene_i = x

        SCORE = gene_phenotype_scores[(pheno_header, gene_i)]

        gene_names = [M.row_names(x)[gene_i] for x in M.row_names()]
        phenotype = pheno_header
        group_names = SCORE["group_names"]
        I = range(len(group_names))
        num_samples = [SCORE["num_samples"][x] for x in I]
        mean_score = [SCORE["mean_score"][x] for x in I]
        relationship = SCORE["relationship"]
        p_value = SCORE["p_value"]

        _fmt = aco._format_list
        x = gene_names + [
            phenotype, _fmt(group_names), _fmt(num_samples), _fmt(mean_score),
            relationship, p_value]
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))


        # Write out plots.
        # Better way to pick gene ID.
        gene_id = aco.get_gene_name(M, gene_i)
        gene_id_h = hashlib.hash_var(gene_id)
        sample_names = M.col_names(arrayio.COL_ID)

        if filestem:
            filename = "%s%s.%s.prism.txt" % (
                filestem, pheno_header, gene_id_h)
            write_prism_file(
                filename, SCORE["scores"], SCORE["phenotypes"],
                SCORE["group_names"])

        if filestem:
            filename = "%s%s.%s.boxplot.png" % (
                filestem, pheno_header, gene_id_h)
            plot_boxplot(
                filename, SCORE["scores"], SCORE["phenotypes"],
                SCORE["group_names"], SCORE["p_value"], gene_id,
                args.box_mar_bottom, args.box_mar_left, args.box_mar_top,
                )
            
        if filestem:
            filename = "%s%s.%s.waterfall.png" % (
                filestem, pheno_header, gene_id_h)
            plot_waterfall(
                filename, SCORE["scores"], SCORE["phenotypes"],
                SCORE["group_names"], sample_names, SCORE["p_value"], gene_id)

            
if __name__ == '__main__':
    main()
