#!/usr/bin/env python


# Functions:
# parse_genes
# parse_gene_sets
# parse_cutoffs
# parse_zscore
# parse_outcomes
# parse_filestem
#
# read_gene_expression
# read_geneset_scores
# read_expression_or_geneset_scores
# read_clinical_annotations
# align_matrix_with_clinical_data
#
# calc_association
# calc_km
# 
# discretize_scores
# discretize_by_value
# discretize_by_zscore
# make_group_names
# make_zscore_names
#
# write_prism_file
# plot_km
# plot_groups
#
# start_and_init_R
# colortuple2hex


def parse_genes(genes):
    # genes is a list of comma-separated gene names or IDs.
    # e.g. ["E2F1,E2F3", "MYC"]
    # Return (potentially empty) list of genes.
    clean = []
    for x in genes:
        clean.extend(x.split(','))
    return clean


def parse_gene_sets(geneset):
    # geneset is a list of gene sets.  Return it unchanged.
    return geneset


def parse_cutoffs(cutoffs):
    # Comma-separated list of breakpoints, e.g. 0.25,0.50.  If None,
    # then default to just 0.50.
    # Return list, e.g. [0.25, 0.50], that:
    # - is sorted
    # - has no duplicates
    # - all numbers > 0 and < 1 (no 0 or 1)
    if not cutoffs:
        return None
    
    cutoffs = cutoffs.split(',')
    cutoffs = [float(x) for x in cutoffs]
    cutoffs.sort()
    # Make sure all numbers between 0 and 1.
    for x in cutoffs:
        assert x >= 0 and x <= 1.0, "cutoff should be between 0 and 1: %s" % x
    assert cutoffs, "no cutoffs"

    # Remove duplicates and remove any numbers too close to 0 and 1.
    DELTA = 0.001
    i = 1
    while i < len(cutoffs):
        if cutoffs[i]-cutoffs[i-1] < DELTA:
            del cutoffs[i]
        else:
            i += 1
    assert cutoffs, "invalid cutoffs"
    if cutoffs[0] < DELTA:
        del cutoffs[0]
    assert cutoffs, "invalid cutoffs"
    if cutoffs[-1] > 1.0-DELTA:
        del cutoffs[-1]
    assert cutoffs, "invalid cutoffs"
    
    return cutoffs


def parse_zscore(zscore):
    # Floating point number, e.g. 1.0, or None.
    if not zscore:
        return None

    zscore = float(zscore)
    assert zscore > 0.10 and zscore < 100
    return zscore


def parse_outcomes(outcomes):
    # List of <time_header>,<dead_header>.  Return a list of tuples:
    # (<time_header>, <dead_header>).
    clean = []
    for x in outcomes:
        x = [x.strip() for x in x.split(",")]
        assert len(x) == 2
        clean.append(x)

    # Make sure there are no duplicate outcomes.
    seen = {}
    for (x1, x2) in clean:
        assert x1 not in seen, "duplicate outcome: %s" % x1
        seen[x1] = 1
        
    return clean


def parse_filestem(filestem):
    # Return an empty string, or a filestem with a '.' at the end.
    if filestem is None:
        filestem = ""
    elif not filestem.endswith("."):
        filestem += "."
    return filestem


def read_gene_expression(filename):
    import os
    import arrayio
    
    assert os.path.exists(filename)
    M = arrayio.read(filename)
    return M


def read_geneset_scores(filename):
    # Read the output from score_geneset.py and return a Matrix
    # object.
    import os
    from genomicode import jmath
    from genomicode import filelib
    from genomicode import Matrix
    
    assert os.path.exists(filename)
    matrix = [x for x in filelib.read_cols(filename)]
    matrix = jmath.transpose(matrix)
    # BUG: Need more checks on size and format of matrix.
    col_names = {}
    col_names['_SAMPLE_NAME'] = matrix[1][1:]
    row_names = {}
    row_names['geneset'] = []
    data = []
    for line in matrix[2:]:
        single_data = [float(i) for i in line[1:]]
        data.append(single_data)
        row_names['geneset'].append(line[0])
    M = Matrix.InMemoryMatrix(data, row_names=row_names, col_names=col_names)
    return M


def read_expression_or_geneset_scores(genes, gene_sets, filename):
    assert not (genes and gene_sets)
    assert genes or gene_sets
    
    if genes:
        # If genes were specified, then the input should be a gene
        # expression data file.
        M = read_gene_expression(filename)
    else:
        M = read_geneset_scores(filename)
    return M


def read_clinical_annotations(M, filename):
    # Read the clinical annotations.
    from genomicode import genesetlib
    
    clinical_annots = {}
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, values = x
        clinical_annots[name] = values

    # Align the gene scores with the clinical annotations.
    x = align_matrix_with_clinical_data(M, clinical_annots)
    M, clinical_annots = x
    
    return M, clinical_annots


def align_matrix_with_clinical_data(M, clinical_dict):
    from genomicode import matrixlib

    assert clinical_dict, "No clinical data."

    sample_name = M._col_names['_SAMPLE_NAME']

    # Figure out the header in the clinical data that contains these
    # samples.
    name2count = {}
    for name, values in clinical_dict.iteritems():
        count = len(set(sample_name).intersection(values))
        name2count[name] = count
    best_name = best_count = None
    for name, count in name2count.iteritems():
        if best_count is None or count > best_count:
            best_name, best_count = name, count
    assert best_count > 0, \
           "I could not align the matrix with the clinical data."
    clinical_name = clinical_dict[best_name]

    # From the clinical data, only keep the samples that occur in the
    # matrix.
    I = [i for (i, name) in enumerate(clinical_name) if name in sample_name]
    clinical_clean = {}
    for name, values in clinical_dict.iteritems():
        values_clean = [values[i] for i in I]
        clinical_clean[name] = values_clean
    clinical_dict = clinical_clean

    # Realign the matrix to match the clinical data.
    x = matrixlib.align_cols_to_annot(
        M, clinical_dict[best_name], reorder_MATRIX=True)
    M, colnames = x
    assert colnames == clinical_dict[best_name]

    return M, clinical_dict


def discretize_scores(scores, cutoffs, zscore, expression_or_score):
    # Discretize into groups based on cutoffs or zscore.  Return tuple
    # of (group_names, groups).  group_names is a list of the names of
    # each group.  groups is a list of integers from 0 to
    # len(group_names)-1, parallel to scores.
    if cutoffs:
        groups = discretize_by_value(scores, cutoffs)
        # groups should range from [0, len(cutoffs)+1].  So if there's
        # one cutoff (e.g. 0.50), groups will be 0 and 1.
        for x in groups:
            assert x >= 0 and x <= len(cutoffs)+1
        group_names = make_group_names(cutoffs, expression_or_score)
    elif zscore:
        groups = discretize_by_zscore(scores, zscore)
        group_names = make_zscore_names(expression_or_score)
    else:
        raise AssertionError
            
    # Make sure there's at least 2 groups.
    assert len({}.fromkeys(groups)) >= 2, "Need at least 2 groups."
    return group_names, groups


def calc_association(survival, dead, scores, cutoffs, zscore,
                     expression_or_score):
    # Return a dictionary with keys:
    # survival             list of <float>
    # dead                 list of <int>
    # scores               list of <float>
    # groups               list of <int>  [0, length(group_names)-1]
    # group_names          list of <string>
    # p_value              <float>
    # num_samples          dict of <group> : <int>
    # mean_score           dict of <group> : <float>
    # surv50               dict of <group> : <float> or None
    # surv90               dict of <group> : <float> or None
    # hi_score_short_surv  <boolean> or None (no difference in surv)
    # relationship         <string>
    from genomicode import jmath
    
    # Select only the samples with both survival, dead, and score
    # information.
    I1 = [i for (i, x) in enumerate(survival) if x]
    I2 = [i for (i, x) in enumerate(dead) if x]
    I3 = [i for (i, x) in enumerate(scores) if x]
    I = sorted(set.intersection(set(I1), set(I2), set(I3)))
    assert I, "No valid samples."

    survival = [float(survival[i]) for i in I]
    dead = [int(dead[i]) for i in I]
    scores = [scores[i] for i in I]

    # Figure out the groupings.
    x = discretize_scores(scores, cutoffs, zscore, expression_or_score)
    group_names, groups = x

    # Calculate the KM model.
    surv = calc_km(survival, dead, groups)

    # Clean up the surv dictionary.  If some groups are missing, some
    # of the members will be missing values.  Fix this.
    for i in range(len(group_names)):
        if i not in surv["num_samples"]:
            surv["num_samples"][i] = 0
        if i not in surv["surv50"]:
            surv["surv50"][i] = None
        if i not in surv["surv90"]:
            surv["surv90"][i] = None

    # Add extra data to the survival dictionary.
    surv["survival"] = survival
    surv["dead"] = dead
    surv["scores"] = scores
    surv["groups"] = groups
    surv["group_names"] = group_names

    # Calculate the mean scores for each group.  If a group is empty,
    # then the mean score is None.
    mean_score = {}
    for group in range(len(group_names)):
        s = [s for (s, g) in zip(scores, groups) if g == group]
        m = None
        if s:
            m = jmath.mean(s)
        mean_score[group] = m
    surv["mean_score"] = mean_score

    # Figure out relationship.
    MAX_SURV = 1E10
    # Compare the time to 50% survival for the low and high scoring
    # groups.
    # ASSUMPTION: lowest group has low scores, while highest group has
    # high scores.
    surv_low = surv["surv50"][min(groups)]    # low score
    surv_high = surv["surv50"][max(groups)]   # high score
    # If neither groups drop to 50% survival, compare the time to 90%
    # survival.
    if surv_low is None and surv_high is None:
        surv_low = surv["surv90"][min(groups)]
        surv_high = surv["surv90"][max(groups)]
    if surv_high is None:
        surv_high = MAX_SURV
    if surv_low is None:
        surv_low = MAX_SURV
    assert surv_low <= MAX_SURV and surv_high <= MAX_SURV
    hi_score_short_surv = None
    if surv_high < surv_low:
        hi_score_short_surv = True
    elif surv_high > surv_low:
        hi_score_short_surv = False
    surv["hi_score_short_surv"] = hi_score_short_surv

    relationship = ""
    if hi_score_short_surv:
        relationship = "High %s has shorter time to outcome." % \
                       expression_or_score.lower()
    elif hi_score_short_surv is not None:
        relationship = "Low %s has shorter time to outcome." % \
                       expression_or_score.lower()
    surv["relationship"] = relationship

    return surv
    

def calc_km(survival, dead, group):
    # Return a dictionary with keys:
    # p_value         <float>
    # num_samples     dict of <group> : <int>
    # surv50          dict of <group> : <float> or None
    # surv90          dict of <group> : <float> or None
    from genomicode import jmath

    assert len(survival) == len(dead)
    assert len(survival) == len(group)
    
    R = start_and_init_R()

    # Suppress a potential warning message:
    # Warning messages:
    # 1: In fitter(X, Y, strats, offset, init, control, weights = weights,  :
    #   Loglik converged before variable  1,2 ; beta may be infinite.
    # 2: In min(abs(s$surv - 0.9)) :
    #   no non-missing arguments to min; returning Inf
    # 3: In min(abs(s$surv - 0.5)) :
    #   no non-missing arguments to min; returning Inf
    jmath.R_equals(survival, 'survival')
    jmath.R_equals(dead, 'dead')
    jmath.R_equals(group, 'group')
    R('ow <- options("warn")')
    R('options(warn=-1)')
    R('x <- calc.km.multi(survival, dead, group)')
    R('options(ow)')

    c = R['x']
    p_value = c.rx2('p.value')[0]
    #hazard_ratio = list(c.rx2('hr'))

    unique_group = [x for x in {}.fromkeys(group)]
    num_samples = {}
    for x in unique_group:
        num_samples[x] = c.rx2('num.samples').rx2(str(x))[0]
    
    surv50 = {}
    surv90 = {}
    for x in unique_group:
        x1 = c.rx2('surv').rx2(str(x)).rx2('surv.50')[0]
        x2 = c.rx2('surv').rx2(str(x)).rx2('surv.90')[0]
        if str(x1) == "NA":
            x1 = None
        if str(x2) == "NA":
            x2 = None
        surv50[x] = x1
        surv90[x] = x2

    SURV = {
        "p_value" : p_value,
        "num_samples" : num_samples,
        "surv50" : surv50,
        "surv90" : surv90,
        }
    return SURV
    

def discretize_by_value(values, breakpoints):
    # Return a list that specifies the group for each member of
    # values.  Groups are specified by numbers starting from 0.  Group
    # 0 are the values lower than the first breakpoint.
    from genomicode import jmath

    R = start_and_init_R()
    jmath.R_equals(values, 'F')
    jmath.R_equals(breakpoints, 'cutoffs')
    R('group <- group.by.value(F, cutoffs)')
    groups = list(R['group'])
    assert len(groups) == len(values)
    return groups


def discretize_by_zscore(values, zscore):
    # Return a list that specifies the group for each member of
    # values.  Group 0 are values <= -zscore, 1 are -zscore < value <
    # zscore, and 2 are values >= zscore.
    from genomicode import jmath

    R = start_and_init_R()

    jmath.R_equals(range(len(values)), 'x')
    jmath.R_equals(values, 'y')
    R('M <- find.outliers(x, y, z.cutoff=1)')
    M = R['M']
    z = list(M.rx2('z'))

    groups = [1] * len(z)
    for i in range(len(z)):
        if z[i] <= -zscore:
            groups[i] = 0
        elif z[i] >= zscore:
            groups[i] = 2
    return groups


def make_group_names(cutoffs, expression_or_score):
    assert len(cutoffs) >= 1
    if len(cutoffs) == 1:
        # Just low vs high.
        # [0, <cutoff>, 1]
        x = ["Low %s" % expression_or_score, "High %s" % expression_or_score]
        return x

    cutoffs = [0] + cutoffs + [1]  # for convenience
    names = []
    for i in range(len(cutoffs)-1):
        x1 = "Middle"
        if i == 0:
            x1 = "Lowest"
        elif i == len(cutoffs)-2:
            x1 = "Highest"
        x2 = "%d%% - %d%%" % ((cutoffs[i])*100, (cutoffs[i+1])*100)
        x = "%s %s %s" % (x1, x2, expression_or_score)
        names.append(x)
    return names


def make_zscore_names(expression_or_score):
    names = [
        "Low Outlier",
        "Normal %s" % expression_or_score,
        "High Outlier",
        ]
    return names


def write_prism_file(filename, survival, dead, group_indexes, group_names):
    from genomicode import jmath
    
    R = start_and_init_R()
    jmath.R_equals(filename, 'filename')
    group = [group_names[i] for i in group_indexes]
    jmath.R_equals(survival, 'survival')
    jmath.R_equals(dead, 'dead')
    jmath.R_equals(group, 'group')
    R('write.km.prism.multi(filename, survival, dead, group)')


def plot_km(filename, survival, dead, group_indexes,
            p_value, gene_id, group_names, xlab, ylab, title):
    from genomicode import colorlib
    from genomicode.jmath import R_equals, R_fn
    
    R = start_and_init_R()
    
    # Set the colors.
    assert len(group_names) >= 2
    colors = ['#1533AD', '#FFB300']
    if len(group_names) > 2:
        x = colorlib.bild_colors(len(group_names))
        x = [colortuple2hex(*x) for x in x]
        colors = x
    # R command:
    # col <- list("<name>"="<color>", ...)
    cmd = []
    for i in range(len(group_names)):
        x = '"%s"="%s"' % (group_names[i], colors[i])
        cmd.append(x)
    cmd = "col <- list(%s)" % (", ".join(cmd))
    R(cmd)

    group = [group_names[i] for i in group_indexes]
    R_equals(survival, 'survival')
    R_equals(dead, 'dead')
    R_equals(group, 'group')
    R_equals('p_value=%.2g' % p_value, 'sub')
    
    R('xlab <- ""')
    R('ylab <- ""')
    R_equals(gene_id, 'title')
    if xlab:
        R_equals(xlab, 'xlab')
    if ylab:
        R_equals(ylab, 'ylab')
    if title:
        R_equals(title, 'title')
        
    R_fn(
        "bitmap", file=filename, type="png256",
         height=1600, width=1600, units="px", res=300)
    # Suppress warning message.  See calc_km.
    R('ow <- options("warn")')
    R('options(warn=-1)')
    R('plot.km.multi(survival, dead, group, '
      'col=col, main=title, sub=sub, xlab=xlab, ylab=ylab)')
    R('options(ow)')
    R('dev.off()')
    

def plot_groups(filename, scores, group_names, groups):
    from genomicode import colorlib
    from genomicode.jmath import R_fn, R_equals, R_var
    
    start_and_init_R()
    #R = start_and_init_R()

    # Set the colors.
    assert len(group_names) >= 2
    colors = ['#1533AD', '#FFB300']
    if len(group_names) > 2:
        x = colorlib.bild_colors(len(group_names))
        x = [colortuple2hex(*x) for x in x]
        colors = x
    col = [colors[x] for x in groups]
    R_equals(col, 'col')
    
    mar = [x+0.1 for x in [5, 6, 4, 2]]
    
    x = range(len(scores))
    y = scores
    R_equals(x, 'x')
    R_equals(y, 'y')

    R_fn(
        "bitmap", filename, type="png256", height=1600, width=1600,
        units="px", res=300)
    R_fn("par", mar=mar, RETVAL="op")
    R_fn(
        "plot", R_var('x'), R_var('y'), type="n", axes=R_var("FALSE"),
        xlab="", ylab="")
    R_fn("par", "usr", RETVAL="usr")
    R_fn(
        "rect", R_var("usr[1]"), R_var("usr[3]"), R_var("usr[2]"),
        R_var("usr[4]"), col="#FFFFFF")
    R_fn("points", R_var('x'), R_var('y'), pch=19, cex=1, col=R_var('col'))
    R_fn("box", lwd=1.5)
    R_fn(
        "axis", 1, lwd=1.5, labels=R_var("FALSE"), tick=R_var("FALSE"),
        **{"cex.axis" : 1.25})
    R_fn("axis", 2, lwd=1.5, **{"cex.axis" : 1.25})
    R_fn("par", mgp=[1.0, 1.5, 0], RETVAL="op2")
    R_fn("title", xlab="Sample", **{"cex.lab" : 1.5})
    R_fn("par", R_var("op2"))
    R_fn("title", main="", xlab="", ylab="Gene Expression", sub="",
         **{"cex.lab":1.5, "cex.sub":1, "col.sub":"#A60400", "cex.main":1.0})
    R_fn(
        "legend", "bottomleft", legend=group_names, fill=colors, inset=0.05,
        cex=1.25, **{"box.lwd":1.5})
    R_fn("par", R_var("op"))
    R_fn("dev.off")
    

GLOBAL_R = None
def start_and_init_R():
    global GLOBAL_R
    import os
    from genomicode import jmath
    from genomicode import config

    if GLOBAL_R is None:
        assert os.path.exists(config.changlab_Rlib)
        km_lib = os.path.join(config.changlab_Rlib, "kaplanmeierlib.R")
        stat_lib = os.path.join(config.changlab_Rlib, "statlib.R")
        assert os.path.exists(km_lib), "File not found: %s" % km_lib
        assert os.path.exists(stat_lib), "File not found: %s" % stat_lib

        R = jmath.start_R()
        R('require(splines, quietly=TRUE)')
        R('source("%s")' % km_lib)
        R('source("%s")' % stat_lib)
        GLOBAL_R = R
    return GLOBAL_R


def colortuple2hex(R, G, B):
    R = hex(int(R * 255))[2:]
    G = hex(int(G * 255))[2:]
    B = hex(int(B * 255))[2:]
    color_hex = '#' + R + G + B
    return color_hex


def _format_list(x):
    for i in range(len(x)):
        if x[i] is None:
            x[i] = ""
    x = map(str, x)
    return ";".join(x)


def main():
    import os
    import sys
    import argparse
    import itertools
    from genomicode import hashlib

    parser = argparse.ArgumentParser(
        description='Associate gene expression patterns with outcomes.')
    
    parser.add_argument(
        'expression_file',
        help='Either a gene expression file (GCT,CDT,PCL format) or gene set '
        'scores from score_geneset.py.')
    parser.add_argument('outcome_file', help='Table of clinical annotations.')
    
    group = parser.add_argument_group(title='Analysis')
    group.add_argument(
        '--outcome', default=[], action='append',
        help='Where to find the outcome information in the clinical '
        'annotation file.  To analyze more than one outcome, use this '
        'parameter multiple times.  Format: <time_header>,<dead_header>')
    group.add_argument(
        '--gene', default=[], action='append',
        help='Name or ID of gene to analyze.  I will search for this gene '
        'in the annotations of the expression_file.  '
        'To analyze more than one gene, use this parameter multiple times.')
    group.add_argument(
        '--geneset', default=[], action='append',
        help='Name of the geneset to analyze. To specify multiple gene sets, '
        'use this parameter multiple times.')
    
    group = parser.add_argument_group(
        title='Discretization', description='We can discretize the '
        'expression data using two different strategies.  '
        'The first approach is to set a simple cutoff based on the '
        'relative magnitude of the expression.  For example, a simple '
        'cutoff of 0.50 would separate the samples with lowest 50% '
        'expression from those with highest 50%.  '
        'The second approach looks for outliers in gene expression '
        'based on z-scores.  So a z-score cutoff of 1.0 would separate '
        'the samples into groups with z <= -1, -1 < z < 1, and z >= 1.  '
        'For the z-score strategy, it is possible that no outliers will '
        'be detected with a specific z-score cutoff.'
        )
    group.add_argument(
        '--cutoff', default=None,
        help='Comma-separated list of breakpoints (between 0 and 1), '
        'e.g. 0.25,0.50,0.75.  Default is to use cutoff of 0.50.  '
        'I will use this strategy unless a --zscore is given.')
    group.add_argument(
        '--zscore', default=None, help='Z-score cutoff, e.g. 1.0.')

    group = parser.add_argument_group(title='Ouput')
    group.add_argument(
        '-o', dest='filestem', default=None,
        help='Prefix used to name files.  e.g. "myanalysis".')
    group.add_argument(
        '--write_prism', action='store_true',
        default=False,
        help='Write a text file that can be imported into GraphPad Prism.')
    group.add_argument(
        '--plot_km', action='store_true', default=False,
        help='Write PNG-formatted Kaplan-Meier plot.')
    group.add_argument(
        '--plot_groups', action='store_true', default=False,
        help='Make a scatter plot showing the groups for each sample '
        '(PNG format).')
    group.add_argument(
        '--xlab', default=None, 
        help='x-axis label for the Kaplan-Meier plot.')
    group.add_argument(
        '--ylab', default=None, 
        help='y-axis label for the Kaplan-Meier plot.')
    group.add_argument(
        '--title', default=None, help='Title for the Kaplan-Meier plot.')

    args = parser.parse_args()

    # Check inputs.
    assert args.expression_file, (
        'Please specify a gene expression or gene set score file.')
    assert os.path.exists(args.expression_file), "File not found: %s" % \
           args.expression_file
    assert args.outcome_file, (
        'Please specify a clinical outcomes file.')
    assert os.path.exists(args.outcome_file), "File not found: %s" % \
           args.outcome_file
    
    assert args.outcome, 'Please specify the clinical outcomes to analyze.'
    assert args.gene or args.geneset, 'Please specify a gene or gene set.'
    assert not (args.gene and args.geneset), (
        'Please specify either a gene or a gene set, not both.')

    assert not (args.cutoff and args.zscore)
    if not args.cutoff and not args.zscore:
        args.cutoff = "0.50"

    # Clean up the input.
    genes = parse_genes(args.gene)
    gene_sets = parse_gene_sets(args.geneset)
    cutoffs = parse_cutoffs(args.cutoff)
    zscore = parse_zscore(args.zscore)
    outcomes = parse_outcomes(args.outcome)
    filestem = parse_filestem(args.filestem)

    # Read the input files.
    M = read_expression_or_geneset_scores(
        genes, gene_sets, args.expression_file)
    x = read_clinical_annotations(M, args.outcome_file)
    M, clinical_annots = x

    # Make sure each of the outcomes are in the clinical annotations.
    for x1, x2 in outcomes:
        assert x1 in clinical_annots, "Missing clinical annotation: %s" % x1
        assert x2 in clinical_annots, "Missing clinical annotation: %s" % x2

    # Select the genes or gene sets of interest.
    x = genes or gene_sets
    M = M.matrix(row=x)
    assert M.nrow(), "I could not find any of the genes or gene sets."

    # Calculate the association of each gene and each outcome.
    # (time_header, dead_header, gene_index) -> returned from calc_association
    
    expression_or_score = "Expression"
    if gene_sets:
        expression_or_score = "Score"

    gene_outcome_scores = {}  
    for x in itertools.product(outcomes, range(M.nrow())):
        (time_header, dead_header), i = x
        survival = clinical_annots[time_header]
        dead = clinical_annots[dead_header]
        scores = M.value(i, None)

        x = calc_association(
            survival, dead, scores, cutoffs, zscore, expression_or_score)
        gene_outcome_scores[(time_header, dead_header, i)] = x

    # Files generated:
    # <filestem>.stats.txt           Or to STDOUT if no <filestem> given.
    # <filestem>.<outcome>.<gene_id>.km.png      K-M plot.
    # <filestem>.<outcome>.<gene_id>.prism.txt   Prism format for K-M analysis.
    # <filestem>.<outcome>.<gene_id>.groups.png  Group for each sample.

    # Write the output in a table with headers:
    # <headers>            # From the expression or gene set file.
    # Outcome
    # Groups               # one for each group
    # Num Samples          # one for each group, separated by semicolon
    # Average Expression   # one for each group, separated by semicolon
    # 90% Survival         # one for each group, separated by semicolon
    # 50% Survival         # one for each group, separated by semicolon
    # Relationship
    # p-value

    outhandle = sys.stdout
    if filestem:
        outhandle = open("%sstats.txt" % filestem, 'w')

    # Figure out the header for the table.
    header = M.row_names() + [
        "Outcome", "Groups", "Num Samples", "Average Expression",
        "90% Survival", "50% Survival", "Relationship", "p-value"]
    print >>outhandle, "\t".join(header)

    # Write out each row of the table.
    for x in itertools.product(outcomes, range(M.nrow())):
        (time_header, dead_header), gene_i = x

        SURV = gene_outcome_scores[(time_header, dead_header, gene_i)]

        gene_names = [M.row_names(x)[gene_i] for x in M.row_names()]
        outcome = time_header
        group_names = SURV["group_names"]
        I = range(len(group_names))
        num_samples = [SURV["num_samples"][x] for x in I]
        mean_score = [SURV["mean_score"][x] for x in I]
        surv90 = [SURV["surv90"][x] for x in I]
        surv50 = [SURV["surv50"][x] for x in I]
        relationship = SURV["relationship"]
        p_value = SURV["p_value"]

        _fmt = _format_list
        x = gene_names + [
            outcome, _fmt(group_names), _fmt(num_samples), _fmt(mean_score),
            _fmt(surv90), _fmt(surv50), relationship, p_value]
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))


        # Write out Prism, Kaplan-Meier curves, etc.
        gene_id = M.row_names(M.row_names()[0])[gene_i]
        gene_id_h = hashlib.hash_var(gene_id)

        if args.write_prism:
            filename = "%s%s.%s.prism.txt" % (
                filestem, time_header, gene_id_h)
            write_prism_file(
                filename, SURV["survival"], SURV["dead"], SURV["groups"],
                SURV["group_names"])
        if args.plot_km:
            filename = "%s%s.%s.km.png" % (
                filestem, time_header, gene_id_h)
            plot_km(
                filename, SURV["survival"], SURV["dead"], SURV["groups"],
                SURV["p_value"], gene_id, SURV["group_names"], 
                args.xlab, args.ylab, args.title)
        if args.plot_groups:
            filename = "%s%s.%s.groups.png" % (
                filestem, time_header, gene_id_h)
            plot_groups(
                filename, SURV["scores"], SURV["group_names"], SURV["groups"])

            
if __name__ == '__main__':
    main()
