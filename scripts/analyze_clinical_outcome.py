#!/usr/bin/env python


# Functions:
# parse_genes
# parse_gene_sets
# parse_cutoffs
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
# group_by_value
# make_group_names
#
# write_prism_file
# plot_km
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
    # Return list of [0.25, 0.50].  List is guaranteed to be sorted
    # with no duplicates.  Numbers indicate internal breakpoints, so 0
    # and 1 will be removed.  All numbers will be > 0 and < 1.  
    
    cutoffs = cutoffs or "0.50"
    
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


def calc_association(survival, dead, scores, cutoffs):
    # Return a dictionary with keys:
    # survival             list of <float>
    # dead                 list of <int>
    # group                list of <int>
    # p_value              <float>
    # num_samples          dict of <group> : <int>
    # mean_score           dict of <group> : <float>
    # surv50               dict of <group> : <float> or None
    # surv90               dict of <group> : <float> or None
    # hi_score_short_surv  <boolean> or None (both the same)

    # <group> will be converted to strings.
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
    groups = group_by_value(scores, cutoffs)
    # groups should range from 0 to len(cutoffs)+1.  So if there's one
    # cutoff (e.g. 0.50), groups will be 0 and 1.
    for x in groups:
        assert x >= 0 and x <= len(cutoffs)+1
    # Make sure there's at least 2 groups.
    assert len({}.fromkeys(groups)) >= 2, "Need at least 2 groups."

    surv = calc_km(survival, dead, groups)
    surv["survival"] = survival
    surv["dead"] = dead
    surv["group"] = groups

    # Calculate the mean scores for each group.  If a group is empty,
    # then the mean score is None.
    mean_score = {}
    unique_groups = sorted({}.fromkeys(groups))
    for group in unique_groups:
        s = [s for (s, g) in zip(scores, groups) if g == group]
        m = None
        if s:
            m = jmath.mean(s)
        mean_score[str(group)] = m
    surv["mean_score"] = mean_score

    # Figure out relationship.
    MAX_SURV = 1E10
    # Compare the time to 50% survival for the low and high scoring
    # groups.
    surv_low = surv["surv50"][str(unique_groups[0])]    # low score
    surv_high = surv["surv50"][str(unique_groups[-1])]  # high score
    # If neither groups drop to 50% survival, compare the time to 90%
    # survival.
    if surv_low is None and surv_high is None:
        surv_low = surv["surv90"][str(unique_groups[0])]
        surv_high = surv["surv90"][str(unique_groups[-1])]
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
    
    return surv
    

def calc_km(survival, dead, group):
    # Return a dictionary with keys:
    # p_value         <float>
    # num_samples     dict of <group> : <int>
    # surv50          dict of <group> : <float> or None
    # surv90          dict of <group> : <float> or None
    # <group> will be converted to strings.
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

    unique_group = [str(x) for x in {}.fromkeys(group)]
    num_samples = {}
    for x in unique_group:
        num_samples[x] = c.rx2('num.samples').rx2(x)[0]
    
    surv50 = {}
    surv90 = {}
    for x in unique_group:
        x1 = c.rx2('surv').rx2(x).rx2('surv.50')[0]
        x2 = c.rx2('surv').rx2(x).rx2('surv.90')[0]
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
    

def group_by_value(values, breakpoints):
    # Return a list that specifies the group for each member of
    # values.  Groups are specified by numbers starting from 0.  Group
    # 0 are the values lower than the first breakpoint.
    from genomicode import jmath

    R = start_and_init_R()
    jmath.R_equals(values, 'F')
    jmath.R_equals(breakpoints, 'cutoffs')
    R('group <- group.by.value(F, cutoffs)')
    group = list(R['group'])
    assert len(group) == len(values)
    return group


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
    from genomicode import jmath
    from genomicode import colorlib
    
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
    jmath.R_equals(survival, 'survival')
    jmath.R_equals(dead, 'dead')
    jmath.R_equals(group, 'group')
    jmath.R_equals(filename, 'filename')
    jmath.R_equals('p_value=%.2g' % p_value, 'sub')
    
    R('xlab <- ""')
    R('ylab <- ""')
    jmath.R_equals(gene_id, 'title')
    if xlab:
        jmath.R_equals(xlab, 'xlab')
    if ylab:
        jmath.R_equals(ylab, 'ylab')
    if title:
        jmath.R_equals(title, 'title')
        
    jmath.R_fn(
        "bitmap", file=jmath.R_var("filename"), type="png256",
         height=1600, width=1600, units="px", res=300)
    # Suppress warning message.  See calc_km.
    R('ow <- options("warn")')
    R('options(warn=-1)')
    R('plot.km.multi(survival, dead, group, '
      'col=col, main=title, sub=sub, xlab=xlab, ylab=ylab)')
    R('options(ow)')
    R('dev.off()')
    

GLOBAL_R = None
def start_and_init_R():
    global GLOBAL_R
    import os
    from genomicode import jmath
    from genomicode import config

    if GLOBAL_R is None:
        kaplanmeierlib_path = config.kaplanmeierlib
        assert os.path.exists(kaplanmeierlib_path), (
            "can not find the kaplanmeierlib script %s" % kaplanmeierlib_path)
    
        R = jmath.start_R()
        R('require(splines, quietly=TRUE)')
        R('source("%s")' % config.kaplanmeierlib)
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
        '--outcome', dest='outcome', default=[], action='append',
        help='Where to find the outcome information in the clinical '
        'annotation file.  To analyze more than one outcome, use this '
        'parameter multiple times.  Format: <time_header>,<dead_header>')
    group.add_argument(
        '--gene', dest='gene', default=[], action='append',
        help='Name or ID of gene to analyze.  I will search for this gene '
        'in the annotations of the expression_file.  '
        'To analyze more than one gene, use this parameter multiple times.')
    group.add_argument(
        '--geneset', dest='geneset', default=[], action='append',
        help='Name of the geneset to analyze. To specify multiple gene sets, '
        'use this parameter multiple times.')
    group.add_argument(
        '--cutoff', dest='cutoff', default=None,
        help='Comma-separated list of breakpoints (between 0 and 1), '
        'e.g. 0.25,0.50')

    group = parser.add_argument_group(title='Ouput')
    group.add_argument(
        '-o', dest='filestem', default=None,
        help='Prefix used to name files.  e.g. "myanalysis".')
    group.add_argument(
        '--write_prism', dest='write_prism', action='store_true',
        default=False,
        help='Write a text file that can be imported into GraphPad Prism.')
    group.add_argument(
        '--plot_km', dest='plot_km', action='store_true', default=False,
        help='Write PNG-formatted Kaplan-Meier plot.')
    group.add_argument(
        '--xlab', dest='xlab', default=None, 
        help='the x label for Kaplan-Meier plot')
    group.add_argument(
        '--ylab', dest='ylab', default=None, 
        help='the y label for Kaplan-Meier plot')
    group.add_argument(
        '--title', dest='title', default=None,
        help='the title for Kaplan-Meier plot')

    args = parser.parse_args()

    # Check inputs.
    assert args.expression_file, (
        'Please specify a gene expression or gene set score file.')
    assert args.outcome_file, (
        'Please specify a clinical outcomes file.')
    assert os.path.exists(args.expression_file), "File not found: %s" % \
           args.expression_file
    assert os.path.exists(args.outcome_file), "File not found: %s" % \
           args.outcome_file
    assert args.outcome, 'Please specify the clinical outcomes to analyze.'

    assert args.gene or args.geneset, 'Please specify a gene or gene set.'
    assert not (args.gene and args.geneset), (
        'Please specify either a gene or a gene set, not both.')

    # Clean up the input.
    genes = parse_genes(args.gene)
    gene_sets = parse_gene_sets(args.geneset)
    cutoffs = parse_cutoffs(args.cutoff)
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
    gene_outcome_scores = {}  
    for x in itertools.product(outcomes, range(M.nrow())):
        (time_header, dead_header), i = x
        survival = clinical_annots[time_header]
        dead = clinical_annots[dead_header]
        scores = M.value(i, None)
        x = calc_association(survival, dead, scores, cutoffs)
        gene_outcome_scores[(time_header, dead_header, i)] = x


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

        # Groups
        expression_or_score = "Expression"
        if gene_sets:
            expression_or_score = "Score"
        groups = make_group_names(cutoffs, expression_or_score)

        # Relationship
        relationship = ""
        if SURV["hi_score_short_surv"]:
            relationship = "High %s has shorter outcome." % \
                           expression_or_score.lower()
        elif SURV["hi_score_short_surv"] is not None:
            relationship = "Low %s has shorter outcome." % \
                           expression_or_score.lower()

        gene_names = [M.row_names(x)[gene_i] for x in M.row_names()]
        outcome = time_header
        group_keys = [str(i) for i in range(len(cutoffs)+1)]
        num_samples = [SURV["num_samples"][x] for x in group_keys]
        mean_score = [SURV["mean_score"][x] for x in group_keys]
        surv90 = [SURV["surv90"][x] for x in group_keys]
        surv50 = [SURV["surv50"][x] for x in group_keys]
        p_value = SURV["p_value"]

        _fmt = _format_list
        x = gene_names + [
            outcome, _fmt(groups), _fmt(num_samples), _fmt(mean_score),
            _fmt(surv90), _fmt(surv50), relationship, p_value]
        assert len(x) == len(header)
        print >>outhandle, "\t".join(map(str, x))


        # Write out Prism, Kaplan-Meier curves, etc.
        gene_id = M.row_names(M.row_names()[0])[gene_i]
        gene_id_h = hashlib.hash_var(gene_id)

        if args.write_prism:
            filename = "%s%s.%s.prism.txt" % (filestem, time_header, gene_id_h)
            write_prism_file(
                filename, SURV["survival"], SURV["dead"], SURV["group"],
                groups)
        if args.plot_km:
            filename = "%s%s.%s.km.png" % (filestem, time_header, gene_id_h)
            plot_km(
                filename, SURV["survival"], SURV["dead"], SURV["group"],
                SURV["p_value"], gene_id, groups, 
                args.xlab, args.ylab, args.title)

            
if __name__ == '__main__':
    main()
