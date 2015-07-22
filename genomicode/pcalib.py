"""

Functions:
svd_project_cols
choose_colors
plot_scatter

select_genes_mv
select_genes_var

"""

import os, sys

def svd_project_cols(X, K):
    # Return (ncol x K matrix, array of percent variance).
    import numpy

    assert len(X) and len(X[0])
    nrow, ncol = len(X), len(X[0])
    nmin = min(nrow, ncol)
    
    # U  nrow x k  columns are principal components
    # V  k x ncol  rows are principal components
    U, s, V = numpy.linalg.svd(X, full_matrices=False)
    S = numpy.zeros((nrow, ncol))
    S[:nmin, :nmin] = numpy.diag(s)
    # Y  k x ncol
    X_hat = numpy.dot(numpy.transpose(U[:,:K]), X)
    X_hat = numpy.transpose(X_hat).tolist()
    # percent variance of each vector.
    x = [float(x*x) for x in s]
    perc_var = [y/sum(x) for y in x]
    return X_hat, perc_var

def choose_colors(group):
    # group should be a list of 0-based integers.  Can be None if no
    # group is assigned to a sample.
    import colorlib

    group_clean = [x for x in group if x is not None]
    if not group_clean:
        # No clusters specified.
        return None
    # If there is only 1 group specified, then color the samples in
    # that group, leaving the others black.
    if max(group_clean) == 0:
        pass
    # Use the middle range of the colorbar, because the colors at the
    # end are dark.
    # Doesn't work.  Middle colors are too pastel.
    #color_fn = colorlib.matlab_colors
    color_fn = colorlib.bild_colors
    palette = color_fn(max(group_clean)+1)
    color = [None] * len(group)
    for i in range(len(group)):
        if group[i] is None:
            continue
        color[i] = palette[group[i]]
    return color

def plot_scatter(
    X, Y, out_file, group=None, color=None, height=None, width=None,
    label=None, xlabel=True, ylabel=True,
    scale_label=None, title=None, pov_file=None, povray=None):
    # group should be a list of 0-N indicating the groupings of the
    # data points.  It should be the same length of X and Y.
    # Returns the output from povray.
    import tempfile
    import graphlib
    import filelib

    if width is None or height is None:
        width, height = 1024, 768
    assert width >= 16 and width < 4096*16
    assert height >= 16 and height < 4096*16
    if scale_label is None:
        scale_label = 1.0

    if not len(X):
        return None
    assert len(X) == len(Y)
    if not group:
        group = [0]*len(X)
    group_clean = [x for x in group if x is not None]
    if not group_clean:
        # All group are None.
        group = group_clean = [0]*len(X)
    assert len(group) == len(Y)
    assert min(group_clean) >= 0 and max(group_clean) < len(X)
    
    # If there is only 1 dataset, then just make everything black.
    if max(group_clean) > 0 and not color:
        color = choose_colors(group)

    is_tempfile = False
    try:
        if pov_file is None:
            is_tempfile = True
            x, pov_file = tempfile.mkstemp(suffix=".pov", dir="."); os.close(x)

        plot_width, plot_height = width, height
        points = zip(X, Y)
        graph = graphlib.scatter(
            points, color=color,
            xtick=True, xtick_label=xlabel, ytick=True, ytick_label=ylabel,
            overpoint_label=label, overpoint_label_size=1.5*scale_label,
            xlabel="Principal Component 1", ylabel="Principal Component 2",
            label_size=1, 
            title=title,
            width=plot_width, height=plot_height)
        output = graph.write(out_file, povray_bin=povray)
        #open(pov_file, 'w').write(graph.draw())
        ## povray -D -J +Opredictions.png -H768 -W1024 +A0.5 predictions.pov
        #r = povraygraph.povray(
        #    pov_file, outfile=out_file,
        #    height=plot_height, width=plot_width, antialias=0.5, quality=9,
        #    povray_bin=povray)
        #output = r.read()
    finally:
        if is_tempfile and pov_file and os.path.exists(pov_file):
            os.unlink(pov_file)
    assert filelib.exists_nz(out_file), "Failed to plot predictions.\n%s" % (
        output)
    return output

def select_genes_mv(X, num_genes_mean=None, num_genes_var=None,
                    test_for_missing_values=False):
    # Return the indexes of the genes that pass both the mean and
    # variance threshold.  num_genes_mean and num_genes_var is the
    # number of genes to keep.
    import jmath

    if num_genes_mean:
        assert type(num_genes_mean) is type(0)
        assert num_genes_mean <= len(X), "%d %d" % (num_genes_mean, len(X))
    if num_genes_var:
        assert type(num_genes_var) is type(0)
        assert num_genes_var <= len(X), "%d %d" % (num_genes_var, len(X))
    
    # BUG: If many genes have the same mean or variance, then may
    # return more than num_genes_mean or num_genes_var.
    #cutoff_mean = min(means)
    #cutoff_var = min(vars)
    #if num_genes_mean:
    #    cutoff_mean = sorted(means)[-num_genes_mean]
    #if num_genes_var:
    #    cutoff_var = sorted(vars)[-num_genes_var]
    #I = [i for i in range(len(vars))
    #     if means[i] >= cutoff_mean and vars[i] >= cutoff_var]

    I_mean = range(len(X))
    I_var = range(len(X))
    if num_genes_mean:
        if not test_for_missing_values:
            means = jmath.mean(X)
        else:
            means = []
            for x in X:
                x = [x for x in x if x is not None]
                m = None  # If all missing values, then mean is None.
                if x:
                    m = jmath.mean(x)
                means.append(m)
        I_mean = jmath.order(means, decreasing=1)[:num_genes_mean]
    if num_genes_var:
        if test_for_missing_values:
            raise NotImplementedError
        vars = jmath.var(X)
        I_var = jmath.order(vars, decreasing=1)[:num_genes_var]
    I_mean = {}.fromkeys(I_mean)
    I_var = {}.fromkeys(I_var)
    # Bug: I_mean and I_var can contain None.  Should handle this case.
    I = [i for i in range(len(X)) if i in I_mean and i in I_var]
    #for i in range(len(X)):
    #    x = i, vars[i], int(i in I)
    #    print "\t".join(map(str, x))
    #import sys; sys.exit(0)
    return I

def select_genes_var(X, num_genes):
    # Return the indexes of the genes with the greatest variance.
    return select_genes_mv(X, num_genes_var=num_genes)
