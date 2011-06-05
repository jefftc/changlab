"""

Functions:
svd_project_cols
plot_scatter

select_genes_mv
select_genes_var

"""

import os, sys

def svd_project_cols(X, K):
    # Return ncol x K matrix.
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
    return X_hat

def choose_colors(group):
    import colorlib

    if not group or max(group) == 0:
        return None
    # Use the middle range of the colorbar, because the colors at the
    # end are dark.
    # Doesn't work.  Middle colors are too pastel.
    #color_fn = colorlib.matlab_colors
    color_fn = colorlib.bild_colors
    palette = color_fn(max(group)+1)
    color = [palette[x] for x in group]
    return color

def plot_scatter(X, Y, out_file, group=None, color=None,
                 pov_file=None, povray=None):
    # group should be a list of 0-N indicating the groupings of the
    # data points.  It should be the same length of X and Y.
    # Returns the output from povray.
    import povraygraph
    import tempfile

    if not len(X):
        return None
    assert len(X) == len(Y)
    if not group:
        group = [0]*len(X)
    assert len(group) == len(Y)
    assert min(group) >= 0 and max(group) < len(X)
    
    # If there is only 1 dataset, then just make everything black.
    if max(group) > 0 and not color:
        color = choose_colors(group)

    is_tempfile = False
    try:
        if pov_file is None:
            is_tempfile = True
            x, pov_file = tempfile.mkstemp(suffix=".pov", dir="."); os.close(x)
            
        plot_width, plot_height = 1024, 768
        x = povraygraph.scatter(
            X, Y, color=color,
            xtick=True, xtick_label=True, ytick=True, ytick_label=True,
            xlabel="Principal Component 1", ylabel="Principal Component 2",
            label_size=1, width=plot_width, height=plot_height)
        open(pov_file, 'w').write(x)
        # povray -D -J +Opredictions.png -H768 -W1024 +A0.5 predictions.pov
        r = povraygraph.povray(
            pov_file, outfile=out_file,
            height=plot_height, width=plot_width, antialias=0.5, quality=9,
            povray_bin=povray)
        output = r.read()
    finally:
        if is_tempfile and pov_file and os.path.exists(pov_file):
            os.unlink(pov_file)
    assert os.path.exists(out_file), "Failed to plot predictions."
    return output

def select_genes_mv(X, num_genes_mean=None, num_genes_var=None):
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
    
    means = jmath.mean(X)
    vars = jmath.var(X)

    cutoff_mean = min(means)
    cutoff_var = min(vars)
    if num_genes_mean:
        cutoff_mean = sorted(means)[-num_genes_mean]
    if num_genes_var:
        cutoff_var = sorted(vars)[-num_genes_var]
    I = [i for i in range(len(vars))
         if means[i] >= cutoff_mean and vars[i] >= cutoff_var]
    return I

def select_genes_var(X, num_genes):
    # Return the indexes of the genes with the greatest variance.
    return select_genes_mv(X, num_genes_var=num_genes)
