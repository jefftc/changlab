from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from matplotlib import cm
        from Betsy import read_label_file
        
        data_node, cls_node = antecedents
        a, b, c = read_label_file.read(cls_node.identifier)
        if len(a) > 1:
            colors = []
            for i in range(5):
                colors.append(cm.hot(i / 5.0, 1))
                colors.append(cm.autumn(i / 5.0, i))
                colors.append(cm.cool(i / 5.0, i))
                colors.append(cm.jet(i / 5.0, i))
                colors.append(cm.spring(i / 5.0, i))
                colors.append(cm.prism(i / 5.0, i))
                colors.append(cm.summer(i / 5.0, i))
                colors.append(cm.winter(i / 5.0, i))
            opts = [colors[int(i)] for i in b]
            legend = [c[int(i)] for i in b]
            plot_pca(data_node.identifier, outfile, opts, legend)
        else:
            plot_pca(data_node.identifier, outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        # TODO: BUG? process or preprocess?  What does this line do?
        #data_node.attributes['process']
        filename = "Pca_%s_%s.png" % (
            original_file, data_node.attributes['preprocess'])
        return filename


def plot_pca(filename, result_fig, opts='b', legend=None):
    import arrayio
    from genomicode import jmath, mplgraph
    from genomicode import filelib

    R = jmath.start_R()
    jmath.R_equals(filename, 'filename')
    M = arrayio.read(filename)
    labels = M._col_names['_SAMPLE_NAME']
    data = M.slice()
    jmath.R_equals(data, 'X')
    R('NUM.COMPONENTS <- 2')
    R('S <- svd(X)')
    R('U <- S$u[,1:NUM.COMPONENTS]')
    R('D <- S$d[1:NUM.COMPONENTS]')
    # Project the data onto the first 2 components.
    R('x <- t(X) %*% U %*% diag(D)')
    x1 = R['x'][0:M.ncol()]
    x2 = R['x'][M.ncol():]

    xlabel = 'Principal Component 1'
    ylabel = 'Principal Component 2'
    
    if len(opts) > 1:
        fig = mplgraph.scatter(
            x1, x2, xlabel=xlabel, ylabel=ylabel, color=opts,
            legend=legend)
    else:
        fig = mplgraph.scatter(
            x1, x2, xlabel=xlabel, ylabel=ylabel, color=opts,
            label=labels)
    fig.savefig(result_fig)
    assert filelib.exists_nz(result_fig), 'the plot_pca.py fails'
