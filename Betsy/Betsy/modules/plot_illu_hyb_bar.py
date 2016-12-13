from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib

        in_data = antecedents
        plot_hyb_bar(in_data.identifier, outfile)
        assert filelib.exists_nz(outfile), (
            'the output file %s for plot_hyb_bar fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'hyb_bar_' + original_file + '.png'
        return filename


def plot_hyb_bar(filename, outfile):
    from genomicode import mplgraph
    from genomicode import filelib
    import math
    import numpy
    
    high = ['ILMN_2038770', 'ILMN_2038769']
    med = ['ILMN_2038768', 'ILMN_2038771']
    low = ['ILMN_1343050', 'ILMN_1343052']
    high_data = []
    med_data = []
    low_data = []
    import arrayio
    M = arrayio.read(filename)
    header = M.row_names()
    for i in range(M.dim()[0]):
        if not M.row_names(header[1])[i] == 'cy3_hyb':
            continue
        if M.row_names(header[0])[i] in high:
            high_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in med:
            med_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in low:
            low_data.extend(M.slice()[i])
    
    
    mean = [numpy.mean(high_data), numpy.mean(med_data), numpy.mean(low_data)]
    flag = [math.isnan(i) for i in mean]
    assert True not in flag, 'input is not a control file'
    std = [numpy.std(high_data), numpy.std(med_data), numpy.std(low_data)]
    fig = mplgraph.barplot(mean, std,
                           ylabel='Signal',
                           box_label=['high', 'med', 'low'])
    fig.savefig(outfile)
    assert filelib.exists_nz(outfile), 'the plot_illu_hyb_bar.py fails'
