from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import mplgraph
        import arrayio
        from Betsy import module_utils
        from genomicode import jmath
        from genomicode import arrayplatformlib
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
        id_ = platforms[0][0]
        platform = platforms[0][1]
        if platform:
            if platform in ['HumanHT_12', 'MouseRef_8', 'HumanHT_12_control',
                            'MouseRef_8_control', 'entrez_ID_human',
                            'entrez_ID_mouse', 'entrez_symbol_human',
                            'entrez_symbol_mouse']:
                import matplotlib.pyplot as plt
                plt.clf()
                plt.plot([0, 0, 0, 0])
                plt.title('no AFFX plot can be generated')
                plt.savefig(outfile)

            else:
                M = arrayio.read(in_data.identifier)
                label = M._col_names['_SAMPLE_NAME']
                row_names = M._row_names[id_]
                index = []
                for i, name in enumerate(row_names):
                    if name.startswith('AFFX-'):
                        index.append(i)
                M_new = M.matrix(index)
                new = M_new.slice()
                a = jmath.mean_matrix(new, byrow=None)
                line = [(i, a[i]) for i in range(len(a))]
                f = mplgraph.lineplot(line,
                                      ylim_min=0,
                                      ylabel='Gene Expression Value',
                                      box_label=label)
                f.savefig(outfile)
            assert module_utils.exists_nz(outfile), (
                'the output file %s for plot_affy_affx_line fails' % outfile
            )




    
    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'control_plot_' + original_file + '.png'
        return filename


