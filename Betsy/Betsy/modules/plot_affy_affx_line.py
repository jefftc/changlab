#plot_affy_affx_line.py
import os
from genomicode import mplgraph, arrayplatformlib, jmath
import arrayio
from Betsy import module_utils, bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    M = arrayio.read(in_data.identifier)
    platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
    id = platforms[0][0]
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
            row_names = M._row_names[id]
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
        out_node = bie3.Data(rulebase.ControlPlot, **out_attributes)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object



def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'control_plot_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)

    return data_node
