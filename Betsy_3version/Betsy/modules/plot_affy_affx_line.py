#plot_affy_affx_line.py
import os
import shutil
from genomicode import mplgraph,arrayplatformlib,jmath
import arrayio 
from Betsy import module_utils,bie3,rulebase


def run(data_node,parameters,user_input, network):
    outfile = name_outfile(data_node,user_input)
    M = arrayio.read(data_node.identifier)
    platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
    id = platforms[0][0]
    platform = platforms[0][1]
    if platform:
        if platform in ['HumanHT_12','MouseRef_8',
                        'HumanHT_12_control','MouseRef_8_control',
                        'entrez_ID_human','entrez_ID_mouse',
                        'entrez_symbol_human',
                        'entrez_symbol_mouse']:
            import matplotlib.pyplot as plt
            plt.clf()
            plt.plot([0,0,0,0])
            plt.title('no AFFX plot can be generated')
            plt.savefig(outfile)
            
        else:
            M=arrayio.read(data_node.identifier)
            label = M._col_names['_SAMPLE_NAME']
            row_names=M._row_names[id]
            index=[]
            for i,name in enumerate(row_names):
                if name.startswith('AFFX-'):
                    index.append(i)
            M_new=M.matrix(index)
            new = M_new.slice()
            a=jmath.mean_matrix(new,byrow=None)
            line=[(i,a[i]) for i in range(len(a))]
            f=mplgraph.lineplot(line,ylim_min=0,
                                ylabel='Gene Expression Value',box_label=label)
            f.savefig(outfile)
        assert module_utils.exists_nz(outfile),(
            'the output file %s for plot_affy_affx_line fails'%outfile)
   	out_node = bie3.Data(rulebase.ControlPlot,**parameters)
        out_object = module_utils.DataObject(out_node,outfile)
        return out_object



def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'control_plot_' + original_file + '.png'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    
    return data_node
