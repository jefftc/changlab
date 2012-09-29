#plot_control.py
import os
from Betsy import module_utils
import shutil
from genomicode import mplgraph,arrayannot,jmath
import arrayio
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    platforms = arrayannot.identify_all_platforms_of_matrix(M)
    id = platforms[0][0]
    platform = platforms[0][1]
    if platform:
        if platform in ['HumanHT_12','MouseRef_8',
                        'HumanHT_12_control','MouseRef_8_control',
                        'entrez_ID_human','entrez_ID_mouse',
                        'entrez_ID_symbol_human',
                        'entrez_ID_symbol_mouse']:
            return None
        else:
            M=arrayio.read(single_object.identifier)
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
            'the output file %s for plot_control fails'%outfile)
            new_objects = get_newobjects(parameters,objects,pipeline)
            module_utils.write_Betsy_parameters_file(
                      parameters,single_object,pipeline,outfile)
            return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'control_plot_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_control does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'control_plot',parameters,objects,single_object)
    return new_objects
