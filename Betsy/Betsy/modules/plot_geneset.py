#plot_geneset.py

import os
from Betsy import module_utils
import shutil
from genomicode import mplgraph, filelib, jmath


def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    matrix = [x for x in filelib.read_cols(single_object.identifier)]
    matrix = [x[1:] for x in matrix]
    matrix = jmath.transpose(matrix)
    sample = matrix[0][1:]
    data = matrix[1:]
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for one_data in data:
        value = one_data[1:]
        value = [float(i) for i in value]
        pair = [(value[i],sample[i]) for i in range(len(value))]
        pair.sort()
        gene_value = [i[0] for i in pair]
        label = [i[1] for i in pair]
        ylabel=one_data[0]
        from genomicode import mplgraph
        fig=mplgraph.barplot(gene_value,box_label=label,xtick_rotation=90,
                             xlabel='sample',ylabel=ylabel)
        output = os.path.join(outfile,ylabel)
        fig.savefig(output+'.png')
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_geneset fails'%outfile)
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
    filename = 'geneset_plot_'+original_file
    outfile = os.path.join(os.getcwd(),filename)
    return outfile


def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'geneset_analysis','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_geneset does not exist'
        %single_object.identifier)
    return single_object


def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'geneset_plot',parameters,objects,single_object)
    return new_objects
