#plot_pca.py

import os
import module_utils
import shutil
import read_label_file
from genomicode import pcalib
import arrayio

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    M = arrayio.read(single_object.identifier)
    X = M._X
    if 'pca_gene_num' in parameters.keys():
        N = int(parameters['pca_gene_num'])
    else:
        N = 500
    N = min(N,M.nrow())
    index = pcalib.select_genes_var(X,N)
    M_new = M.matrix(index,None)
    tmp = 'tmp'
    f = file(tmp,'w')
    arrayio.tab_delimited_format.write(M_new,f)
    f.close()
    a,b,c = read_label_file.read(label_file.identifier)
    if len(a)>1:
        colors = ['r','b','g','y']
        opts = [colors[int(i)] for i in b]
        legend = [c[int(i)] for i in b]
        module_utils.plot_pca(tmp,outfile,opts,legend)
    else:
        module_utils.plot_pca(tmp,outfile)
    os.remove(tmp)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for pca_plot fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = parameters['objecttype']+'_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for pca_plot does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'pca_plot',parameters,objects,single_object)
    return new_objects



