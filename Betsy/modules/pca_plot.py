#pca_plot.py

import os
import module_utils
import shutil
def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    plot_pca(identifier,outfile)
    assert module_utils.exists_nz(outfile),'the output file %s\
                                for pca_plot fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input file %s\
                    for pca_plot does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'pca_plot',parameters,objects,single_object)
    return new_objects

def plot_pca(filename,result_fig):
    from genomicode import jmath
    import arrayio
    R=jmath.start_R()
    jmath.R_equals("\'"+filename+"\'",'filename')
    M = arrayio.read(filename)
    data = M.slice()
    jmath.R_equals(data,'X')
    R('NUM.COMPONENTS <- 2')
    R('S <- svd(X)')
    R('U <- S$u[,1:NUM.COMPONENTS]')
    R('D <- S$d[1:NUM.COMPONENTS]')
    # Project the data onto the first 2 components.
    R('x <- t(X) %*% U %*% diag(D)')
    R('x.all <- t(X) %*% S$u %*% diag(S$d)')
    pca_file = "\'" + result_fig + "\'"
    jmath.R_equals(pca_file,'pca_file')
    R('library(R.utils)')
    command='png2(pca_file)'
    R(command)
    R('plot(x, xlab="", ylab="")')
    R('dev.off()')
    assert module_utils.exists_nz(result_fig),'the plot_pca.py fails'
