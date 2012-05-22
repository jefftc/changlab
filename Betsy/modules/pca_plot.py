#pca_plot.py

import os
import module_utils
import shutil
import read_label_file

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file,obj = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    if label_file:
        a,b,c=read_label_file.read(label_file)
        colors = ['"red"','"blue"','"green"','"yellow"']
        opts = [colors[int(i)] for i in b]
        plot_pca(identifier,outfile,opts)
    else:
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
        parameters,objects,'signal_file','contents,preprocess',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(identifier),'the input file %s\
                    for pca_plot does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'pca_plot',parameters,objects,single_object)
    return new_objects

def plot_pca(filename,result_fig,opts='"red"'):
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
    jmath.R_equals(opts,'opts')
    R(command)
    R('plot(x, xlab="", ylab="",col=opts)')
    R('dev.off()')
    assert module_utils.exists_nz(result_fig),'the plot_pca.py fails'

def find_object(parameters,objects,objecttype,attribute):
    identifier = None
    single_object = None
    attributes = attribute.split(',')
    for i in range(len(attributes)):  #consider the content is [unknown]
            if '[' in attributes[i]:
                attribute = attributes
            else:
                attribute = parameters[attributes[i]]
    compare_attribute = [parameters[i] for i in attributes]
    
    for single_object in objects:
        flag = True
        if objecttype in single_object.objecttype:
            for i in compare_attribute:
                if i not in single_object.attributes:
                    flag=False
            if flag:
                identifier=single_object.identifier
                break
    return identifier,single_object