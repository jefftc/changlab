#plot_intensity.py
import os
import module_utils
import shutil

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    plot_intensity(identifier,outfile)
    assert module_utils.exists_nz(outfile),'the output\
                        file %s for plot intensity fails'%outfile
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
    assert os.path.exists(identifier),'the input\
                file %s for plot_intensity does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'intensity_plot',parameters,objects,single_object)
    return new_objects

def plot_intensity(filename,outfile):
    import arrayio
    from genomicode import jmath
    R=jmath.start_R()
    M = arrayio.read(filename)
    jmath.R_equals_matrix(M.slice(),'ma')
    R('library(R.utils)')
    command = 'png2("' + outfile + '")'
    R(command)
    R('boxplot(ma, main="signal intensity",\
       xlab="sample", ylab="signal")')
    R('dev.off()')
    assert module_utils.exists_nz(outfile),'plot_intensity fails'
