#plot_intensity.py
import os
import module_utils
import shutil

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    plot_intensity(identifier,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId','intensity_plot',pipeline)
    
def get_identifier(parameters,objects):
    return module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')

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
