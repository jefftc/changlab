#plot_intensity.py
import os
import module_utils
import shutil
import math
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    plot_intensity(single_object.identifier,outfile)
    assert module_utils.exists_nz(outfile),'the output\
                        file %s for plot intensity fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'intensity_'+original_file+'.png'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),'the input\
                file %s for plot_intensity does not exist'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'intensity_plot',parameters,objects,single_object)
    return new_objects

def plot_intensity(filename,outfile):
    import arrayio
    from genomicode import jmath
    R=jmath.start_R()
    M = arrayio.read(filename)
    jmath.R_equals_matrix(M.slice(),'ma',by_row=True)
    R('library(R.utils)')
    command = 'png2("' + outfile + '")'
    R(command)
    R('boxplot(ma, main="signal intensity",\
       xlab="sample", ylab="signal",axes=F)')
    label = ['""']*M.ncol()
    name = M._col_names.keys()[0]
    if M.ncol()<=12:
        label = ['"'+M._col_names[name][i]+'"' for i in range(M.ncol())]
    else:
        index = [int(round(M.ncol()/12.0*i)) for i in range(12)]
        for i in range(12):
            label[index[i]] = '"'+M._col_names[name][index[i]]+'"'
    data = jmath.transpose(M._X)
    max_value = max(data[0])
    for i in range(1,len(data)):
       max_value = max(max_value,max(data[i]))
    jmath.R_equals(M.ncol(),'ncol')
    jmath.R_equals(label,'label')
    jmath.R_equals(max_value,'max_value')
    R('axis(1,at=seq(1,ncol,by=1),lab=label,cex.axis=0.7,las=3)')
    R('axis(2,at=seq(0,ceiling(max_value)),cex.axis=0.8)')
    R('dev.off()')
    assert module_utils.exists_nz(outfile),'plot_intensity fails'
