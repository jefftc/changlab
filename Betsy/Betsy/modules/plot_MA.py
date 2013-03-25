#plot_MA.py
from Betsy import module_utils
import shutil
import os
from genomicode import jmath
from time import strftime,localtime

def run(parameters,objects,pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    cwd = os.getcwd()
    R = jmath.start_R()
    R('require(limma,quietly=TRUE)')
    R('library(marray)')
    os.chdir(single_object.identifier)
    try:
        R('dir<-getwd()')
        R('files<-list.files(dir)')
        R('x.read<-read.Agilent(files)')
    finally:
        os.chdir(cwd)
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    R('xnorm.loc <- maNorm(x.read, norm = "loess")')
    R('x.norm <- maNormScale(xnorm.loc, norm = "p")')
    jmath.R_equals(os.path.join(outfile,'before_boxplot.png'),'filename')
    R('bitmap(file=filename,type="png256")')
    R('boxplot(x.read, main="Before Normalization",xlab="samples", ylab="ratio")')
    R('dev.off()')
    jmath.R_equals(os.path.join(outfile,'after_boxplot.png'),'filename')
    R('bitmap(file=filename,type="png256")')
    R('boxplot(x.norm, main="After Normalization",xlab="samples", ylab="ratio")')
    R('dev.off()')
    jmath.R_equals(os.path.join(outfile,'before_plot.png'),'filename')
    R('bitmap(file=filename,type="png256")')
    R('plot(x.read, main="Before Normalization",xlab="samples", ylab="ratio")')
    R('dev.off()')
    jmath.R_equals(os.path.join(outfile,'after_plot.png'),'filename')
    R('bitmap(file=filename,type="png256")')
    R('plot(x.norm, main="After Normalization",xlab="samples", ylab="ratio")')
    R('dev.off()')
    assert module_utils.exists_nz(outfile),(
        'the output file %s for plot_MA' %outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
    parameters,single_object,pipeline,outfile,starttime,user,jobname)
    return new_objects
    

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'MA_plots_' + original_file
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
            parameters,objects,'agilent_files','contents',['agilent'])
    assert os.path.exists(single_object.identifier),(
        'the input file %s for plot_MA does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = module_utils.DataObject('ma_plot',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
