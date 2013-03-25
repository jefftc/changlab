#normalize_sampels_with_shiftscale.py

import os
from Betsy import module_utils
import shutil
from Betsy import read_label_file
from genomicode import shiftscalenorm
import arrayio
from time import strftime,localtime


def run(parameters,objects,pipeline,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),(
        'cannot find label_file %s'%label_file.identifier)
    M = arrayio.read(single_object.identifier)
    result,label_line,second_line=read_label_file.read(label_file.identifier)
    assert len(result) == 2, 'for shiftscale,there should be only 2 classes'
    index1=result[0][0]
    index2=result[1][0]
    M_1=M.matrix(None,index1)
    M_2=M.matrix(None,index2)
    M_y = shiftscalenorm.normalize(M_1,M_2)
    for i in range(M_y.dim()[0]):
        for j in range(M_y.dim()[1]):
            if str(M_y._X[i][j]) == 'nan':         
                M_y._X[i][j] = M_2._X[i][0]
    for j in range(M.nrow()):
        for i in range(len(index1)):
            M._X[j][index1[i]]=M_y._X[j][i]
        
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(M,f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for shiftscale fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,
                                             [single_object,label_file],
                                             pipeline,outfile,starttime,user,jobname)
    return new_objects



def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_shiftscale_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
                    parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for shiftscale does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
