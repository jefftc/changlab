#module_utils.py
import hash_method
import arrayio
from genomicode import binreg,Matrix,jmath
import rule_engine
import os
import read_label_file
import Betsy_config

"""contain some functions that are called by many modules"""
def make_unique_hash(parameters,objects,objecttype,attribute):
    identifier,single_object = find_object(parameters,objects,objecttype,attribute)
    filename = os.path.split(identifier)[-1]
    if '_BETSYHASH_' in filename: 
        original_file = '_'.join(filename.split('_')[:-2])
    else:
        original_file = filename
    new_parameters = parameters.copy()
    new_parameters['filename'] = original_file
    hash_result = hash_method.hash_parameters(**new_parameters)
    return hash_result

def get_outfile(parameters,objects,in_objecttype,in_attribute,out_objecttype):
    identifier,single_object = find_object(parameters,objects,in_objecttype,in_attribute)
    old_filename = os.path.split(identifier)[-1]
    if '_BETSYHASH_' in old_filename: 
        original_file = '_'.join(old_filename.split('_')[:-2])
    else:
        original_file = old_filename
    hash_string = make_unique_hash(parameters,objects,in_objecttype,in_attribute)
    filename = original_file +'_BETSYHASH_'+ hash_string
    outfile = os.path.join(os.getcwd(),filename)
    if 'Status' in parameters.keys():
        del parameters['Status']
    attributes = parameters.values()
    new_object = rule_engine.DataObject(out_objecttype,attributes,outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return outfile,new_objects

    
def merge_two_files(A_file,B_file,handle):
    """input two files and merge,write the output to handle"""
    M_A = arrayio.read(A_file)
    M_B = arrayio.read(B_file)
    assert arrayio.tab_delimited_format.is_matrix(M_A)
    assert arrayio.tab_delimited_format.is_matrix(M_B)
    [M_A,M_B] = binreg.align_rows(M_A,M_B)
    X = []
    for i in range(M_A.dim()[0]):
        x = M_A._X[i]+M_B._X[i]
        X.append(x)
    row_names = M_A._row_names
    row_order = M_A._row_order
    col_names = {}
    for name in M_A._col_names:
        if name not in M_B._col_names:
            continue                                                                                                         
        x = M_A._col_names[name] + M_B._col_names[name]
        col_names[name] = x
    M_c = Matrix.InMemoryMatrix(X,row_names,col_names,row_order)
    #M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
    arrayio.pcl_format.write(M_c,handle)
    

        
def format_convert(X):
    data = []
    for i in range(X.dim()[1]):
        data.append(X.value(None,i))
    return data

def write_Betsy_parameters_file(parameters,single_object):
    f = file(os.path.join(os.getcwd(),'Betsy_parameters.txt'),'w')
    f.write('Module input:\n')
    if not isinstance(single_object,list):
        f.write(single_object.objecttype + '\n') 
        f.write(single_object.identifier + '\n')
    else:
        for one_object in single_object:
            f.write(one_object.objecttype + '\n') 
            f.write(one_object.identifier + '\n')
    f.write('Module output:\n')
    for i in parameters.keys():
         f.write(i+':'+parameters[i]+'\n')
    #f.write('Pipeline module sequence:\n')
    #for analysis in pipeline:
    #    f.write(analysis.name+'\n')
    f.close()
    

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
def run_gp_module(module_name,parameters):
    """given the module_name and the module parameters
       in dict, call module in Genepatttern"""
    R=jmath.start_R()
    username='\"'+Betsy_config.USERNAME+'\"'
    password='\"'+Betsy_config.PASSWORD+'\"'
    servername='\"'+Betsy_config.SERVERNAME+'\"'
    jmath.R_equals(password,'password')
    jmath.R_equals(servername,'servername')
    jmath.R_equals(username,'username')
    command="\'"+module_name+"\'"
    for key in parameters.keys():
        command=command+','+key+'='+'\"'+parameters[key]+'\"'
    fullcommand='result<-run.analysis(gp.client,'+command+')'
    R('library(GenePattern)')
    R('gp.client <- gp.login(servername, username, password)')
    R(fullcommand)
    R('download.directory <- job.result.get.job.number(result)')
    R('download.directory <- as.character(download.directory)')
    R('job.result.download.files(result, download.directory)')
    download_directory=os.path.realpath(R('download.directory')[0])
    result_files = os.listdir(download_directory)
    if 'stderr.txt' in result_files:
        return None
    else:
        return download_directory


