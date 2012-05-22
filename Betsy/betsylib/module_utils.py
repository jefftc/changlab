#module_utils.py
import hash_method
import arrayio
from genomicode import binreg,Matrix,jmath,matrixlib
import rule_engine
import os
import read_label_file
import Betsy_config
import rule_engine

def get_inputid(identifier):
    old_filename = os.path.split(identifier)[-1]
    if '_BETSYHASH1_' in old_filename: 
        inputid = '_'.join(old_filename.split('_')[:-2])
    else:
        inputid = old_filename
    return inputid

"""contain some functions that are called by many modules"""
def make_unique_hash(identifier,pipeline,parameters):
    parameters = renew_parameters(parameters,['status'])
    inputid = get_inputid(identifier)
    new_parameters = parameters.copy()
    hash_result = hash_method.hash_parameters(
        inputid,pipeline,**new_parameters)
    return hash_result

def get_outfile(parameters,objects,in_objecttype,in_attribute,pipeline):
    identifier,single_object = find_object(
        parameters,objects,in_objecttype,in_attribute)
    original_file = get_inputid(identifier)
    hash_string = make_unique_hash(identifier,pipeline,parameters)
    filename = original_file + '_BETSYHASH1_' + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_newobjects(outfile,out_objecttype,parameters,objects,single_object):
    parameters = renew_parameters(parameters,['status'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject(out_objecttype,attributes,outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return new_objects

    
def merge_two_files(A_file,B_file,handle):
    """input two files and merge,write the output to handle"""
    M_A = arrayio.read(A_file)
    M_B = arrayio.read(B_file)
    assert arrayio.tab_delimited_format.is_matrix(M_A)
    assert arrayio.tab_delimited_format.is_matrix(M_B)
    [M_A,M_B] = matrixlib.align_rows(M_A,M_B)
    assert M_A.nrow() > 0, 'there is no common genes between two files'
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
        newsample_list = []
        for sample in M_B._col_names[name]:
            if sample in M_A._col_names[name]:
                newsample = sample + '_2'
            else:
                newsample = sample
            newsample_list.append(newsample)
        #x = M_A._col_names[name] + M_B._col_names[name]
        x = M_A._col_names[name] + newsample_list 
        col_names[name] = x
    M_c = Matrix.InMemoryMatrix(X,row_names,col_names,row_order)
    #M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
    arrayio.pcl_format.write(M_c,handle)
    
def which(program):
    
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    
    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate
    return None
        
def format_convert(X):
    data = []
    for i in range(X.dim()[1]):
        data.append(X.value(None,i))
    return data

def write_Betsy_parameters_file(parameters,single_object,pipeline):
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
    f.write('Pipeline module sequence:\n')
    f.write('\t'.join(pipeline))
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

def exists_nz(filename):
    if os.path.exists(filename):
        size = os.path.getsize(filename)
        if size > 0:
            return True
        else:
            return False
    else:
        return False
    
def plot_R(filename,keywords,outfile):
    from genomicode import jmath
    R = jmath.start_R()
    R('library(R.utils)')
    command = 'png2("' + outfile + '")'
    R(command)
    row = len(keywords)
    jmath.R_equals(row,'row')
    R('par(mfrow=c(row,1))')
    M = arrayio.read(filename)
    header = M.row_names()
    for keyword in keywords:
        data= []
        legend_name = []
        for i in range(M.dim()[0]):
            if M.row_names(header[1])[i] == keyword:
                data.append(M.slice()[i])
                legend_name.append('"'+M.row_names(header[0])[i]+'"')
        assert len(data)>0,'input is not a control file'
        max_value = max(data[0])
        for i in range(1,len(data)):
            max_value = max(max_value,max(data[i]))
        R('opts = c("red","blue","yellow","green","black","orangered","indianred",\
          "deepskyblue","pink")')
        jmath.R_equals(max_value,'max_value')
        jmath.R_equals_vector(data[0],'data')
        jmath.R_equals('"'+keyword+'"','keyword')
        R('plot(data,ylim=c(0,max_value),\
          xlab = "sample",ylab = keyword,type="l",col=opts[1])')
        for i in range(1,len(data)):
            jmath.R_equals(i,'i')
            jmath.R_equals_vector(data[i],'data')
            R('lines(data,col=opts[i+1])')
        jmath.R_equals_vector(legend_name,'legend_name')
        R('legend("bottomleft", legend_name,col=opts, lty=1,cex=0.8)')
    R('dev.off()')
    assert exists_nz(outfile),'the plot_R fails'
    

def renew_parameters(parameters,key_list):
    newparameters = parameters.copy()
    for key in key_list:
        if key in newparameters.keys():
            del newparameters[key]
    return newparameters
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def check_rename_file(filename):
    f=file(filename,'rU')
    text=f.read()
    f.close()
    lines=text.split('\n')
    lines=[line for line in lines if len(line)>0]
    for line in lines:
        assert len(line.split('\t'))==2,'the format of %s is not correct'%filename
    
def read_rename_file(filename):
    file_dict = dict()
    f=file(filename,'rU')
    text=f.read()
    f.close()
    lines=text.split('\n')
    lines=[line for line in lines if len(line)>0]
    for line in lines:
        line = line.split('\t')
        file_dict[line[0]]=line[1]
    return file_dict