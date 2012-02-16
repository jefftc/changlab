#convert_res_gct.py
import os
import module_utils

def run(parameters,objects):
    """convert res signal file to gct format"""
    import arrayio
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    f=file(identifier,'r')
    text=f.readlines()
    f.close()
    first_line=text[0].split('\t')
    new_first_line=[i for i in first_line if len(i)>0]
    newtext=[]
    for single_line in text[3:]:
        single_line=single_line.split('\t')
        new_single_line=[i for i in single_line if i not in ['A','P','M']]
        newtext.append(new_single_line)
    f=file('tmp.txt','w')
    f.write('\t'.join(new_first_line))
    for i in newtext:
        f.write('\t'.join(i))
    f.close()
    os.remove('tmp.txt')
    M=arrayio.read('tmp.txt') 
    f=file(outfile,'w')
    M_c = arrayio.convert(M,to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c,f)
    f.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
