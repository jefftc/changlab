#convert_v3_4.py

import os
import module_utils
import rule_engine
import shutil
import gzip

def run(parameters,objects):
    """convert from cc to v3_4"""
    from genomicode import affyio
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    filenames = os.listdir(identifier)
    os.mkdir(outfile)
    for filename in filenames:
        fileloc = os.path.join(identifier,filename)
        cel_v = affyio.guess_cel_version(fileloc)
        if fileloc.endswith('.gz'):
            newcelfname =  os.path.splitext(filename)[0]
            #unzip the gz data
            cel_file = fileloc[:-3]
            fileObj = gzip.GzipFile(fileloc, 'rb');
            fileObjOut = file(cel_file, 'wb');
            while 1:
                line = fileObj.readline()
                if line == '':
                    break
                fileObjOut.write(line)
            fileObj.close()
            fileObjOut.close()
            assert os.path.exists(cel_file),'the unzip fails'
        else:
                cel_file = fileloc
                newcelfname = filename
        if cel_v == 'cc1':
            f = file(os.path.join(outfile,newcelfname),'w')
            affyio.convert_cel_cc1_to_3(cel_file,f)
            f.close()
        elif cel_v in ('v3','v4'):
            shutil.copyfile(cel_file,os.path.join(outfile,newcelfname))
        if fileloc.endswith('.gz'):
            os.remove(cel_file)
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects

    
def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'geo_dataset','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'geo_dataset','Contents,DatasetId','geo_dataset')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
    
