#check_v3_4.py
import os
import hash_method
import shutil
import rule_engine
import module_utils
def run(parameters,objects,pipeline):
    """check if all the cel file are v3_4"""
    from genomicode import affyio
    identifier,single_object = get_identifier(parameters,objects)
    assert os.path.exists(identifier),'folder %s does not exit.' % identifier
    assert os.path.isdir(identifier),"input is not a folder"
    outfile,new_objects = get_outfile(parameters,objects)
    cel_vdict = dict()
    for i in os.listdir(identifier):
        if i == '.DS_Store':
            pass
        else:
            cel_v = affyio.guess_cel_version(os.path.join(identifier,i))
            cel_vdict[cel_v] = cel_vdict.get(cel_v,0) + 1
    cel_vlist = cel_vdict.keys()
    if ((len(cel_vlist) == 1 and cel_vlist[0] in ('v3','v4')) or
        (len(cel_vlist) == 2 and 'v3' in  cel_vlist and 'v4' in cel_vlist)):
        os.mkdir(outfile)
        for i in os.listdir(identifier):
            if i != '.DS_Store':
                old_file = os.path.join(identifier,i)
                new_file = os.path.join(outfile,i)
                shutil.copyfile(old_file,new_file)
        module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
        return new_objects
    else:
         return None
        
    
def make_unique_hash(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    hash_profile={'filename':os.path.split(identifier)[-1],
                   'version': 'v3_4',
                   'number of files':str(len(os.listdir(identifier)))}
    hash_result=hash_method.hash_parameters(**hash_profile)
    return hash_result


def get_outfile(parameters,objects):
    identifier,single_object = get_identifier(parameters,objects)
    old_filename = os.path.split(identifier)[-1]
    if '_BETSYHASH1_' in old_filename: 
        original_file = '_'.join(old_filename.split('_')[:-2])
    else:
        original_file = old_filename
    hash_string = make_unique_hash(parameters,objects)
    filename = original_file + hash_string
    outfile = os.path.join(os.getcwd(),filename)
    attributes = parameters.values()
    objecttype = 'geo_dataset'
    new_object = rule_engine.DataObject(objecttype,attributes,outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return outfile,new_objects

def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'geo_dataset','Contents,DatasetId')
 
