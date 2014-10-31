#extract_illumina_idat_files.py

import shutil
import os
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(data_node, parameters, user_input,network):
    outfile = name_outfile(data_node,user_input)
    directory = module_utils.unzip_if_zip(data_node.identifier)
    illumina_file = []
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in [f for f in filenames if f.endswith(".idat")]:
            illumina_file.append(os.path.join(dirpath, filename))
    if illumina_file:
        os.mkdir(outfile)
        for filename in illumina_file:
            newfilename = os.path.split(filename)[-1]
            if newfilename[:-5].endswith('_Grn'):
                newfilename = newfilename[:-9] + newfilename[-5:]
            new_file = os.path.join(outfile, newfilename)
            shutil.copyfile(filename, new_file)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_illumina_idat_files fails'
            % outfile)
        out_node = bie3.Data(rulebase.IDATFiles,**parameters)
    	out_object = module_utils.DataObject(out_node,outfile)
    	return out_object
    else:
        print 'There is no illumina idat file in the input.'
        return None

        
##def run(data_node, parameters, user_input,network):
##    outfile = name_outfile(data_node,user_input)
##    directory = module_utils.unzip_if_zip(data_node.identifier)
##    illumina_file = []
##    filenames = os.listdir(directory)
##    assert filenames, 'The input folder or zip file is empty.'
##    for filename in filenames:
##        if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
##            continue
##        if filename.endswith('.idat'):
##            illumina_file.append(filename)
##    if illumina_file:
##        os.mkdir(outfile)
##        for filename in illumina_file:
##            if filename[:-5].endswith('_Grn'):
##                newfilename = filename[:-9] + filename[-5:]
##            else:
##                newfilename = filename
##            old_file = os.path.join(directory, filename)
##            new_file = os.path.join(outfile, newfilename)
##            shutil.copyfile(old_file, new_file)
##        assert module_utils.exists_nz(outfile), (
##            'the output file %s for extract_illumina_idat_files fails'
##            % outfile)
##        out_node = bie3.Data(rulebase.IDATFiles,**parameters)
##    	out_object = module_utils.DataObject(out_node,outfile)
##    	return out_object
##    else:
##        print 'There is no illumina idat file in the input.'
##        return None
        


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(
        data_node.identifier)
    filename = 'idat_files_' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters,data_node):
    return parameters


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node
