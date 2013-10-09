#extract_agilent_files.py
#from Betsy
import module_utils
import shutil
import os
from Betsy import hash_method
import bie
import rulebase

def run(data_node, parameters):
    outfile = name_outfile(data_node)
    directory = module_utils.unzip_if_zip(data_node.attributes['filename'])
    agilent_files = []
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    for filename in filenames:
        if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
            continue
        if os.path.isdir(os.path.join(directory,filename)):
            continue
        postag = []
        fline = []
        f = open(os.path.join(directory, filename), 'r')
        for i in range(10):
            line = f.readline()
            words = line.split()
            if len(words) > 0:
                postag.append(words[0])
                if words[0] == 'FEATURES':
                    fline = set(words)
        f.close()
        signal_tag = set(['gProcessedSignal','rProcessedSignal'])
        if signal_tag.issubset(fline):
            if postag == ['TYPE', 'FEPARAMS', 'DATA', '*', 'TYPE', 'STATS',
                          'DATA', '*', 'TYPE', 'FEATURES']:
                agilent_files.append(filename)
                
    if agilent_files:
        if not os.path.exists(outfile):
            os.mkdir(outfile)
        for filename in agilent_files:
            old_file = os.path.join(directory, filename)
            new_file = os.path.join(outfile, filename)
            shutil.copyfile(old_file, new_file)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_agilent_files fails' % outfile)
        new_parameters = parameters.copy()
        new_parameters['filename'] = os.path.split(outfile)[-1]
        out_node = bie.Data(rulebase.AgilentFiles, **new_parameters)
        return out_node
    else:
        print 'There is no agilent file in the input.'
        return None

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def get_out_attributes(parameters,data_node):
    return parameters


def name_outfile(data_node):
    original_file = module_utils.get_inputid(data_node.attributes['filename'])
    filename = 'agilent_files' + original_file 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes)
    return data_node
