#extract_agilent_files.py

import shutil
import os
from Betsy import hash_method
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    directory = module_utils.unzip_if_zip(in_data.identifier)
    agilent_files = []
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    for filename in filenames:
        if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
            continue
        if os.path.isdir(os.path.join(directory, filename)):
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
        signal_tag = set(['gProcessedSignal', 'rProcessedSignal'])
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
            'the output file %s for extract_agilent_files fails' % outfile
        )
        out_node = bie3.Data(rulebase.AgilentFiles, **out_attributes)
        out_object = module_utils.DataObject(out_node, outfile)
        return out_object
    else:
        assert ValueError('There is no agilent file in the input.')



def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'agilent_files' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node
