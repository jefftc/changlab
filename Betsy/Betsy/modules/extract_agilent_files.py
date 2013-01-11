#extract_agilent_files.py
from Betsy import module_utils
import shutil
import os
from Betsy import hash_method


def run(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    directory = module_utils.unzip_if_zip(single_object.identifier)
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
        new_objects = get_newobjects(parameters, objects, pipeline)
        module_utils.write_Betsy_parameters_file(parameters, single_object,
                                                 pipeline, outfile)
        return new_objects
    else:
        print 'There is no agilent file in the input.'
        return None


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'agilent_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'agilent_files', parameters, objects, single_object)
    return new_objects


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'agilent_files', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for extract_agilent_files does not exist'
        % single_object.identifier)
    return single_object
