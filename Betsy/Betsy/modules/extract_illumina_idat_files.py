#extract_illumina_idat_files.py
from Betsy import module_utils
import shutil
import os


def run(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    directory = module_utils.unzip_if_zip(single_object.identifier)
    illumina_file = []
    filenames = os.listdir(directory)
    assert filenames, 'The input folder or zip file is empty.'
    for filename in filenames:
        if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
            continue
        if filename.endswith('.idat'):
            illumina_file.append(filename)
    if illumina_file:
        os.mkdir(outfile)
        for filename in illumina_file:
            if filename[:-5].endswith('_Grn'):
                newfilename = filename[:-9] + filename[-5:]
            else:
                newfilename = filename
            old_file = os.path.join(directory, filename)
            new_file = os.path.join(outfile, newfilename)
            shutil.copyfile(old_file, new_file)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_illumina_idat_files fails'
            % outfile)
        new_objects = get_newobjects(parameters, objects, pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters, single_object, pipeline, outfile)
        return new_objects
    else:
        print 'There is no illumina idat file in the input.'
        return None


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'idat_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'idat_files', parameters, objects, single_object)
    return new_objects


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'idat_files', 'contents')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for extract_illumina_idat_files does not exist'
        % single_object.identifier)
    return single_object
