#extract_gpr_files.py
import os
import gzip
import shutil
from Betsy import module_utils
from Betsy import gpr_module


def run(parameters, objects, pipeline, options=None):
    """extract the files that are gpr format"""
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    check = []
    newfiles = []
    files = os.listdir(single_object.identifier)
    for i in files:
        if i == '.DS_Store':
            pass
        else:
            gpr = gpr_module.check_gpr(
                os.path.join(single_object.identifier, i))
            check.append(gpr)
            newfiles.append(i)
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    for i in range(len(check)):
        if check[i]:
            old_file = os.path.join(single_object.identifier, newfiles[i])
            new_file = os.path.join(outfile, newfiles[i])
            shutil.copyfile(old_file, new_file)
    if True in check:
        assert module_utils.exists_nz(outfile), (
            'the output file %s for extract_gpr fails' % outfile)
        new_objects = get_newobjects(parameters, objects, pipeline)
        module_utils.write_Betsy_parameters_file(
            parameters, single_object, pipeline, outfile)
        return new_objects
    else:
        return None


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'gpr_files_' + original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'gpr_files', parameters, objects, single_object)
    return new_objects


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, 'gpr_files', 'contents')
    assert os.path.exists(single_object.identifier), (
        'folder %s for are_all_gpr does not exit.' % single_object.identifier)
    assert os.path.isdir(single_object.identifier), (
        "input %s for are_all_gpr is not a folder" % single_object.identifier)
    return single_object
