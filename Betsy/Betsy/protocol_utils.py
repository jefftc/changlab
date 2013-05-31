#protocol_utils.py
import os
import shutil
import config
import hash_method
import time
import imghdr
import module_utils

COMMON = 'Common Parameters'
NORMALIZE = 'Normalize Parameters'
OPTIONAL = 'Optional Parameters'
ILLUMINA = 'Illumina Normalize Parameters'
CLASS_NEIGHBORS='Class Neighbor Parameters'

class Parameter:
    def __init__(self,name,pretty_name,default_value,
                 type,choices,category,description):
        self.name = name
        self.pretty_name = pretty_name
        self.default = default_value
        self.type = type
        self.choices = choices
        self.category = category
        self.description = description

        
def get_result_folder(protocol, outfiles, parameters, pipeline, foldername):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0][0])
    folder_string = hash_method.hash_parameters(
        inputid, pipeline[0][0], **parameters[0][0])
    folder_name = foldername + '_BETSYHASH1_' + folder_string
    result_folder = os.path.join(OUTPUTPATH, folder_name)
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    final_outfiles = []
    for i in range(len(outfiles)):
        result_files = []
        new_outfiles = []
        for j in range(len(outfiles[i])):
                final_output = os.path.split(outfiles[i][j])[-1]
                output_folder = os.path.split(
                    os.path.split(outfiles[i][j])[-2])[-1]
                new_filename = output_folder + '_' + final_output
                new_outfiles.append(new_filename)
                result_file = os.path.join(result_folder, new_filename)
                if not os.path.exists(result_file):
                    if os.path.isdir(outfiles[i][j]):
                        shutil.copytree(outfiles[i][j], result_file)
                    else:
                        shutil.copyfile(outfiles[i][j], result_file)
        final_outfiles.append(new_outfiles)
    summarize_report(protocol, final_outfiles,
                     result_folder, parameters, pipeline)
    print 'Report:' + os.path.join(result_folder, 'report.html')


def format_prolog_query(predicate, parameters, modules):
    parameters = [str(i) for i in parameters]
    str_parameters = ','.join(parameters)
    output = str('[' + str_parameters + '],' + modules)
    query = predicate + '(' + output + ')'
    return query


def import_protocol(protocol):
    protocol_name = 'Betsy.protocols.' + protocol
    mod = __import__(protocol_name,  globals(), locals(),
                     [protocol], -2)
    return mod


def pretty_hostname():
    import subprocess
    cmd = "hostname"
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    wh, r = p.stdin, p.stdout
    wh.close()
    hostname = r.read().strip()
    assert hostname, "I could not get the hostname."
    return hostname


##def check_parameters(parameters):
##    allow_type = ['string', 'float', 'integer', 'list']
##    for value in parameters.itervalues():
##        if not isinstance(value, list):
##            assert value in allow_type, '%s is not a allow type' % value
##    return True

def check_parameters(PARAMETERS):
    allow_type = ['string', 'float', 'integer', 'list']
    for parameter in PARAMETERS:
        if parameter.type:
            assert parameter.type in allow_type, '%s is not a allow type' % parameter.type
    return True


##def check_default(default, parameters):
##    for key in default.keys():
##        assert key in parameters.keys(), '%s is not a valid parameter' % key
##        value = str(default[key])
##        if isinstance(parameters[key], list):
##            assert value in parameters[key], (
##                '%s is not correct value for %s' % (value, key))
##        elif parameters[key] == 'float':
##            assert module_utils.isnumber(value), (
##                '%s is not a float number for %s' % (value, key))
##        elif parameters[key] == 'string':
##            assert isinstance(value, str), '%s is not a string for %s' % (
##                str(value), key)
##        elif parameters[key] == 'integer':
##            assert value.isdigit(), '%s is not a digit for %s' % (
##                str(value), key)
##        elif parameters[key] == 'list':
##            assert isinstance(value, str), (
##                '%s is not a string which will convert to a list for %s' % (
##                    str(value), key))
##    return True


def check_default(PARAMETERS):
    for parameter in PARAMETERS:
        if parameter.default:
            if parameter.type:
                if parameter.type == 'float':
                    assert module_utils.is_number(parameter.default), (
                    '%s is not a float number for %s' % (parameter.default, parameter.name))
                elif parameter.type == 'string':
                    assert isinstance(parameter.default,str), (
                    '%s is not a string for %s' % (parameter.default, parameter.name))
                elif parameter.type == 'integer':
                    assert parameter.default.isdigit(), (
                    '%s is not a digit for %s' % (parameter.default, parameter.name))
                elif parameter.type == 'list':
                    assert isinstance(parameter.default,str), (
                    '%s is not a string which will convert to a list for %s'
                    % (parameter.default, parameter.name))
            else:
                assert parameter.default in parameter.choices,('%s is not correct value for %s'
                                                              % (parameter.default, parameter.name))
    return True

        
