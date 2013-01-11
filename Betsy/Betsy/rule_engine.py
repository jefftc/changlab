#! /usr/bin/env python
#rule_engine.py
import os
import swipl
import shutil
import sys
import logging
import json
import config
import module_utils
import re
import time
import tempfile


def plstring2dataobject(pl_inputs, identifiers):
    assert len(pl_inputs) == len(identifiers), 'the length\
        of pl_inputs and identifiers should be the same'
    dataobject_list = []
    for i in range(len(pl_inputs)):
        assert '(' in pl_inputs[i] and pl_inputs[i].endswith(')')
        start_index = pl_inputs[i].index('(')
        objecttype = pl_inputs[i][0:start_index]
        parameters = '[' + pl_inputs[i][start_index + 1:-1] + ']'
        parameters = parameters.replace('\'', '')
        newline = re.sub(r"([\w/\.:\\]+)\b(?!['\"])", r'"\1"', parameters)
        jline = json.loads(newline)
        #the output of json.loads is in unicode format,convert to string
        jline = convert_unicode_str(jline)
        attributes = []
        for j in range(len(jline[0]) / 2):
            if isinstance(jline[0][j * 2 + 1], list):
                tmp = str(jline[0][j * 2 + 1]).replace('\'', '')
                tmp = tmp.replace(' ', '')
                attributes.append(tmp)
            else:
                attributes.append(jline[0][j * 2 + 1])
        identifier = identifiers[i]
        dataobject_list.append(module_utils.DataObject(
            objecttype, attributes, identifier))
    return dataobject_list


def dataobject2plstring(dataobject):
    identifier = dataobject.identifier
    pl_input = dataobject.objecttype
    if dataobject.attributes.strip():
        pl_input = pl_input + '(' + dataobject.attributes + ')'
    return pl_input, identifier


def convert_unicode_str(text):
    """convert a unicode string to string(u'DatasetId' to 'DatasetId')"""
    if isinstance(text, list):
        final = []
        for subtext in text:
            final.append(convert_unicode_str(subtext))
    else:
        final = str(text)
    return final


def convert_one_pipeline(jline):
    modules = []
    parameters = []
    pipeline = []
    for i in range(len(jline) / 2):
        modules.append(jline[i * 2])
        parameters.append(jline[i * 2 + 1])
    for i in range(len(parameters)):
        para_dict = dict()
        for j in range(len(parameters[i]) / 2):
            b = str(parameters[i][j * 2 + 1])
            # str(['neg', pos']), in the modules,need
            # to compare the parameters of the input objects
            # and the pipeline, string is more convient
            # to process
            if isinstance(parameters[i][j * 2 + 1], list):   # ['neg', pos']
                b = b.replace("\'", '')
                b = b.replace(' ', '')    # '[neg,pos]'
            para_dict[parameters[i][j * 2]] = b
        if 'status' in para_dict.keys():
            del para_dict['status']
        analysis = module_utils.Analysis(modules[i], para_dict)
        pipeline.append(analysis)
    return pipeline


def parse_text_pipeline(line):
    """parse a text line to a pipeline"""
    #convert to JSON object format
    line = line.replace('\'', '')
    if line.endswith('.'):
        line = line[0:-1]
    newline = re.sub(r"([\w/\.:\\]+)\b(?!['\"])", r'"\1"', line)
    jlines = json.loads(newline)
    #the output of json.loads is in unicode format,convert to string
    jlines = convert_unicode_str(jlines)
    pipelines = []
    if [True] * len(jlines) == [isinstance(jline, list) for jline in jlines]:
        for jline in jlines:
            #generate the module list and parameters list for this pipeline
            pipeline = convert_one_pipeline(jline)
            pipelines.append(pipeline)
    else:
        pipeline = convert_one_pipeline(jlines)
        pipelines.append(pipeline)
    return pipelines


def make_module_wd_name(module, module_name, pipeline_sequence, analysis, objects):
    single_object = module.get_identifier(analysis.parameters, objects)
    hash_string = module.make_unique_hash(
        single_object.identifier, pipeline_sequence,
        analysis.parameters)
    working_dir = module_name + '_BETSYHASH1_' + hash_string
    return working_dir


def make_module_wd(working_dir):
    os.mkdir(working_dir)


def copy_result_folder(working_dir, temp_dir, temp_outfile):
    """copy result files in temp folder to result folder,
      if fails, delete result folder"""
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    try:
        shutil.copyfile(os.path.join(temp_dir,'Betsy_parameters.txt'),
                        os.path.join(working_dir,'Betsy_parameters.txt'))
        if os.path.isdir(temp_outfile):
            shutil.copytree(os.path.join(temp_dir,temp_outfile),
                                     os.path.join(working_dir,temp_outfile))
        else:
            shutil.copyfile(os.path.join(temp_dir,temp_outfile),
                                     os.path.join(working_dir,temp_outfile))
        shutil.copyfile(os.path.join(temp_dir,'stdout.txt'),
                                     os.path.join(working_dir,'stdout.txt'))
    finally:
        if not os.path.isfile(os.path.join(working_dir,'stdout.txt')):
            shutil.rmtree(working_dir)
            raise


def run_module(analysis,objects,pipeline_sequence,clean_up=True):
    output_path = config.OUTPUTPATH
    assert os.path.exists(output_path), (
            'the output_path %s does not exist' % output_path)
    # get module
    module_name = analysis.name
    module = __import__('modules.' + module_name, globals(),
                        locals(), [module_name], -1)
    # make directory name
    working_dir = os.path.join(output_path, make_module_wd_name(
        module, module_name, pipeline_sequence, analysis, objects))
    # make name of outfile
    outfile = module.get_outfile(
            analysis.parameters, objects, pipeline_sequence)
    final_outfile = os.path.join(working_dir,os.path.split(outfile)[-1])
    temp_dir = ''
    temp_outfile = ''
    if not os.path.exists(os.path.join(working_dir,'stdout.txt')):
        #if no result has been generated, create temp folder and run analysis
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        try:
            temp_objects = module.run(
                analysis.parameters, objects, pipeline_sequence)
            if temp_objects:
                f = file(os.path.join(temp_dir, 'stdout.txt'),'w')
                f.write('This module runs successully.')
                f.close()
        except:
            logging.exception('Got exception on %s .run()'%module_name)
            raise
        temp_outfile = os.path.split(outfile)[-1]
    try:
        # found if the same analysis has been run and wait for it finished
        try:
            os.mkdir(working_dir)
        except OSError:
            try :
                module_timeout = int(config.MODULE_TIMEOUT)
            except AttributeError:
                module_timeout = 60
            i = 0
            while i <= module_timeout:
                if os.path.isfile(os.path.join(working_dir,'stdout.txt')):
                    break
                i = i + 1
                time.sleep(1)
        # if previous analysis not working, copy the current results to working_dir
        if (module_utils.exists_nz(temp_outfile) and not
            os.path.exists(os.path.join(working_dir,'stdout.txt'))):
            copy_result_folder(working_dir, temp_dir, temp_outfile)
    finally:
        if clean_up and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    os.chdir(working_dir)
    new_objects = module.get_newobjects(
        analysis.parameters, objects, pipeline_sequence)
    return final_outfile, new_objects


def run_pipeline(pipeline, objects, clean_up=True):
    """run a single pipeline"""
    OUTPUT = None
    outfile_list = []
    LOG_FILENAME = os.path.join(config.OUTPUTPATH, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)
    try:
        for  i in range(len(pipeline)):
            pipeline_sequence = [analysis.name for
                                 analysis in pipeline[0:i + 1]]
            analysis = pipeline[i]
            module_name = analysis.name
            print str(i + 1) + '.[' + time.strftime('%l:%M%p') + '] ' + module_name
            outfile,new_objects = run_module(
                analysis,objects,pipeline_sequence,clean_up=clean_up)
            outfile_list.append(outfile)
            objects = new_objects[:]
            if analysis == pipeline[-1]:
                if module_utils.exists_nz(outfile):
                    OUTPUT = outfile_list
                    print ('['+ time.strftime('%l:%M%p') +
                           '] Completed successfully and generated a file:')
                    print  OUTPUT[-1] + '\r'
                    print '\r'
                else:
                    print 'This pipeline has completed unsuccessfully'
                    raise ValueError(
                        'there is no output for this pipeline')
    except Exception, x:
            raise
    return OUTPUT


def make_pipelines(pl_output, pl_inputs):
    cwd = os.getcwd()
    pl_lib = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'rules')
    try:
        os.chdir(pl_lib)
        p = swipl.connect()
        file_names = os.listdir(pl_lib)
        for file_name in file_names:
            if file_name.endswith('.pl') and not file_name.startswith('._'):
                swipl.source_code(p, file_name)
        for one_input in pl_inputs:
            more_command = 'asserta(' + one_input + ').'
            swipl.send_query(p, more_command)
        results = swipl.send_query(p, pl_output)
        swipl.disconnect(p)
    finally:
        os.chdir(cwd)
    #get the variable name which show the pipeline information
    variable_name = pl_output.split(',')[-1][:-1].strip()
    pipelines = []
    for single_pipeline in results:
        text, variables = swipl.parse_solution(single_pipeline)
        assert variable_name in variables, ('the variable name %s does not '
                                             'exsit in the pipeline information'
                                             %variable_name)
        
        line = variables[variable_name]
        pipeline = parse_text_pipeline(line)
        pipelines.extend(pipeline)
    return pipelines
