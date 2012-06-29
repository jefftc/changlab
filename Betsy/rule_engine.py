#! /usr/bin/env python
#rule_engine.py
import os
import swipl
import shutil
import sys
import logging
import json
import argparse
import Betsy_config
import module_utils
import re
class Analysis:
    def __init__(self,name,parameters):
        self.name = name
        self.parameters = parameters
        
class DataObject:
    def __init__(self,objecttype,attributes,identifier):
        self.objecttype = objecttype
        self.attributes = attributes
        self.identifier = identifier

          
def plstring2dataobject(pl_inputs,identifiers):
    assert len(pl_inputs) == len(identifiers),'the length\
           of pl_inputs and identifiers should be the same'
    dataobject_list = []
    for i in range(len(pl_inputs)):
        assert '(' in pl_inputs[i] and pl_inputs[i].endswith(')')
        start_index = pl_inputs[i].index('(')
        objecttype = pl_inputs[i][0:start_index]
        parameters = '['+pl_inputs[i][start_index+1:-1]+']'
        parameters = parameters.replace('\'','')
        newline = re.sub(r"([\w/\.:\\]+)\b(?!['\"])", r'"\1"', parameters)
        jline = json.loads(newline)
        #the output of json.loads is in unicode format,convert to string
        jline=convert_unicode_str(jline)
        attributes = []
        for j in range(len(jline[0])/2):
            if isinstance(jline[0][j*2+1],list):
                tmp = str(jline[0][j*2+1]).replace('\'','')
                tmp = tmp.replace(' ','')
                attributes.append(tmp)
            else:
                attributes.append(jline[0][j*2+1])
        identifier = identifiers[i]
        dataobject_list.append(DataObject(objecttype,attributes,identifier))
    return dataobject_list

def dataobject2plstring(dataobject):
    identifier = dataobject.identifier
    pl_input = dataobject.objecttype
    if dataobject.attributes.strip():
        pl_input = pl_input + '(' + dataobject.attributes + ')'
    return pl_input,identifier

def convert_unicode_str(text):
    """convert a unicode string to string(u'DatasetId' to 'DatasetId')"""
    if isinstance(text,list):
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
    for i in range(len(jline)/2):
        modules.append(jline[i*2])
        parameters.append(jline[i*2+1])
    for i in range(len(parameters)):
        para_dict = dict()
        for j in range(len(parameters[i])/2):
            b=str(parameters[i][j*2+1])    #str(['neg', pos']), in the modules,need
                                           #to compare the parameters of the input objects 
                                           #and the pipeline, string is more convient
                                           #to process
                                          
            if isinstance(parameters[i][j*2+1],list): #['neg', pos']
                b = b.replace("\'",'')
                b = b.replace(' ','')    #'[neg,pos]'
            para_dict[parameters[i][j*2]] = b
        if 'status' in para_dict.keys():
            del para_dict['status']
        analysis = Analysis(modules[i],para_dict)
        pipeline.append(analysis)
    return pipeline

def parse_text_pipeline(line):
    """parse a text line to a pipeline"""
    #convert to JSON object format
    line = line.replace('\'','')
    newline = re.sub(r"([\w/\.:\\]+)\b(?!['\"])", r'"\1"', line)
    jlines = json.loads(newline)
    #the output of json.loads is in unicode format,convert to string
    jlines=convert_unicode_str(jlines)
    pipelines =[]
    if [True]*len(jlines)==[isinstance(jline,list) for jline in jlines]:
        for jline in jlines:
            #generate the module list and parameters list for this pipeline
            pipeline = convert_one_pipeline(jline)
            pipelines.append(pipeline)
    else:
        pipeline = convert_one_pipeline(jlines)
        pipelines.append(pipeline)
    return pipelines

def make_module_wd_name(module_name,hash_string):
    working_dir = module_name + '_BETSYHASH1_' + hash_string
    return working_dir

def make_module_wd(working_dir):
    os.mkdir(working_dir)
    
    
def run_pipeline(pipeline,objects):
    """run a single pipeline"""
    cwd = os.getcwd()
    OUTPUT = None
    outfile_list = []
    try:
        output_path = Betsy_config.OUTPUTPATH
        assert os.path.exists(output_path), 'the output_path does not exist'
        for  i in range(len(pipeline)):
            pipeline_sequence = [analysis.name for analysis in pipeline[0:i+1]]
            analysis = pipeline[i]
            module_name = analysis.name
            print str(i+1) + '.' + module_name
            module = __import__(module_name)
            single_object = module.get_identifier(analysis.parameters,objects)
            hash_string = module.make_unique_hash(single_object.identifier,pipeline_sequence,analysis.parameters)
            working_dir = os.path.join(output_path,make_module_wd_name(
                                        module_name,hash_string))
            if not os.path.exists(working_dir):
                make_module_wd(working_dir)
            try:
                os.chdir(working_dir)
                outfile = module.get_outfile(
                    analysis.parameters,objects,pipeline_sequence)
                if module_utils.exists_nz(outfile):
                    new_objects = module.get_newobjects(
                        analysis.parameters,objects,pipeline_sequence)
                elif not module_utils.exists_nz(outfile):
                    new_objects = module.run(
                        analysis.parameters,objects,pipeline_sequence)
                    if not new_objects:
                        break
                outfile_list.append(outfile)
                objects = new_objects[:]
                if analysis == pipeline[-1]:
                    if module_utils.exists_nz(outfile):
                        OUTPUT = outfile_list
                        print 'File: ', OUTPUT[-1] + '\r'
                        print '\r'
                    else:
                        raise ValueError('there is no output for this pipeline')
            finally:
                os.chdir(cwd)
    
    except Exception,x:
            raise
    
    return OUTPUT

def make_pipelines(pl_output,pl_inputs):
    cwd = os.getcwd()
    pl_lib = os.path.join(os.path.dirname(__file__),'rules')
    try:
        os.chdir(pl_lib)
        p = swipl.connect()
        swipl.source_code(p,'signal_file.pl')
        swipl.source_code(p,'signal_raw.pl')
        swipl.source_code(p,'signal_clean.pl')
        swipl.source_code(p,'signal_norm1.pl')
        swipl.source_code(p,'signal_norm2.pl')
        swipl.source_code(p,'prolog_utils.pl')
        swipl.source_code(p,'clustering.pl')
        swipl.source_code(p,'classification.pl')
        swipl.source_code(p,'class_label_file.pl')
        swipl.source_code(p,'differentail_expressed_gene_analysis.pl')
        swipl.source_code(p,'signature_analysis.pl')
        swipl.source_code(p,'comparison_report.pl')
        for one_input in pl_inputs:
            more_command = 'asserta(' + one_input +').'
            swipl.send_query(p,more_command)
        results = swipl.send_query(p,pl_output)
        swipl.disconnect(p)
    finally:
        os.chdir(cwd)
    #get the variable name which show the pipeline information
    variable_name = pl_output.split(',')[-1][:-1]  
    pipelines = []
    for single_pipeline in results:
        text,variables = swipl.parse_solution(single_pipeline)
        line = variables[variable_name]
        pipeline = parse_text_pipeline(line)
        pipelines.extend(pipeline)
    return pipelines

def main():
    
      parser = argparse.ArgumentParser(description='Run the engine')
      
      parser.add_argument('--plin',dest='plin',default=None,action='append',type=str,
                          help='prolog command for pipeline start')
      parser.add_argument('--plout',dest='plout',default=None,type=str,
                          help='prolog command for querying')
      parser.add_argument('--id',dest='id',default=None,action='append',type=str,
                          help='sepecify the file or database id')
      parser.add_argument('--dry_run',dest='dry_run',action='store_const',
                          const=True,default=False,
                          help='show the test case procedure')
      
      args = parser.parse_args()
      assert len(args.id) == len(args.plin),'the length of id and plin should be the equal'
      pl_inputs = args.plin
      pl_output = args.plout
      identifiers = args.id
      for i in range(len(identifiers)):
          if os.path.exists(identifiers[i]):
              identifiers[i] = os.path.realpath(identifiers[i])
      objects = plstring2dataobject(pl_inputs,identifiers)
      pipelines = make_pipelines(pl_output,pl_inputs)
      
      if args.dry_run:
            print len(pipelines)
            for pipeline in pipelines:
                  for analysis in pipeline:
                        print 'module',analysis.name,'\r'
                        print 'parameters',analysis.parameters
                  print '------------------------'
            print len(pipelines)
      else:
            k = 1
            for pipeline in pipelines:
                print  'pipeline' + str(k) + ':','\r'
                run_pipeline(pipeline,objects)
                print '\r'
                k = k + 1
                
                
                  
if __name__=='__main__':
    main()
