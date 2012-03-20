#!/usr/bin/env python

#protocol_engine.py
import rule_engine
import argparse
import os
import protocol_utils



def filter_pipelines(protocol, inputs, output_info,
                     in_dataset_ids, in_contents, parameters=None):
    """given the Inputs and Output and Parameters dictionary,
       return a list of pipelines,if Parameters is None,
       then use the default parameters"""
    module = protocol_utils.import_protocol(protocol)
    pl_inputs = []
    assert len(in_dataset_ids) == len(inputs)
    assert len(in_contents) == len(inputs)
    for i in range(len(inputs)):
        x = module.predicate2arguments[inputs[i]]
        assert len(x) == 2
        in_parameters, modules = x
        query = protocol_utils.format_prolog_query(
            inputs[i], in_dataset_ids[i],
            in_contents[i], in_parameters,modules)
        pl_inputs.append(query)
    new_parameters = dict()
    if parameters:
        for key in parameters.keys():
            new_parameters[key.lower()] = parameters[key].lower()
    else:
        new_parameters = module.DEFAULT
    parameter_list = []
    for key in new_parameters.keys(): 
        parameter_list.extend([key,new_parameters[key]])
    test_dataset=None
    test_content=None
    if len(output_info) == 3:
        [output, out_dataset_id, out_content] = output_info
    if len(output_info) == 5:
        [output, out_dataset_id, out_content, test_dataset,test_content] = output_info
    pl_output = protocol_utils.format_prolog_query(
        output, out_dataset_id, out_content, parameter_list,'Modules',test_content)
    pipelines = rule_engine.make_pipelines(pl_output,pl_inputs)
    return pipelines
   
def run_protocol(protocol, inputs, output_info, identifiers,
                 in_dataset_ids, in_contents, parameters=None):
    """given the Inputs and Output and Parameters dictionary,
       run the pipelines,
       return a list of final result file,
      if Parameters is None, then use the default parameters"""
    
    module = protocol_utils.import_protocol(protocol)
    pipelines = filter_pipelines(protocol, inputs, output_info,
                                 in_dataset_ids, in_contents, parameters)
    pl_inputs = []
    for i in range(len(inputs)):
        x = module.predicate2arguments[inputs[i]]
        assert len(x) == 2
        in_parameters, modules = x
        query = protocol_utils.format_prolog_query(
            inputs[i], in_dataset_ids[i], in_contents[i],
            in_parameters, modules)
        pl_inputs.append(query)
    objects = rule_engine.plstring2dataobject(pl_inputs,identifiers)
    output_files_all = []
    parameters_all = []
    pipeline_sequence_all = []
    for pipeline in pipelines:
          out_files = rule_engine.run_pipeline(pipeline,objects)
          if out_files:
              pipeline_sequence = [analysis.name for analysis in pipeline]
              output_files_all.append(out_files[-1])
              parameters_all.append(pipeline[-1].parameters)
              pipeline_sequence_all.append(pipeline_sequence)
              
    return output_files_all,parameters_all,pipeline_sequence_all

def main():
    parser = argparse.ArgumentParser(
        description = 'run the protocol engine')
    parser.add_argument('--protocol',
                        dest = 'protocol',type = str,
                        help = 'The name of the protocol,eg. cluster_genes')
    parser.add_argument('--input',dest = 'input',type = str,
                        action = 'append',default = None,
                        help ='input:datasetid:contents:identifier_or_file\
                        or input:identifier')
    parser.add_argument('--output',dest = 'output',type=str,
                        help = 'datasetid:contents;datasetid:contents (train set;test set')
    parser.add_argument('--parameters',dest = 'parameters',
                        action = 'append',default = [],
                        type = str,help='key:value')
    parser.add_argument('--dry_run',dest = 'dry_run',const = True,
                        default = False,action = 'store_const',
                        help ='only shows the pipelines')
    parser.add_argument('--describe_protocol',
                        dest = 'describe_protocol',
                        action = 'store_const',default = False,
                        const = True,
                        help = 'shows the protocol details')

    args = parser.parse_args()
    if not args.protocol:
        raise parser.error('please specify the protocol')
    if not args.input:
        raise parser.error('please specify the input')
    len_inputs = []
    for i in args.input:
        len_input = len(i.split(':'))
        assert len_input in [2,4],'input is length 2 or 4'
        len_inputs.append(len_input)
        if len(i.split(':')) == 4:
            assert args.output,'please specify the output'
            if ';' in args.output:
                outputs = args.output.split(';')
                for output in outputs:
                    assert len(output.split(':')) == 2,'the\
                      format of output is datasetid:contents\
                      or datasetid:contents;datasetid:contents'
            else:
                assert len(args.output.split(':')) == 2,'the\
                      format of output is datasetid:contents'
                      
    assert len_inputs == [len_inputs[0]]*len(len_inputs),'the format of all inputs do not match'
    for i in args.parameters:
        assert len(i.split(':')) == 2, 'parameters should format like key:value'
        
    module = protocol_utils.import_protocol(args.protocol)
    
    inputs = []
    identifiers = []
    in_dataset_ids = []
    in_contents = []
    for i in args.input:
        inpair = i.split(':')
        if len(inpair) == 4:
            predicate, dataset_id, content, identifier = inpair
            in_contents.append(content)
            in_dataset_ids.append(dataset_id)
        elif len(inpair) == 2:
            predicate, identifier = inpair
            in_contents.append('not_applicable')
            in_dataset_ids.append('not_applicable')
        else:
            raise ValueError('the length of input should be 2 or 4')
        assert predicate in module.INPUTS,("%s is not recognized in %s)"
                                        %(str(predicate),args.protocol))
        inputs.append(predicate)
        identifiers.append(identifier)
    if args.output:
        if ';' in args.output:
            output = []
            outputs = args.output.split(';')
            for out in outputs:
                output.extend(out.split(':'))
        else:
            output = args.output.split(':')
    else:
        output = ['not_applicable','not_applicable']
 
    parameters = module.DEFAULT
    if args.parameters:
        for i in args.parameters:
            parpair = i.split(':')
            key,value = parpair
            assert key in module.PARAMETERS.keys(),(
                '%s is not a valid parameter key in %s'%(key,args.protocol))
            if not value.isdigit():
                assert value in module.PARAMETERS[key],(
                '   %s is not a valid parameter value in %s'%(value,args.protocol))
            parameters[key.lower()] = value.lower()
            
    if args.describe_protocol:
        print 'INPUTS', module.INPUTS
        print 'OUTPUTS', module.OUTPUTS
        print 'PARAMETERS', module.PARAMETERS
        print 'DEFAULTS', module.DEFAULT
        
    if args.dry_run:
        for output_file in module.OUTPUTS:
            output_info = [output_file]
            output_info.extend(output)
            pipelines = filter_pipelines(
            args.protocol, inputs,output_info,
            in_dataset_ids,in_contents,parameters)
            print len(pipelines)
            for pipeline in pipelines:
                for analysis in pipeline:
                    print analysis.name
                    print analysis.parameters
                print '-------------------------'
    else:
        final_output = []
        final_parameters = []
        final_pipeline_sequence = []
        for output_file in module.OUTPUTS:
            output_info = [output_file]
            output_info.extend(output)
            output_files,parameters_all,pipeline_sequence_all = run_protocol(
                args.protocol,inputs,output_info,
                identifiers,in_dataset_ids,in_contents,parameters)
            final_output.append(output_files)
            final_parameters.append(parameters_all)
            final_pipeline_sequence.append(pipeline_sequence_all)
        print final_output
        protocol_utils.get_result_folder(
            args.protocol,final_output,final_parameters,final_pipeline_sequence)
        
if __name__=='__main__':
    main()
    
