#protocol_engine.py

import rule_engine
import argparse


def filter_pipelines(protocol, inputs, output_info,
                     in_dataset_ids, in_contents, parameters=None):
    """given the Inputs and Output and Parameters dictionary,
       return a list of pipelines,if Parameters is None,
       then use the default parameters"""
    module = __import__(protocol)
    pl_inputs = []
    assert len(in_dataset_ids) == len(inputs)
    assert len(in_contents) == len(inputs)
    for i in range(len(inputs)):
        x = module.predicate2arguments[inputs[i]]
        assert len(x) == 2
        in_parameters, modules = x
        query = module.format_prolog_query(
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
    [output, out_dataset_id, out_content] = output_info
    pl_output = module.format_prolog_query(
        output, out_dataset_id, out_content, parameter_list,'Modules')
    pipelines = rule_engine.make_pipelines(pl_output,pl_inputs)
    return pipelines
   
def run_protocol(protocol, inputs, output_info, identifiers,
                 in_dataset_ids, in_contents, parameters=None):
    """given the Inputs and Output and Parameters dictionary,
       run the pipelines,
       return a list of final result file,
      if Parameters is None, then use the default parameters"""
    module = __import__(protocol)
    pipelines = filter_pipelines(protocol, inputs, output_info,
                                 in_dataset_ids, in_contents, parameters)
    pl_inputs = []
    for i in range(len(inputs)):
        x = module.predicate2arguments[inputs[i]]
        assert len(x) == 2
        in_parameters, modules = x
        query = module.format_prolog_query(
            inputs[i], in_dataset_ids[i], in_contents[i],
            in_parameters, modules)
        pl_inputs.append(query)
    objects = rule_engine.plstring2dataobject(pl_inputs,identifiers)
    output_files = []
    for pipeline in pipelines:
          output_file = rule_engine.run_pipeline(pipeline,objects)
          if output_file:
              output_files.append(output_file)
    return output_files

def main():
    parser = argparse.ArgumentParser(
        description = 'run the protocol engine')
    parser.add_argument('--protocol_file',
                        dest = 'protocol_file',type = str,
                        help = 'input the protocol file')
    parser.add_argument('--input',dest = 'input',type = str,
                        action = 'append',default = None,
                        help ='input:datasetid:contents:identifier_or_file\
                        or input:identifier')
    parser.add_argument('--output',dest = 'output',type=str,
                        help = 'the output predicate:datasetid:contents\
                        or the output predicate')
    parser.add_argument('--parameters',dest = 'parameters',
                        action = 'append',default = None,
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
    if not args.protocol_file:
        raise ValueError('please specify the protocol_file')
    module = __import__(args.protocol_file)
    inputs = []
    identifiers = []
    in_dataset_ids = []
    in_contents = []
    output_info = []
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
        inputs.append(predicate)
        identifiers.append(identifier)
    
    output_info = args.output.split(':')
    assert len(output_info) == len(inpair) - 1
    if len(output_info) == 1:
        output_info.extend(['not_applicable','not_applicable'])
        
    parameters = module.DEFAULT
    if args.parameters:
        for i in args.parameters:
            parpair = i.split(':')
            key,value = parpair
            parameters[key.lower()] = value.lower()
            
    if args.describe_protocol:
        print 'INPUTS', module.INPUTS
        print 'OUTPUTS', module.OUTPUTS
        print 'PARAMETERS', module.PARAMETERS
        print 'DEFAULTS', module.DEFAULT
        
    if args.dry_run:
        pipelines = filter_pipelines(
        args.protocol_file, inputs,output_info,
        in_dataset_ids,in_contents,parameters)
        print len(pipelines)
        for pipeline in pipelines:
            for analysis in pipeline:
                print analysis.name
                print analysis.parameters
            print '-------------------------'
    else:
        output_files = run_protocol(
            args.protocol_file,inputs,output_info,
            identifiers,in_dataset_ids,in_contents,parameters)
        print output_files
    
    
        
if __name__=='__main__':
    main()
    
