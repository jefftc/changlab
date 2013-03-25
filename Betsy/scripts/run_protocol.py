#!/usr/bin/env python
#run_protocol.py
from Betsy import rule_engine
import argparse
import os
from Betsy import protocol_utils
from Betsy import module_utils
import getpass

def filter_pipelines(protocol, inputs, in_contents, output, parameters):
    """given the Inputs and Output and Parameters dictionary,
       return a list of pipelines,if Parameters is None,
       then use the default parameters"""
    module = protocol_utils.import_protocol(protocol)
    pl_inputs = []
    assert len(in_contents) == len(inputs)
    for i in range(len(inputs)):
        x = module.predicate2arguments[inputs[i]]
        assert len(x) == 2
        in_parameters = x[0][:]
        in_parameters.extend(['contents', '[' + in_contents[i] + ']'])
        modules = x[1]
        query = protocol_utils.format_prolog_query(
            inputs[i], in_parameters, modules)
        pl_inputs.append(query)
    parameter_list = []
    for key in parameters.keys():
        parameter_list.extend([key, parameters[key]])
    pl_output = protocol_utils.format_prolog_query(
        output, parameter_list, 'Modules')
    print 'Generating pipelines...'
    pipelines = rule_engine.make_pipelines(pl_output, pl_inputs)
    assert pipelines, ('No pipeline has been generated, '
                                'please check your command.')
    print '%d pipelines have been generated.\n' %len(pipelines)
    return pipelines



def run_protocol(protocol, inputs, output, identifiers,
                 in_contents, parameters,user=getpass.getuser(),job_name='',clean_up=True):
    """given the Inputs and Output and Parameters dictionary,
       run the pipelines,
       return a list of final result file,
      if Parameters is None, then use the default parameters"""
    module = protocol_utils.import_protocol(protocol)
    pipelines = filter_pipelines(protocol, inputs, in_contents, output,
                                 parameters)
    pl_inputs = []
    for i in range(len(inputs)):
        x = module.predicate2arguments[inputs[i]]
        assert len(x) == 2
        in_parameters = x[0][:]
        in_parameters.extend(['contents', '[' + in_contents[i] + ']'])
        modules = x[1]
        query = protocol_utils.format_prolog_query(
            inputs[i], in_parameters, modules)
        pl_inputs.append(query)
    objects = rule_engine.plstring2dataobject(pl_inputs, identifiers)
    output_files_all = []
    parameters_all = []
    pipeline_sequence_all = []
    k = 1
    for pipeline in pipelines:
        pipeline_label = pipeline[-1].parameters['label']
        pipeline_label = pipeline_label.replace('_',' ')
        pipeline = pipeline[:-1]
        print  'Pipeline' + str(k) + ': '+ pipeline_label + '.\r'
        out_files = rule_engine.run_pipeline(pipeline, objects,user=user,job_name=job_name,clean_up=clean_up)
        k = k + 1
        if out_files:
            pipeline_sequence = [analysis.name for analysis in pipeline]
            output_files_all.append(out_files[-1])
            parameters_all.append(pipeline[-1].parameters)
            pipeline_sequence_all.append(pipeline_sequence)
    print 'All pipelines have completed successfully.'
    return output_files_all, parameters_all, pipeline_sequence_all


def main():
    parser = argparse.ArgumentParser(description='run the protocol engine')
    parser.add_argument('--protocol',
                        dest='protocol', type=str,
                        help='The name of the protocol,eg. cluster_genes')
    parser.add_argument('--input', dest='input', type=str,
                        action='append', default=None,
                        help='input:contents:identifier_or_file\
                        or input:identifier')
    parser.add_argument('--parameters', dest='parameters',
                        action='append', default=[],
                        type=str, help='key:value')
    parser.add_argument('--dry_run', dest='dry_run', const=True,
                        default=False, action='store_const',
                        help='only shows the pipelines')
    parser.add_argument('--describe_protocol',
                        dest='describe_protocol',
                        action='store_const', default=False,
                        const=True,
                        help='shows the protocol details')
    parser.add_argument('--user',
                        dest='user', default=getpass.getuser(),
                        type=str,
                        help='the username who run the command')
    parser.add_argument('--job_name',
                        dest='job_name', default='',
                        type=str,
                        help='the name of this job')
    parser.add_argument('--dont_cleanup',
                        dest='clean_up',
                        action='store_const', default=True,
                        const=False,
                        help='do not clean up the temp folder')
    args = parser.parse_args()
    if not args.protocol:
        raise parser.error('please specify the protocol')
    if not args.input:
        raise parser.error('please specify the input')
    len_inputs = []
    for i in args.input:
        len_input = len(i.split(':'))
        assert len_input in [2, 3], 'input is length 2 or 3'
        len_inputs.append(len_input)
    for i in args.parameters:
        assert len(i.split(':')) == 2, (
            'parameters should format like key:value')
    module = protocol_utils.import_protocol(args.protocol)
    protocol_utils.check_parameters(module.PARAMETERS)
    protocol_utils.check_default(module.DEFAULT, module.PARAMETERS)
    inputs = []
    identifiers = []
    in_dataset_ids = []
    in_contents = []
    for i in args.input:
        inpair = i.split(':')
        if len(inpair) == 3:
            predicate, content, identifier = inpair
        elif len(inpair) == 2:
            predicate, identifier = inpair
            content = 'unknown'
        else:
            raise ValueError('--input is in the wrong format')
        assert predicate in module.INPUTS, (
            "%s is not a recognized input for the %s protocol"
            % (str(predicate), args.protocol))
        if os.path.exists(identifier):
            identifier = os.path.realpath(identifier)
        inputs.append(predicate)
        in_contents.append(content)
        identifiers.append(identifier)
    parameters = module.DEFAULT
    if args.parameters:
        for i in args.parameters:
            parpair = i.split(':')
            key, value = parpair
            assert key in module.PARAMETERS.keys(), (
                '%s is not a valid parameter key in %s' % (key, args.protocol))
            if module.PARAMETERS[key] == 'float':
                assert module_utils.is_number(value), (
                    '%s is not a number' % value)
                parameters[key] = str(value)
            elif module.PARAMETERS[key] == 'integer':
                assert value.isdigit(), '%s is not a number' % value
                parameters[key] = str(value)
            elif module.PARAMETERS[key] == 'list':
                parameters[key] = '[' + value + ']'
            elif module.PARAMETERS[key] == 'string':
                parameters[key] = '\'' + value + '\''
            else:
                assert value in module.PARAMETERS[key], (
                    '%s is not a valid parameter value in %s'
                    % (value, args.protocol))
                parameters[key] = str(value)
    if args.describe_protocol:
        print 'INPUTS', module.INPUTS
        print 'OUTPUTS', module.OUTPUTS
        print 'PARAMETERS', module.PARAMETERS
        print 'DEFAULTS', module.DEFAULT
    if args.dry_run:
        p_n = 1
        output_file = module.OUTPUTS
        pipelines = filter_pipelines(
            args.protocol, inputs, in_contents, output_file,
            parameters)
        for pipeline in pipelines:
            print 'pipeline' + str(p_n), '\r'
            p_n = p_n + 1
            k = 1
            for analysis in pipeline:
                print str(k) + '.', analysis.name, '\r'
                print 'parameters', analysis.parameters
                k = k + 1
            print '------------------------'
    else:
        final_output = []
        final_parameters = []
        final_pipeline_sequence = []
        output = module.OUTPUTS
        output_file, parameters_all, pipeline_sequence_all = run_protocol(
            args.protocol, inputs, output, identifiers, in_contents,
            parameters,user=args.user,job_name = args.job_name,clean_up=args.clean_up)
        assert output_file, ('No output file is generated and '
                             'no report can be generated')
        module1 = protocol_utils.import_protocol(args.protocol)
        print module1.OUTPUTS
        module = __import__('Betsy.modules.' + module1.OUTPUTS, globals(),
                            locals(), [module1.OUTPUTS], -2)
        module.run(output_file, parameters_all,pipeline_sequence_all, args.user,args.job_name)

        
if __name__ == '__main__':
    main()
