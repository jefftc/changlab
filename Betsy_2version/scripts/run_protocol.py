#!/usr/bin/env python
#run_protocol.py
from Betsy import rule_engine_bie
import argparse
import os
from Betsy import protocol_utils
from Betsy import module_utils,rulebase,bie
import getpass


def main():
    parser = argparse.ArgumentParser(description='run the protocol engine')
    parser.add_argument('--protocol',
                        dest='protocol', type=str,
                        help='The name of the protocol,eg. cluster_genes')
    parser.add_argument('--input', dest='input', type=str,
                        action='append', default=None,
                        help='input:identifier_or_file\
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
        assert len_input == 2, 'input is length 2'
        len_inputs.append(len_input)
    for i in args.parameters:
        assert len(i.split(':')) == 2, (
            'parameters should format like key:value')
    module = protocol_utils.import_protocol(args.protocol)
    protocol_utils.check_parameters(module.PARAMETERS)
    protocol_utils.check_default(module.PARAMETERS)
    inputs = []
    identifiers = []
    in_dataset_ids = []
    in_contents = []
    for i in args.input:
        inpair = i.split(':')
        predicate, identifier = inpair
        assert predicate in module.INPUTS, (
            "%s is not a recognized input for the %s protocol"
            % (str(predicate), args.protocol))
        if os.path.exists(identifier):
            identifier = os.path.realpath(identifier)
        inputs.append(predicate)
        identifiers.append(identifier)
    parameters = dict()
    PARAMETERS_dict = dict()
    for parameter in module.PARAMETERS:
        PARAMETERS_dict[parameter.name]=parameter
        if parameter.default:
            parameters[parameter.name]=parameter.default
    keys = PARAMETERS_dict.keys()
    if args.parameters:
        for i in args.parameters:
            parpair = i.split(':')
            key, value = parpair
            assert key in keys, (
                '%s is not a valid parameter key in %s' % (key, args.protocol)) 
            if PARAMETERS_dict[key].type == 'float':   
                assert module_utils.is_number(value), (
                    '%s is not a number' % value)
                parameters[key] = str(value)
            elif PARAMETERS_dict[key].type == 'integer':
                assert value.isdigit(), '%s is not a number' % value
                parameters[key] = str(value)
            elif PARAMETERS_dict[key].type == 'list':
                parameters[key] = '[' + value + ']'
            elif PARAMETERS_dict[key].type == 'string':
                parameters[key] = '\'' + value + '\''
            else:
                assert value in PARAMETERS_dict[key].choices, (
                    '%s is not a valid parameter value in %s'
                    % (value, args.protocol))
                parameters[key] = str(value)
    in_datas = []
    for i,in_data in enumerate(inputs):
       in_datas.append(eval('rulebase.'+in_data+'(filename="'+identifiers[i]+'")'))
    goal_datatype = eval('rulebase.'+module.OUTPUTS)
    goal_attributes = parameters
    if args.describe_protocol:
        print 'INPUTS', module.INPUTS   
        print 'OUTPUTS', module.OUTPUTS  
        print 'PARAMETERS', module.PARAMETERS
    print 'Generating network...'
    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_datas)
    assert network, ('No pipeline has been generated, '
                       'please check your command.')
    if args.dry_run:
        bie._print_network(network)
    else:
        rule_engine_bie.run_pipeline(network,in_datas)
        print 'All pipelines have completed successfully.'

        
if __name__ == '__main__':
    main()
