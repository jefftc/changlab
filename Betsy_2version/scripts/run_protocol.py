#!/usr/bin/env python
#run_protocol.py
from Betsy import rule_engine_bie
import argparse
import os
import sys
from Betsy import protocol_utils
from Betsy import module_utils,rulebase,bie
import getpass

def main():
    parser = argparse.ArgumentParser(description='run the protocol engine')
    parser.add_argument('--protocol',
                        dest='protocol', type=str,
                        help='The name of the protocol,eg. cluster_genes')
    parser.add_argument('--in_datatype', dest='in_datatype', default=None, action='append',
                        type=str, help='input data type for network')
    parser.add_argument('--in_param', dest='param', default=[], type=str,action='append',
                        help='parameter in "key=value" format')
    parser.add_argument('--out_param', dest='out_param', default=[], type=str,action='append',
                        help='parameter in "key=value" format')
    
    parser.add_argument('--dry_run', dest='dry_run', const=True,
                        default=False, action='store_const',
                        help='only shows the pipelines')
    parser.add_argument('--network', dest='network',type=str, default=None,
                        help='generate the new network png file')
    parser.add_argument('--network_text', dest='network_text',type=str, default=None,
                        help='generate the output network text file')
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
    if not args.in_datatype:
        raise parser.error('please specify the in_datatype')
    out_param = args.out_param
    module = protocol_utils.import_protocol(args.protocol)
    protocol_utils.check_parameters(module.PARAMETERS)
    protocol_utils.check_default(module.PARAMETERS)
    assert args.in_datatype, 'please specify the in_datatype'
    for in_datatype in args.in_datatype:
        assert in_datatype in module.INPUTS, (
            "%s is not a recognized input for the %s protocol"
            % (str(in_datatype), args.protocol))
    goal_attributes = {}    
    parameters_dict = {}
    for parameter in module.PARAMETERS:
        parameters_dict[parameter.name]=parameter
        if parameter.default:
            goal_attributes[parameter.name]=parameter.default
    keys = parameters_dict.keys()   
    if args.out_param:
        for out in out_param:
             key,value = out.split('=')
             assert key in keys, (
                '%s is not a valid parameter key in %s' % (key, args.protocol)) 
             if parameters_dict[key].choices:
                 assert value in parameters_dict[key].choices, (
                    '%s is not a valid parameter value in %s'
                    % (value, args.protocol))  
                 goal_attributes[key] =  str(value)

    in_datatypes = []
    in_parameters = {}
    
    for i, arg in enumerate(sys.argv):
        if arg == "--in_datatype":
            assert len(sys.argv) > i+1
            in_datatypes.append(sys.argv[i+1])
        elif arg == "--in_param":
            assert len(sys.argv) > i+1
            assert in_datatypes
            index = len(in_datatypes)-1
            if index not in in_parameters:
                in_parameters[index] = {}
            x = sys.argv[i+1].split('=')
            key = x[0]
            value = x[1] 
            in_parameters[index][key]=value
    in_data = []
    for i,in_datatype in enumerate(in_datatypes):
        fn = getattr(rulebase,in_datatype)
        in_data.append(fn(**in_parameters[i]))
    goal_datatype = getattr(rulebase,module.OUTPUTS)
    if args.describe_protocol:
        print 'INPUTS', module.INPUTS   
        print 'OUTPUTS', module.OUTPUTS  
        print 'PARAMETERS', module.PARAMETERS
    print 'Generating network...'
    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    assert network, ('No network has been generated, '
                       'please check your command.')
    if args.network:
        print args.network
        bie.print_network(network)
        bie.plot_network_gv(args.network, network)
    if args.network_text:
        handle = file(args.network_text,'w')
        try:
            bie.print_network(network, handle)
        finally:
            handle.close()
    if args.dry_run:
        bie.print_network(network)
    else:
        rule_engine_bie.run_pipeline(network,in_data,args.user,args.job_name)
        print 'The network has completed successfully.'
       
        
if __name__ == '__main__':
    main()
