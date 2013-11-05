#!/usr/bin/env python
#run_protocol.py
from Betsy import rule_engine_bie
import argparse
import os
import sys
from Betsy import protocol_utils
from Betsy import module_utils,rulebase,bie
import getpass

def get_indices(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices

def main():
    parser = argparse.ArgumentParser(description='run the protocol engine')
    parser.add_argument('--protocol',
                        dest='protocol', type=str,
                        help='The name of the protocol,eg. cluster_genes')
    parser.add_argument('--in_datatype', dest='in_datatype', default=None, action='append',
                        type=str, help='input data type for network')
    parser.add_argument('--param', dest='param', default=[], type=str,action='append',
                        help='parameter in "key=value" format')
    parser.add_argument('--out_param', dest='out_param', default=[], type=str,action='append',
                        help='parameter in "key=value" format')
    
    parser.add_argument('--dry_run', dest='dry_run', const=True,
                        default=False, action='store_const',
                        help='only shows the pipelines')
    parser.add_argument('--network', dest='network',type=str, default=None,
                        help='generate the new network file')
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
        
    parameters = {}    
    PARAMETERS_dict = dict()
    for parameter in module.PARAMETERS:
        PARAMETERS_dict[parameter.name]=parameter
        if parameter.default:
            parameters[parameter.name]=parameter.default
    keys = PARAMETERS_dict.keys()   
    if args.out_param:
        for out in out_param:
             key,value = out.split('=')
             assert key in keys, (
                '%s is not a valid parameter key in %s' % (key, args.protocol)) 
             if PARAMETERS_dict[key].choices:
                 assert value in PARAMETERS_dict[key].choices, (
                    '%s is not a valid parameter value in %s'
                    % (value, args.protocol))  
                 parameters[key] =  str(value)
    indices = get_indices('--in_datatype',sys.argv)
    indices.sort()
    data_dict = {}
    for i in range(len(indices)):
        p_list = []
        if i == len(indices)-1:
            p_list = sys.argv[indices[i]:]
        else:
            p_list = sys.argv[indices[i]:indices[i+1]]
        param_indices = get_indices('--param',p_list)
        in_parameters = [p_list[j+1] for j in param_indices]
        in_parameters = [ k.replace('=','="')+'"' for k in in_parameters]
        data_dict[indices[i]]=(sys.argv[indices[i]+1],','.join(in_parameters))
    in_data = []
    for key in data_dict:
        datatype, attributes = data_dict[key]
        in_data.append(eval('rulebase.'+datatype +'('+attributes+')'))
    goal_datatype = eval('rulebase.'+module.OUTPUTS)
    goal_attributes = parameters
    if args.describe_protocol:
        print 'INPUTS', module.INPUTS   
        print 'OUTPUTS', module.OUTPUTS  
        print 'PARAMETERS', module.PARAMETERS
    print 'Generating network...'
##    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
##    network = bie.optimize_network(network)
##    network = bie.prune_network_by_start(network, in_data)
    import pickle
    f=file('network_classify_report','rb')
    network = pickle.load(f)
    f.close()
    assert network, ('No network has been generated, '
                       'please check your command.')
    if args.network:
        print args.network
        bie._print_network(network)
        bie._plot_network_gv(args.network, network)
    if args.dry_run:
        bie._print_network(network)
    else:
        rule_engine_bie.run_pipeline(network,in_data)
        print 'The network has completed successfully.'
       
        
if __name__ == '__main__':
    main()
