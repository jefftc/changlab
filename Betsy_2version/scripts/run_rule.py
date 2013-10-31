#! /usr/bin/env python
#run_rule.py
import os
import sys
import argparse
from Betsy import rule_engine_bie
from Betsy import bie
from Betsy import rulebase

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
    parser = argparse.ArgumentParser(description='Run the engine')
    parser.add_argument('--in_datatype', dest='in_datatype', default=None, action='append',
                        type=str, help='input data type for network')
    parser.add_argument('--out_datatype', dest='out_datatype', default=None, type=str,
                        help='out_datatype')
    parser.add_argument('--param', dest='param', default=[], type=str,action='append',
                        help='parameter in "key=value" format')
    parser.add_argument('--dry_run', dest='dry_run', action='store_const',
                        const=True, default=False,
                        help='show the test case procedure')
    parser.add_argument('--network', dest='network',type=str, default=None,
                        help='generate the new network file')
    args = parser.parse_args()
    in_datatype = args.in_datatype
    out_datatype = args.out_datatype
    assert in_datatype,'please specify the in_datatype'
    assert out_datatype,'please specify the out_datatype'
    param = args.param
    param = [ i.replace('=','="')+'"' for i in param]
    indices = []
    in_indices = get_indices('--in_datatype',sys.argv)
    out_index = get_indices('--out_datatype',sys.argv)
    indices= in_indices + out_index
    indices.sort()
    data_dict = {}
    for i in range(len(indices)):
        p_list = []
        if i == len(indices)-1:
            p_list = sys.argv[indices[i]:]
        else:
            p_list = sys.argv[indices[i]:indices[i+1]]
        param_indices = get_indices('--param',p_list)
        parameters = [p_list[j+1] for j in param_indices]
        parameters = [ k.replace('=','="')+'"' for k in parameters]
        data_dict[indices[i]]=(sys.argv[indices[i]+1],','.join(parameters))
    in_data = []
    for key in data_dict:
        datatype, attributes = data_dict[key]
        if key in in_indices:
            in_data.append(eval('rulebase.'+datatype +'('+attributes+')'))
        elif key in out_index:
            goal_datatype = eval('rulebase.' + datatype)
            goal_attributes = eval('dict(' + attributes + ')')
    print 'Generating network...'
    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)
    assert network, ('No pipeline has been generated, '
                       'please check your command.')
    if args.network:
        bie._plot_network_gv(args.network, network)
    if args.dry_run:
        bie._print_network(network)
    else:
        rule_engine_bie.run_pipeline(network,in_data)
        print 'All pipelines have completed successfully.'
     


if __name__ == '__main__':
    main()

