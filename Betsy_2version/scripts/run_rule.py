#! /usr/bin/env python
#run_rule.py
import os
import sys
import argparse
from Betsy import rule_engine_bie
from Betsy import bie
from Betsy import rulebase


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
                        help='generate the output network png file')
    parser.add_argument('--network_text', dest='network_text',type=str, default=None,
                        help='generate the output network text file')
    args = parser.parse_args()
    assert args.in_datatype,'please specify the in_datatype'
    assert args.out_datatype,'please specify the out_datatype'
    goal_datatype = getattr(rulebase, args.out_datatype)
    in_datatypes = []
    in_parameters = {}
    goal_attributes = {}
    flag = None
    for i, arg in enumerate(sys.argv):
        if arg == "--in_datatype":
            assert len(sys.argv) > i+1
            in_datatypes.append(sys.argv[i+1])
            flag = 'in'
        if arg == "--out_datatype":
            assert len(sys.argv) > i+1
            flag = 'out'
        elif arg == "--param":
            assert len(sys.argv) > i+1
            x = sys.argv[i+1].split('=')
            key = x[0]
            value = x[1] 
            if flag == 'in':
                assert in_datatypes
                index = len(in_datatypes)-1
                if index not in in_parameters:
                    in_parameters[index] = {}
                in_parameters[index][key]=value
            elif flag == 'out':
                assert goal_datatype
                goal_attributes[key]=value
    in_data = []
    for i,in_datatype in enumerate(in_datatypes):
        fn = getattr(rulebase,in_datatype)
        in_data.append(fn(**in_parameters[i]))
    print 'Generating network...'
    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.select_start_node(network, in_data)
    network = bie.optimize_network(network)
    assert network, ('No pipeline has been generated, '
                       'please check your command.')
    if args.network:
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
        rule_engine_bie.run_pipeline(network,in_data)
        print 'All pipelines have completed successfully.'
     


if __name__ == '__main__':
    main()

