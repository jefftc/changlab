#! /usr/bin/env python
#run_rule.py
import os
import sys
import argparse
from Betsy import rule_engine_bie3
from Betsy import bie3
from Betsy import rulebase


def main():
    parser = argparse.ArgumentParser(description='Run the engine')
    parser.add_argument('--in_datatype', dest='in_datatype', default=None, action='append',
                        type=str, help='input data type for network')
    parser.add_argument('--identifier', dest='identifier', default=None, action='append',
                        type=str, help='input data identifier')
    parser.add_argument('--user_input', dest='user_input', default=None, action='append',
                        type=str, help='user input in key=value format')
    parser.add_argument('--out_datatype', dest='out_datatype', default=None, type=str,
                        help='out_datatype')
    parser.add_argument('--attr', dest='param', default=[], type=str,action='append',
                        help='attribute for datatype in "datatype,key=value" or "key=value" format')
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
    Attributes = []
    user_inputs = {}
    identifiers = []
    flag = None
    for i, arg in enumerate(sys.argv):
        if arg == "--in_datatype":
            assert len(sys.argv) > i+1
            in_datatypes.append(sys.argv[i+1])
            flag = 'in'
        elif arg == "--out_datatype":
            assert len(sys.argv) > i+1
            flag = 'out'
        elif arg == "--attr":
            assert len(sys.argv) > i+1
            if flag == 'in':
                assert in_datatypes
                x = sys.argv[i+1].split('=')
                key = x[0]
                value = x[1] 
                index = len(in_datatypes)-1
                if index not in in_parameters:
                    in_parameters[index] = {}
                in_parameters[index][key]=value
            elif flag == 'out':
                sub_datatype = sys.argv[i+1].split(',')[0]
                attr = ','.join(sys.argv[i+1].split(',')[1:])
                x = attr.split('=')
                key = x[0]
                value = x[1]
                fn = getattr(rulebase,sub_datatype)
                Attributes.append(bie3.Attribute(fn,key,value))
        elif arg == '--identifier':
            if len(in_datatypes)==len(identifiers)+1:
                identifiers.append(sys.argv[i+1])
            else:
                identifiers.extend(['']*(len(in_datatypes)-1))
                identifiers.append(sys.argv[i+1])
        elif arg == '--user_input':
            x = sys.argv[i+1].split('=')
            key = x[0]
            value = x[1]
            user_inputs[key]=value
    in_objects = []
    for i,in_datatype in enumerate(in_datatypes):
        fn = getattr(rulebase,in_datatype)
        in_data = fn.input()
        if in_parameters:
            if i in in_parameters:
                in_data = fn.input(**in_parameters[i])
        in_object = rule_engine_bie3.DataObject(in_data)
        if identifiers:
            in_object = rule_engine_bie3.DataObject(in_data,identifiers[i])
        in_objects.append(in_object)
    print 'Generating network...'
    network = bie3.backchain(rulebase.all_modules, goal_datatype, *Attributes)
    network = bie3.optimize_network(network)
    assert network, ('No pipeline has been generated, '
                       'please check your command.')
    if args.network:
        bie3.plot_network_gv(args.network, network)
    if args.network_text:
        handle = file(args.network_text,'w')
        try:
            bie3.print_network(network, handle)
        finally:
            handle.close()
    if args.dry_run:
        bie.print_network(network)
    else:
        rule_engine_bie3.run_pipeline(network,in_objects,user_inputs)
     


if __name__ == '__main__':
    main()

