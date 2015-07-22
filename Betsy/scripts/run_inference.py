#! /usr/bin/env python
#run_inference.py
import os
import sys
import argparse
from Betsy import bie3
from Betsy import rulebase
import itertools

# parents_of
# combinations
# get_input_nodes
# print_start_nodes

def parents_of(network,node_id):
    parent_nodes = []
    for key in network.transitions:
       if node_id in network.transitions[key]:
            parent_nodes.append(key)
    return parent_nodes

def combinations(*start_ids):
    start_nodes = []
    for ids in start_ids:
        if len(start_nodes)>=len(ids):
            start_nodes.extend([zip(x,ids)for x in itertools.permutations(
                start_nodes,len(ids))])
        else:
            start_nodes.extend([zip(x,start_nodes) for x in itertools.permutations(
                        ids,len(start_nodes))])
    return start_nodes


def get_input_nodes(network,out_id,start_nodes):
    node = network.nodes[out_id]
    parent_nodes = parents_of(network,out_id)
    if isinstance(node,bie3.Data):
        if out_id not in start_nodes:
            start_nodes.append(out_id)
        for module_id in parent_nodes:
            get_input_nodes(network,module_id,start_nodes)
    elif isinstance(node,bie3.Module):
        in_data_num = len(network.nodes[out_id].in_datatypes)
        if in_data_num == 1:
            for parent_id in parent_nodes:
                get_input_nodes(network,parent_id,start_nodes)
        elif in_data_num >= 2:
            #need to consider the case when there are more than the required input point to the module
            assert len(parent_nodes)==in_data_num,'input data number is different from the module require'
            start_ids = []
            for data_id in parent_data:
                ids = []
                get_input_nodes(network,data_id,ids)
                start_ids.append(ids)
            assert len(start_ids)>=2
            temp_parent_start_nodes = combinations(*start_ids)
            start_nodes.extend(temp_parent_start_nodes)

            

def print_start_nodes(network):
    start_nodes = [0]
    get_input_nodes(network,0,start_nodes)
    start_nodes.sort()
    for i in start_nodes:
        if isinstance(i,int):
            print i,network.nodes[i]
        elif isinstance(i,tuple):
            for node_id in i:
                print node_id,network.nodes[node_id]
    return start_nodes

    
def main():
    parser = argparse.ArgumentParser(description='Run the inference')
    parser.add_argument('--out_datatype', dest='out_datatype', default=None, type=str,
                        help='out_datatype')
    parser.add_argument('--in_datatype', dest='in_datatype', default=None, action='append',
                        type=str, help='input data type for network')
    parser.add_argument('--attr', dest='attr', default=[], type=str,action='append',
                        help='attribute for datatype in "datatype,key=value" or "key=value" format')
    parser.add_argument('--no_optimization', dest='no_optimization', action='store_const',
                        const=True, default=False,
                        help='do not optimization the network')
    parser.add_argument('--plot_network', dest='network',type=str, default=None,
                        help='generate the output network png file')
    parser.add_argument('--print_start_nodes', dest='print_start_nodes',action='store_const',
                        default=False,const=True,
                        help='print start nodes')
    args = parser.parse_args()
    assert args.out_datatype,'please specify the out_datatype'
    goal_datatype = getattr(rulebase, args.out_datatype)
    in_datatypes = []
    in_parameters = {}
    flag = None
    Attributes = []
    for i, arg in enumerate(sys.argv):
        if arg == "--in_datatype":
            assert len(sys.argv) > i+1
            in_datatypes.append(sys.argv[i+1])
            flag = 'in'
        if arg == "--out_datatype":
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
                attr = sys.argv[i+1].split(',')[1]
                x = attr.split('=')
                key = x[0]
                value = x[1]
                fn = getattr(rulebase,sub_datatype)
                Attributes.append(bie3.Attribute(fn,key,value))
    network = bie3.backchain(rulebase.all_modules,goal_datatype,*Attributes)
 
    
    in_data = []
    if in_data:
        for i,in_datatype in enumerate(in_datatypes):
            fn = getattr(rulebase,in_datatype)
            in_data.append(fn(**in_parameters[i]))
 
    if args.in_datatype:
        pass
        #network = bie3.select_start_node(network, in_data)
    if not args.no_optimization:
        network = bie3.optimize_network(network)
        
    assert network, ('No pipeline has been generated, '
                       'please check your command.')
    
    
    if args.network:
        bie3.plot_network_gv(args.network, network)
    if args.print_start_nodes:
        print_start_nodes(network)
    else: 
        bie3.print_network(network)
    


if __name__ == '__main__':
    main()

