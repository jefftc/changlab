#! /usr/bin/env python
#run_rule.py
import sys
import os
import shutil
import getpass
import argparse
from Betsy import rule_engine_bie3
from Betsy import bie3
from Betsy import rulebase
from Betsy import userfile




def get_all_user_inputs():
    modules = rulebase.all_modules
    user_inputs = []
    for module in modules:
        for user_input in module.user_inputs:
            if user_input.name not in user_inputs:
                user_inputs.append(user_input.name)
    return user_inputs
                

def _break_into_lines(one_long_line, width=72, indent1=0, indento=20):
    assert width > 0
    assert indent1 >= 0
    assert indento >= 0
    assert width > indent1
    assert width > indento

    assert "\n" not in one_long_line
    assert "\r" not in one_long_line
    assert "\t" not in one_long_line
    
    lines = []
    while 1:
        ind = " "*indent1
        if lines:
            ind = " "*indento
        if ind:
            one_long_line = one_long_line.lstrip()  # no leading spaces
        one_long_line = ind + one_long_line

        if len(one_long_line) < width:
            lines.append(one_long_line)
            break

        # Try to split on a space.
        w = width
        i = one_long_line.rfind(" ", len(ind), w)
        if i > 0:
            w = i
        x = one_long_line[:w]
        one_long_line = one_long_line[w:]
        lines.append(x)
    return lines
    

def list_datatypes(rulebase):
    # Make a list of the DataType objects.
    x = [getattr(rulebase, x) for x in dir(rulebase)]
    x = [x for x in x if isinstance(x, bie3.DataType)]
    datatypes = x

    # Make a list of the modules
    modules = rulebase.all_modules

    # Print each DataType object.
    for data in datatypes:
        print "DATATYPE %s:" % data.name
        for attr in data.attributes:
            x1 = "%-20s" % attr.name
            x2 = []
            for val in attr.values:
                if val == attr.default_in:
                    val = val + " (in)"
                if val == attr.default_out:
                    val = val + " (out)"
                x2.append(val)
            x2 = ", ".join(x2)
            x = x1 + x2
            lines = _break_into_lines(x)
            for line in lines:
                print line
        print

    # Print the user input from each module.
    for module in modules:
        if not module.user_inputs:
            continue
        print "MODULE %s:" % module.name
        for user_input in module.user_inputs:
            x1 = "%-20s" % user_input.name
            default = ""
            if user_input.default:
                default = user_input.default
            x2 = str(default)
            x = x1 + x2
            lines = _break_into_lines(x)
            for line in lines:
                print line
        print


def print_attribute(data):
    data = getattr(rulebase, data)
    if isinstance(data, bie3.DataType):
        attributes = data.attributes
        for attribute in attributes:
            print (data.name + '\t' + attribute.name + '\t' +
                   str(attribute.values) + '\t' + ' default_in\t' +
                   attribute.default_in + ' \t default_out\t' +
                   attribute.default_out)
        if attributes:
            print '-------------------------------'

def print_user_input(modules):
    print 'user_inputs:'
    user_inputs = {}
    for module in modules:
        for user_input in module.user_inputs:
            if user_input.name not in user_inputs:
                user_inputs[user_input.name] = user_input.default
                default = str(', default\t ' + str(user_input.default)
                              if user_input.default else '')
                print user_input.name + default


def assign_args(args):
    in_datatypes = []
    in_parameters = {}
    identifiers = []
    out_identifier = None
    outtype = None
    Attributes = []
    flag = None
    for i, arg in enumerate(args):
        if arg == "--intype":
            assert len(args) > i + 1
            in_datatypes.append(args[i + 1])
            flag = 'in'
        elif arg == "--outtype":
            assert len(args) > i + 1
            outtype = args[i+1]
            flag = 'out'
        elif arg == "--attr":
            assert len(args) > i + 1
            if flag == 'in':
                assert in_datatypes
                x = args[i + 1].split('=')
                key, value = x
                index = len(in_datatypes) - 1
                if index not in in_parameters:
                    in_parameters[index] = []
                in_parameters[index].append([key,value])
            elif flag == 'out':
                sub_datatype = args[i + 1].split(',')[0]
                attr = ','.join(args[i + 1].split(',')[1:])
                x = attr.split('=')
                key = x[0]
                value = x[1]
                Attributes.append([sub_datatype,key,value])
        elif arg == '--input':
            if flag == 'in':
                if not len(in_datatypes) == len(identifiers) + 1:
                    identifiers.extend([None] * (len(in_datatypes) - len(identifiers) - 1))
                identifiers.append(args[i + 1])  
            elif flag == 'out':
                out_identifier = args[i + 1]
    result = []
    for i, in_datatype in enumerate(in_datatypes):
        param = []
        if in_parameters:
            if i in in_parameters:
                param = in_parameters[i]
        result.append([in_datatype,identifiers[i],param])
    out = [outtype, out_identifier, Attributes]
    return result, out


def main():
    parser = argparse.ArgumentParser(description='Run the engine')
    group = parser.add_argument_group(title="Input/Output Nodes")
    group.add_argument(
        '--intype',  default=None, action='append',
        type=str, help='input data type for network')
    group.add_argument(
        '--input', default=None, action='append',
        type=str, help='input data path')
    group.add_argument(
        '--user_input',  default=[], action='append',
        type=str, help='user input in key=value format')
    group.add_argument(
        '--outtype',  default=None, type=str,
        help='outtype')
    group.add_argument(
        '--attr', default=[], type=str, action='append',
        help='attribute for datatype in "datatype,key=value" or "key=value" '
        'format')
    group.add_argument(
        '--output', type=str, default=None,
        help='file or folder of output result')
    
    group = parser.add_argument_group(title="Outfiles")
    group.add_argument(
        '--png_file',  type=str, default=None,
        help='generate the output network png file')
    group.add_argument(
        '--text_file', type=str, default=None,
        help='generate the output network text file')
    group.add_argument(
        '--json_file', type=str, default=None,
        help='generate the output network json file')
    
    group = parser.add_argument_group(title="Introspection")
    group.add_argument(
        '--list_datatypes', action='store_const',
        const=True, default=False,
        help='show all the possbile datatypes and their attributes')
    group.add_argument(
        '--list_attributes_for_network', 
        action='store_const', const=True, default=False,
        help='show all the possbile datatypes and attributes for the network')
    group.add_argument(
        '--diagnose', action='store_const',
        const=True, default=False,
        help='diagnose the input data')
    
    parser.add_argument(
        '--clobber', action='store_const',
        const=True, default=False,
        help='overwrite the output_data if it already exists')
    parser.add_argument(
        '--dry_run',  action='store_const',
        const=True, default=False,
        help='generate the network, do not run the network')
    # parse
    args = parser.parse_args()
    parser_list, out_list = assign_args(sys.argv)
    outtype, out_identifier, out_attributes = out_list
    # test
    assert not out_identifier,' --input is not for outtype'
    for x in parser_list:
        intype, identifier, attributes = x
        if identifier:
            assert os.path.exists(identifier),'input %s does not exists' %identifier
    # test outtype and build Attributes
    if outtype:
        goal_datatype = getattr(rulebase, outtype)
        Attributes = []
        for x in out_attributes:
            subtype, key, value = x
            fn = getattr(rulebase, subtype)
            Attributes.append(bie3.Attribute(fn, key, value))
    # test intype attributes and build objects
    in_objects = []
    for x in parser_list:
         intype, identifier, attributes = x
         fn = getattr(rulebase, intype)
         in_data = fn.input()
         if attributes:
             for i in attributes:
                 key, value = i
                 parameters[i][key] = value
                 in_data = fn.input(**parameters[i])
         if identifier:
             in_object = rule_engine_bie3.DataObject(in_data,identifier)
         in_objects.append(in_object)
    # test user_input
    all_inputs = get_all_user_inputs()
    user_inputs = {}
    for i in args.user_input:
        assert '=' in i
        key, value = i.split('=')
        assert key in all_inputs,'user input %s is not valid' % i
        user_inputs[key] = value
    # test introspection
    if args.output:
        realpath = os.path.realpath(args.output)
        if os.path.exists(args.output):
            if not args.clobber:
                raise ValueError('the output path %s is already exisit,\
                                 please use --clobber option to overwrite'
                                 % args.output)
    if args.list_datatypes:
        assert not (args.list_attributes_for_network or args.diagnose)
        assert not parser_list, 'no intype should be given'
        assert not outtype, 'no outtype should be given'
    if args.list_attributes_for_network:
        assert not (args.list_datatypes or args.diagnose)
        assert outtype,'an outtype should be given'
    if args.diagnose:
        assert not (args.llist_datatypes or list_attributes_for_network)
        assert out_list,'an outtype should be given'
    if not args.list_datatypes:
        assert args.outtype, 'please specify the outtype'
    # Action
    if args.list_datatypes:
         list_datatypes(rulebase)
         return
    print 'Generating network...'
    network = bie3.backchain(rulebase.all_modules, goal_datatype, *Attributes)
    network = bie3.complete_network(network)
    network = bie3.optimize_network(network)
    assert network, ('No pipeline has been generated,\
                      please check your command.')
    if args.png_file:
        bie3.plot_network_gv(args.png_file, network)

    if args.text_file:
        handle = file(args.text_file, 'w')
        try:
            bie3.print_network(network, handle)
        finally:
            handle.close()
    if args.json_file:
        bie3.write_network(args.json_file, network)
    if args.list_attributes_for_network:
        assert network, 'no network generated'
        network_datas = [i for i in network.nodes
                         if isinstance(i, bie3.Data)]
        network_modules = [i for i in network.nodes
                           if isinstance(i, bie3.Module)]
        datatype_names = list(set([i.datatype.name for i in network_datas]))
        for datatype_name in datatype_names:
            print_attribute(datatype_name)
        print_user_input(network_modules)
    if args.diagnose:
        in_datas = [i.data for i in in_objects]
        bie3.diagnose_start_node(network, in_datas)
    if args.dry_run:
        return
    assert in_objects, 'please specify the intype'
    for i in in_objects:
        start_node = bie3._find_start_nodes(network,i.data)
        assert start_node, 'intype %s is not matched any node in the network' % i.data.datatype.name
        store_file = userfile.set(getpass.getuser(), i.identifier)
        i.identifier = store_file
    output_file = rule_engine_bie3.run_pipeline(
        network, in_objects, user_inputs)
    if args.output:
        if os.path.exists(args.output) and args.clobber:
            if os.path.isdir(args.output):
                shutil.rmtree(args.output)
        if os.path.isdir(output_file):
            shutil.copytree(output_file, realpath)
        else:
            shutil.copy(output_file, realpath)


if __name__ == '__main__':
    main()
