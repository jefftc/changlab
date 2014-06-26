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

def parse_args(argv_list):
    in_datatypes = []
    in_parameters = {}
    Attributes = []
    user_inputs = {}
    identifiers = []
    flag = None
    for i, arg in enumerate(argv_list):
        if arg == "--intype":
            assert len(argv_list) > i + 1
            in_datatypes.append(argv_list[i + 1])
            flag = 'in'
        elif arg == "--outtype":
            assert len(argv_list) > i + 1
            flag = 'out'
        elif arg == "--attr":
            assert len(argv_list) > i + 1
            if flag == 'in':
                assert in_datatypes
                x = argv_list[i + 1].split('=')
                key = x[0]
                value = x[1]
                index = len(in_datatypes) - 1
                if index not in in_parameters:
                    in_parameters[index] = {}
                in_parameters[index][key] = value
            elif flag == 'out':
                sub_datatype = argv_list[i + 1].split(',')[0]
                attr = ','.join(argv_list[i + 1].split(',')[1:])
                x = attr.split('=')
                key = x[0]
                value = x[1]
                fn = getattr(rulebase, sub_datatype)
                Attributes.append(bie3.Attribute(fn, key, value))
        elif arg == '--input':
            if not len(in_datatypes) == len(identifiers) + 1:
                identifiers.extend([''] * (len(in_datatypes) - len(identifiers) - 1))
            assert os.path.exists(argv_list[i + 1]),'input %s does not exists' %argv_list[i + 1]
            store_file = userfile.set(getpass.getuser(), argv_list[i + 1])
            identifiers.append(store_file)
        elif arg == '--user_input':
            x = argv_list[i + 1].split('=')
            key = x[0]
            value = x[1]
            user_inputs[key] = value
    in_objects = []
    for i, in_datatype in enumerate(in_datatypes):
        fn = getattr(rulebase, in_datatype)
        in_data = fn.input()
        if in_parameters:
            if i in in_parameters:
                in_data = fn.input(**in_parameters[i])
        in_object = rule_engine_bie3.DataObject(in_data)
        if identifiers:
            in_object = rule_engine_bie3.DataObject(in_data, identifiers[i])
        in_objects.append(in_object)
    return in_objects, user_inputs, Attributes


def main():
    parser = argparse.ArgumentParser(description='Run the engine')
    parser.add_argument(
        '--intype', dest='in_datatype', default=None, action='append',
        type=str, help='input data type for network')
    parser.add_argument(
        '--input', dest='identifier', default=None, action='append',
        type=str, help='input data path')
    parser.add_argument(
        '--user_input', dest='user_input', default=None, action='append',
        type=str, help='user input in key=value format')
    parser.add_argument(
        '--outtype', dest='out_datatype', default=None, type=str,
        help='outtype')
    parser.add_argument(
        '--attr', dest='param', default=[], type=str, action='append',
        help='attribute for datatype in "datatype,key=value" or "key=value" '
        'format')
    parser.add_argument(
        '--list_datatypes', dest='all_datatypes', action='store_const',
        const=True, default=False,
        help='show all the possbile datatypes and their attributes')
    parser.add_argument(
        '--list_attributes_for_network', dest='network_attributes',
        action='store_const', const=True, default=False,
        help='show all the possbile datatypes and attributes for the network')
    parser.add_argument(
        '--output', dest='output', type=str, default=None,
        help='file or folder of output result')
    parser.add_argument(
        '--clobber', dest='clobber', action='store_const',
        const=True, default=False,
        help='overwrite the output_data if it already exists')
    parser.add_argument(
        '--dry_run', dest='dry_run', action='store_const',
        const=True, default=False,
        help='generate the network, do not run the network')
    parser.add_argument(
        '--png_file', dest='png_file', type=str, default=None,
        help='generate the output network png file')
    parser.add_argument(
        '--text_file', dest='text_file', type=str, default=None,
        help='generate the output network text file')
    parser.add_argument(
        '--json_file', dest='json_file', type=str, default=None,
        help='generate the output network json file')
    parser.add_argument(
        '--diagnose', dest='diagnose', action='store_const',
        const=True, default=False,
        help='diagnose the input data')
    args = parser.parse_args()
    if args.all_datatypes:
        list_datatypes(rulebase)
    if args.output:
        realpath = os.path.realpath(args.output)
        if os.path.exists(args.output):
            if not args.clobber:
                raise ValueError('the output path %s is already exisit,\
                                 please use --clobber option to overwrite'
                                 % args.output)
    if not args.out_datatype:
        if args.all_datatypes:
            return
        else:
            assert args.out_datatype, 'please specify the outtype'
    goal_datatype = getattr(rulebase, args.out_datatype)
    x = parse_args(sys.argv)
    in_objects, user_inputs, Attributes = x
    all_inputs = get_all_user_inputs()
    for user_input in user_inputs:
        assert user_input in all_inputs,'user input %s is not valid' % user_input
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
    if args.network_attributes:
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
    assert args.in_datatype, 'please specify the in_datatype'
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
