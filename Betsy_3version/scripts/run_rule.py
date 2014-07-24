#! /usr/bin/env python
#run_rule.py
import sys
import os
import shutil
import getpass
import argparse
from Betsy import rule_engine_bie3
#from Betsy import rule_engine_bie3_paralleling as rule_engine_bie3
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


def pretty_print_datatype(datatype, handle=None):
    handle = handle or sys.stdout
    
    print >>handle, "DATATYPE %s:" % datatype.name
    for attr in datatype.attributes:
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
            print >>handle, line


def pretty_print_module(module, handle=None):
    handle = handle or sys.stdout

    print >>handle, "MODULE %s:" % module.name
    for user_input in module.user_inputs:
        x1 = "%-20s" % user_input.name
        default = ""
        if user_input.default:
            default = user_input.default
        x2 = str(default)
        x = x1 + x2
        lines = _break_into_lines(x)
        for line in lines:
            print >>handle, line
            
    

def list_datatypes(rulebase):
    # Make a list of the DataType objects.
    x = [getattr(rulebase, x) for x in dir(rulebase)]
    x = [x for x in x if isinstance(x, bie3.DataType)]
    datatypes = x

    # Make a list of the modules
    modules = rulebase.all_modules

    # Print each DataType object.
    for dt in datatypes:
        pretty_print_datatype(dt)
        print

    # Print the user input from each module.
    for module in modules:
        if not module.user_inputs:
            continue
        pretty_print_module(module)
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
def get_necessary_user_input(modules):
    user_inputs = []
    for module in modules:
        for user_input in module.user_inputs:
            if user_input.name not in user_inputs:
                if not user_input.default:
                    user_inputs.append(user_input.name)
    return user_inputs

def assign_args(args):
    intypes = []
    in_parameters = {}
    in_identifiers = []
    outtype = None
    out_parameters = []
    out_identifier = None
    flag = None
    for i, arg in enumerate(args):
        if arg == "--input":
            assert len(args) > i + 1
            intypes.append(args[i + 1])
            flag = 'in'
        elif arg == "--output":
            assert len(args) > i + 1
            outtype = args[i+1]
            flag = 'out'
        elif arg == "--dattr":
            assert len(args) > i + 1
            if flag == 'in':
                assert intypes
                x = args[i + 1].split('=')
                key, value = x
                index = len(intypes) - 1
                if index not in in_parameters:
                    in_parameters[index] = []
                in_parameters[index].append([key,value])
            elif flag == 'out':
                sub_datatype = args[i + 1].split(',')[0]
                attr = ','.join(args[i + 1].split(',')[1:])
                x = attr.split('=')
                key = x[0]
                value = x[1]
                out_parameters.append([sub_datatype,key,value])
        elif arg == '--input_file':
            if flag == 'in':
                if not len(intypes) == len(in_identifiers) + 1:
                    in_identifiers.extend([None] * (len(intypes) - len(in_identifiers) - 1))
                in_identifiers.append(args[i + 1])  
            elif flag == 'out':
                out_identifier = args[i + 1]
    in_result = []
    if not in_identifiers:
        in_identifiers = ['']*len(intypes)
    for i, intype in enumerate(intypes):
        param = []
        if in_parameters:
            if i in in_parameters:
                param = in_parameters[i]
        in_result.append([intype,in_identifiers[i],param])
    out_result = [outtype, out_identifier, out_parameters]
    return in_result, out_result


def main():
    parser = argparse.ArgumentParser(description='Run the engine')
    group = parser.add_argument_group(title="Input/Output Nodes")
    group.add_argument(
        '--input',  default=None, action='append',
        type=str, help='input data type for network')
    group.add_argument(
        '--input_file', default=None, action='append',
        type=str, help='input data path')
    group.add_argument(
        '--mattr',  default=[], action='append',
        type=str, help='module attribute in key=value format')
    group.add_argument(
        '--output',  default=None, type=str,
        help='output type')
    group.add_argument(
        '--dattr', default=[], type=str, action='append',
        help='attribute for datatype in "datatype,key=value" or "key=value" '
        'format')
    group.add_argument(
        '--output_file', type=str, default=None,
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
    if not args.output and not args.input:
        list_datatypes(rulebase)
        return
    if not args.output and args.input:
        raise ValueError("output is expected")
    assert not out_identifier,' --input_file is not for outtype'
    for x in parser_list:
        intype, identifier, attributes = x
        if identifier:
            assert os.path.exists(identifier),'input_file %s does not exists' %identifier
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
             parameters = {}
             for i in attributes:
                 key, value = i
                 parameters[key] = value
                 in_data = fn.input(**parameters)
         in_object = rule_engine_bie3.DataObject(in_data,identifier)
         in_objects.append(in_object)
    # test mattr are valid
    all_inputs = get_all_user_inputs()
    user_inputs = {}
    for i in args.mattr:
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
        
    if args.output and not in_objects:
        assert network, 'no network generated'
        print "Possible Inputs"
        inputs = bie3.get_inputs(network)
        dt2inputs = bie3.group_inputs_by_datatype(network, inputs)
        for i, dt in enumerate(sorted(dt2inputs)):
            x = [x.name for x in dt]
            print "%d.  %s" % (i+1, ", ".join(x))
            for j, inputs in enumerate(dt2inputs[dt]):
                for k, inp in enumerate(inputs):
                    node = network.nodes[inp]
                    assert isinstance(node, bie3.Data)
                    print node.datatype.name
                    for name in sorted(node.attributes):
                        print "%s%s=%s" % (" "*5, name, node.attributes[name])
                print
            print
        return 
    if args.dry_run:
        return
    if not network or len(network.nodes)==1:
        in_datas = [i.data for i in in_objects]
        bie3.diagnose_start_node(network, in_datas)
    for i in in_objects:
        start_node = bie3._find_start_nodes(network,i.data)
        assert start_node, 'input %s is not matched any node in the network' % i.data.datatype.name
        if os.path.exists(i.identifier):
            store_file = userfile.set(getpass.getuser(), i.identifier)
            i.identifier = store_file
    #test mattr are given when necessary
    network_modules = [i for i in network.nodes
                           if isinstance(i, bie3.Module)]
    necessary_user_inputs = get_necessary_user_input(network_modules)
    for user_input in necessary_user_inputs:
        assert user_input in user_inputs, 'mattr %s should be given' % user_input
        
    output_file = rule_engine_bie3.run_pipeline(
        network, in_objects, user_inputs)
    
    if args.output_file:
        if os.path.exists(args.output_file) and args.clobber:
            if os.path.isdir(args.output_file):
                shutil.rmtree(args.output_file)
        if os.path.isdir(output_file):
            shutil.copytree(output_file, realpath)
        else:
            shutil.copy(output_file, realpath)


if __name__ == '__main__':
    main()

