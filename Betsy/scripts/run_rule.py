#! /usr/bin/env python


def get_all_option_names():
    # Return a list of the names for all OptionDef objects in the
    # rulebase.
    from Betsy import rulebase
    
    names = {}
    for module in rulebase.all_modules:
        for option in module.option_defs:
            names[option.name] = 1
    return sorted(names)


def get_required_option_names(modules):
    # From a list of modules, return a dictionary where the keys are
    # option name and the values are a list of the modules that
    # require it.

    option2modules = {}
    for module in modules:
        for option in module.option_defs:
            if option.default is not None:
                continue
            on, mn = option.name, module.name
            if on not in option2modules:
                option2modules[on] = []
            option2modules[on].append(mn)
    return option2modules
                

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
    import sys
    
    handle = handle or sys.stdout

    LW = 72
    print >>handle, "="*LW
    print >>handle, "DataType: %s" % datatype.name
    if datatype.help:
        for x in _break_into_lines(datatype.help, indent1=10, indento=10):
            print >>handle, x
    print >>handle, "-"*LW
    for attr in datatype.attribute_defs.itervalues():
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
    import sys
    
    handle = handle or sys.stdout

    LW = 72
    print >>handle, "="*LW
    print >>handle, "ModuleNode: %s" % module.name
    if module.help:
        for x in _break_into_lines(module.help, indent1=8, indento=8):
            print >>handle, x
    print >>handle, "-"*LW
    for option in module.option_defs:
        x1 = "%-20s" % option.name
        x2 = ""
        if option.help:
            x2 = "%s " % option.help
        default = ""
        if option.default:
            default = "(default %s)" % option.default
        x3 = str(default)
        x = x1 + x2 + x3
        x = x.strip()
        lines = _break_into_lines(x)
        for line in lines:
            print >>handle, line
            
    

def list_datatypes(rulebase):
    from Betsy import bie3 

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
        print

    # Print the options from each module.
    for module in modules:
        if not module.option_defs:
            continue
        pretty_print_module(module)
        print
        print


def print_attribute(data):
    from Betsy import rulebase
    from Betsy import bie3 

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

def print_option(modules):
    print 'options:'
    all_options = []
    for module in modules:
        for option in module.option_defs:
            all_options.append(option)
            
    seen = {}
    for option in all_options:
        if option.name in seen:
            continue
        seen[option.name] = 1
        
        default = str(', default\t ' + str(option.default)
                      if option.default else '')
        print option.name + default
        print option.help
                

def assign_args(args):
    inputs = []            # list of names of datatypes
    in_parameters = {}     # index (into inputs) -> list of (key, value)
    in_identifiers = {}    # index (into inputs) -> identifier (e.g. filename)
    
    output = None          # name of datatype
    out_parameters = []    # list of (key, value)
    out_identifier = None  # identifier

    input_or_output = None
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == "--input":
            assert len(args) >= i+1
            inputs.append(args[i + 1])
            input_or_output = arg
            i += 2
        elif arg == "--output":
            assert len(args) >= i+1
            assert not output, "multiple --output not allowed"
            output = args[i+1]
            input_or_output = arg
            i += 2
        elif arg == "--dattr" and input_or_output == "--input":
            assert len(args) >= i+1
            dattr = args[i+1]
            i += 2

            # Format: <datatype>,<key>=<value>
            assert inputs
            x = dattr.split(",", 1)
            assert len(x) == 2, \
                   "--dattr should be <datatype>,<key>=<value>"
            datatype, x = x
            x = x.split("=", 1)
            assert len(x) == 2, \
                   "--dattr should be <datatype>,<key>=<value>"
            key, value = x
            index = len(inputs) - 1
            assert datatype == inputs[index], (
                "Datatype mismatch: --dattr %s and --input %s.  Wrong "
                "order?" % (datatype, inputs[index]))
            if index not in in_parameters:
                in_parameters[index] = []
            in_parameters[index].append((key,value))
        elif arg == "--dattr" and input_or_output == "--output":
            assert len(args) >= i+1
            dattr = args[i+1]
            i += 2
            
            # Format: <datatype>,<key>=<value>
            x = dattr.split(",", 1)
            assert len(x) == 2, \
                   "--dattr should be <datatype>,<key>=<value>"
            datatype, x = x
            x = x.split('=', 1)
            assert len(x) == 2, \
                   "--dattr should be <datatype>,<key>=<value>"
            key, value = x
            out_parameters.append((datatype, key, value))
        elif arg == "--dattr":
            raise AssertionError, "--dattr before --input or --output"
        elif arg == '--input_file':
            assert input_or_output == "--input"
            assert len(args) >= i+1
            filename = args[i+1]
            i += 2
            index = len(inputs) - 1
            assert index not in in_identifiers, \
                   "only one --input_file per --input"
            in_identifiers[index] = filename
        elif arg == '--output_file':
            assert input_or_output == "--output", \
                   "--output must precede --output_file"
            assert len(args) >= i+1
            filename = args[i+1]
            i += 2
            assert not out_identifier
            out_identifier = filename
        else:
            i += 1

    in_results = []
    for i, input_ in enumerate(inputs):
        x = input_, in_identifiers.get(i, ""), in_parameters.get(i, [])
        in_results.append(x)
    out_result = output, out_identifier, out_parameters
    return in_results, out_result

def print_possible_inputs(network, user_attributes):
    from Betsy import bie3 

    print "Possible Inputs"
    inputs = bie3.get_inputs(network, user_attributes)
    dt2inputs = bie3.group_inputs_by_datatype(network, inputs)
    for i, dt in enumerate(sorted(dt2inputs)):
        x = [x.name for x in dt]
        print "%d.  %s" % (i+1, ", ".join(x))
        for j, inputs in enumerate(dt2inputs[dt]):
            for k, inp in enumerate(inputs):
                node = network.nodes[inp]
                assert isinstance(node, bie3.DataNode)
                print node.datatype.name
                print node.datatype.help
                for name in sorted(node.attributes):
                    print "%s%s=%s" % (" "*5, name, node.attributes[name])
            print
        print
        
def is_complete_input(network, user_attributes, start_node_ids):
    # Return a boolean indicating whether this is a valid set of start
    # nodes.
    from Betsy import bie3 

    start_node_ids = sorted(start_node_ids)
    inputs = bie3.get_inputs(network, user_attributes)
    for input_ in inputs:
        # Bug: should check if input_ is subset of start_node_ids.
        if start_node_ids == sorted(input_):
            return True
    return False
    #dt2inputs = bie3.group_inputs_by_datatype(network, inputs)
    #for i, dt in enumerate(sorted(dt2inputs)):
    #    x = [x.name for x in dt]
    #    for j, dtinput in enumerate(dt2inputs[dt]):
    #        required = sorted(list(set(dtinput)))
    #        if start_node_ids == required:
    #            return True
    #return False

## def print_missing_inputs(network, user_attributes, input_nodes):
##     """print the missing inputs which is required to generate a sub network"""
##     from Betsy import bie3 

##     input_nodes = sorted(input_nodes)
##     inputs = bie3.get_inputs(network, user_attributes)
##     dt2inputs = bie3.group_inputs_by_datatype(network, inputs)
##     all_possible_required = []
##     for i, dt in enumerate(sorted(dt2inputs)):
##         x = [x.name for x in dt]
##         for j, dtinput in enumerate(dt2inputs[dt]):
##             required = sorted(list(set(dtinput)))
##             if required not in all_possible_required:
##                 all_possible_required.append(required)
##     print 'Please provide the following input to generate the network'
##     for required in all_possible_required:
##         if set(input_nodes).issubset(required):
##             missing_nodes = list(set(required)-set(input_nodes))
##             for inp in missing_nodes:
##                 node = network.nodes[inp]
##                 assert isinstance(node, bie3.DataNode)
##                 print node.datatype.name
##                 print node.datatype.help
##                 for name in sorted(node.attributes):
##                     print "%s%s=%s" % (" "*5, name, node.attributes[name])
##             print     


def main():
    import os
    import sys
    import argparse
    import shutil
    
    import getpass
    from Betsy import config
    from Betsy import module_utils
    from Betsy import rule_engine_bie3
    from Betsy import userfile
    from Betsy import rulebase
    from Betsy import bie3 
    
    parser = argparse.ArgumentParser(description='Run the engine')
    parser.add_argument(
        "--run",  action="store_const", const=True,
        help="Perform the analysis.")
    parser.add_argument(
        '--user', default=getpass.getuser(),
        help='the username who run the command')
    parser.add_argument(
        '--job_name', default='', help='the name of this job')
    parser.add_argument(
        '--num_cores', default=8, type=int,
        help='number of cores used in the processing')
    #parser.add_argument(
    #    "--diagnose", action="store_true",
    #    help="Try to figure out why a start node can't be found.")
    
    group = parser.add_argument_group(title="Input/Output Nodes")
    group.add_argument(
        '--input',  action='append', help='DataType of the input')
    group.add_argument(
        '--input_file', action='append',
        help='File corresponding to the previous --input.')
    group.add_argument(
        '--mattr', default=[], action='append',
        help='Set the option for a module.  Format should be: '
        '<key>=<value>.')
    
    group.add_argument(
        '--output',  help='Desired DataType for the output.')
    group.add_argument(
        '--dattr', default=[], action='append',
        help="Specify a Datatype's attribute.  Since there may be multiple "
        "inputs or outputs with the same Datatype, the position of "
        "this argument is important.  --dattr should be given "
        "immediately following the Datatype that it refers to.  "
        "No --dattr should be given before the first --input.  "
        "--dattr that refers to nodes not in the network should be "
        "given after the --output.  Format: <datatype>,<key>=<value>.")
    group.add_argument(
        '--output_file', help='file or folder of output result')
    
    group = parser.add_argument_group(title="Outfiles")
    group.add_argument(
        '--network_png', help='generate the output network png file')
    group.add_argument(
        '--unpruned_network_png', help='generate the output network png file')
    group.add_argument(
        '--sparse_network_png', action="store_true",
        help="Leave out details in network plot.")
    group.add_argument(
        '--unpruned_network_text', help='')
    group.add_argument(
        '--network_text', help='generate the output network text file')
    group.add_argument(
        '--network_json', help='generate the output network json file')
    #parser.add_argument(
    #    '--clobber', action='store_const', const=True, default=False,
    #    help='overwrite the output_data if it already exists')
    #parser.add_argument(
    #    '--dry_run',  action='store_const', const=True, default=False,
    #    help='generate the network, do not run the network')

    print "Starting rule engine."
    sys.stdout.flush()
    
    # Parse the arguments.
    args = parser.parse_args()
    input_list, x = assign_args(sys.argv)
    outtype, out_identifier, out_attributes = x

    args.clobber = True
    assert args.num_cores > 0,'num_cores should be greater than 0'
    
    # Make sure configuration directory exists.
    if not os.path.exists(config.OUTPUTPATH):
        print "Making BETSY working path: %s." % config.OUTPUTPATH
        os.mkdir(config.OUTPUTPATH)

    # Make sure files exist.
    if input_list or args.output_file:
        print "Looking for files."
    for x in input_list:
        intype, identifier, attributes = x
        if identifier:
            assert os.path.exists(identifier), \
                   "File not found: %s" % identifier
    if args.output_file:
        args.output_file = os.path.realpath(args.output_file)
        if os.path.exists(args.output_file):
            assert args.clobber, "Output already exists: %s" % args.output_file

    # Making the IdentifiedDataNode objects.
    if input_list:
        print "Parsing input options."
    in_data_nodes = []
    for x in input_list:
        intype, identifier, attributes = x
        assert hasattr(rulebase, intype), "Unknown datatype: %s" % intype
        fn = getattr(rulebase, intype)
        params = {}
        for x in attributes:
            key, value = x
            params[key] = value
        in_data = fn.input(**params)
        x = module_utils.IdentifiedDataNode(in_data, identifier)
        in_data_nodes.append(x)
         
    # test outtype and build the list of user_attributes.
    user_attributes = []  # List of bie3.Attribute objects.
    if outtype:
        print "Parsing output options."
        assert hasattr(rulebase, outtype), "Unknown datatype: %s" % outtype
        out_datatype = getattr(rulebase, outtype)
        for x in out_attributes:
            subtype, key, value = x
            assert hasattr(rulebase, subtype), "Unknown datatype: %s" % subtype
            fn = getattr(rulebase, subtype)
            user_attributes.append(bie3.Attribute(fn, key, value))

    #fn = getattr(rulebase, "SignalFile")
    #x = bie3.Attribute(fn, "preprocess", "rma")
    #user_attributes.append(x)
    #fn = getattr(rulebase, "SignalFile")
    #x = bie3.Attribute(fn, "quantile_norm", "yes")
    #user_attributes.append(x)

    # Store files in the user's directory.
    for i_data_node in in_data_nodes:
        filename = i_data_node.identifier
        if not os.path.exists(filename):
            continue
        x = userfile.set(getpass.getuser(), filename)
        i_data_node.identifier = x
            
    # test mattr are valid
    print "Parsing module attributes."
    all_mattrs = get_all_option_names()
    user_options = {}
    for x in args.mattr:
        assert '=' in x, "--mattr should be in format: <option>=<value>"
        key, value = x.split('=', 1)
        assert key in all_mattrs, "Unknown module attribute: %s" % key
        user_options[key] = value

    # Possible cases.
    # 1.  No inputs and no outputs.
    #     Show the whole rulebase.
    # 2.  No inputs and has output.  (requires network)
    #     Show the possible inputs.
    # 3.  Input and no outputs.
    #     ERROR
    # 4.  Input and output.  (requires network)
    #     Can run the network.
            
    # Case 1.  No inputs and no outputs.
    if not args.output and not args.input:
        # If there's no inputs and no outputs, just print the whole
        # rulebase.
        list_datatypes(rulebase)
        return
    # Case 3.  Input and no outputs.
    if args.input and not args.output:
        raise AssertionError, "Missing --output"

    print "Generating a network that produces a %s." % outtype
    sys.stdout.flush()
    network = bie3.backchain(
        rulebase.all_modules, out_datatype, user_attributes)
    network = bie3.complete_network(network, user_attributes)
    network = bie3.optimize_network(network, user_attributes)
    # When does this happen?  What can go wrong with the command?
    assert network, 'No network generated.  Please check your command.'

    if args.unpruned_network_png:
        print "Plotting unpruned network: %s." % args.unpruned_network_png
        bie3.plot_network_gv(
            args.unpruned_network_png, network, options=user_options,
            verbose=(not args.sparse_network_png))
    if args.unpruned_network_text:
        print "Writing detailed network: %s." % args.unpruned_network_text
        bie3.print_network(network, outhandle=args.unpruned_network_text)

    # Case 2.  No input and has output.
    if args.output and not in_data_nodes:
        # If there's an output and no inputs, then show the possible
        # inputs.
        print "No inputs given.  Here are the possibilities."
        print
        print_possible_inputs(network, user_attributes)
        return 

    # Case 4.  Has input and output.
    assert args.output and in_data_nodes
    
    # See if I have a complete set of inputs.
    # If this is taking too long, can speed it up by simplifying the
    # network with more constraints.
    print "Searching for each input in the network."
    sys.stdout.flush()
    start_node_ids = []
    missing = []  # list of indexes into in_data_nodes that are not found.
    for i, node in enumerate(in_data_nodes):
        x = bie3._find_start_nodes(network, node.data)
        start_node_ids.extend(x)
        if not x:
            missing.append(i)

    # Re-plot the network now that we know which start nodes are matched.
    if start_node_ids and args.unpruned_network_png:
        print "Plotting unpruned network: %s." % args.unpruned_network_png
        bie3.plot_network_gv(
            args.unpruned_network_png, network, options=user_options,
            highlight_node_ids=start_node_ids,
            verbose=(not args.sparse_network_png))
        
    # If any node_ids are missing, try to diagnose it.
    if missing:
        x = [in_data_nodes[i] for i in missing]
        x = [x.data.datatype.name for x in x]
        node_or_nodes = "node"
        if len(x) > 1:
            node_or_nodes = "nodes"
        print "Input %s not found: %s." % (node_or_nodes, ", ".join(x))
        print ("It is likely there is an incompatibility in the attributes "
               "somewhere.")
        print
        bie3.diagnose_start_node(network, in_datas)
        return
        
    if 0:  # This takes way too long.
        if not is_complete_input(network, user_attributes, start_node_ids):
            print_possible_inputs(network, user_attributes)
            return

    ## If the network has only one node, try to figure out why there's
    ## no pipeline.  Probably means incompatibility with start node
    ## somewhere.
    #if args.diagnose or len(network.nodes) == 1:
    #    in_datas = [x.data for x in in_data_nodes]
    #    bie3.diagnose_start_node(network, in_datas)
    #    return

    # Prune network.
    x = [x.data for x in in_data_nodes]
    x = [x.datatype.name for x in x]
    x = sorted({}.fromkeys(x))
    x_str = ", ".join(x)
    input_or_inputs = "input"
    if len(x) > 1:
        input_or_inputs = "inputs"
    print "Pruning the network for paths that start from %s: %s." % (
        input_or_inputs, x_str)
    x = [x.data for x in in_data_nodes]
    network = bie3.select_start_node(network, x)
    assert network.nodes, \
           "No network.  Most likely not enough inputs to generate result."
   
    if args.network_png:
        print "Plotting network: %s." % args.network_png
        bie3.plot_network_gv(
            args.network_png, network, options=user_options)
    if args.network_text:
        print "Writing detailed network: %s." % args.network_text
        bie3.print_network(network, outhandle=args.network_text)
    if args.network_json:
        print "Writing network in json format: %s." % args.network_json
        bie3.write_network(args.network_json, network)
        
    # Make sure all required mattr are given.  This can only be run
    # after the final network is generated.
    # TODO: "option" should be "module attribute"
    modules = [
        x for x in network.nodes if isinstance(x, bie3.ModuleNode)]
    opt2mods = get_required_option_names(modules)
    x = opt2mods.keys()
    x = [x for x in x if x not in user_options]
    for on in x:
        for mn in opt2mods[on]:
            print 'Missing --mattr: %s requires attribute "%s".' % (mn, on)
    if x:
        return
     
    if not args.run:
        print "Completed network operations.  Use --run to run the analysis."
        return
    print "Running the analysis."
    output_file = rule_engine_bie3.run_pipeline(
        network, in_data_nodes, user_attributes, user_options, args.user,
        args.job_name, args.num_cores)
    if not output_file:
        return
    if not args.output_file:
        return
    
    # Remove the previous files if they exist.
    if args.clobber and os.path.exists(args.output_file):
        if os.path.isdir(args.output_file):
            shutil.rmtree(args.output_file)
        else:
            os.unlink(args.output_file)
    # Copy the file or directory.
    if os.path.isdir(output_file):
        shutil.copytree(output_file, args.output_file)
    else:
        shutil.copy(output_file, args.output_file)


if __name__ == '__main__':
    main()
