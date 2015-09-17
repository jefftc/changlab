#! /usr/bin/env python

# Functions:
# print_rulebase            To help select output DataType.
# _pretty_print_datatype
# _pretty_print_module
# print_input_datatypes     To help select input DataType(s).
# print_input_nodes         To help configure input attributes.
# print_diagnose            Ready to run, but something is wrong with input.
# plot_network
#
# is_complete_input           SLOW
# get_all_option_names        Get the names of all OptionDef in the rulebase.
# get_required_option_names   Get the required options in a list of modules.
#
# _prune_unwanted_input_nodes
# _list_differences_in_nodelists
# _merge_nodelists
#
# _parse_args                 Parse the arguments from the user.


# TODO: Rename this script.
# make_it_so_betsy.py
# analyze_it.py
# pipeline.py
# betsy_help_me.py


def print_rulebase(rulebase):
    from Betsy import bie3

    # Make a list of the DataType objects.
    x = [getattr(rulebase, x) for x in dir(rulebase)]
    x = [x for x in x if isinstance(x, bie3.DataType)]
    datatypes = x

    # Make a list of the modules
    modules = rulebase.all_modules

    # Print each DataType object.
    for dt in datatypes:
        # Skip the private datatypes.
        if dt.name.startswith("_"):
            continue
        _pretty_print_datatype(dt)
        print
        print

    # Print the options from each module.
    for module in modules:
        if not module.option_defs:
            continue
        _pretty_print_module(module)
        print
        print


def _pretty_print_datatype(datatype, handle=None):
    import sys
    from genomicode import parselib

    handle = handle or sys.stdout

    LW = 72
    print >>handle, "="*LW
    print >>handle, "DataType: %s" % datatype.name
    if datatype.help:
        for x in parselib.linesplit(datatype.help, prefix1=10, prefixn=10):
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
        lines = parselib.linesplit(x, prefixn=20)
        for line in lines:
            print >>handle, line


def _pretty_print_module(module, handle=None):
    import sys
    from genomicode import parselib

    handle = handle or sys.stdout

    LW = 72
    print >>handle, "="*LW
    print >>handle, "Module: %s" % module.name
    if module.help:
        for x in parselib.linesplit(module.help, prefix1=8, prefixn=8):
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
        lines = parselib.linesplit(x)
        for line in lines:
            print >>handle, line


def print_input_datatypes(
    network, required_datatypes, excluded_datatypes, user_attributes,
    max_inputs=None, outhandle=None):
    # required_datatypes is a list of the names of datatypes that must
    # be in this combination.
    import sys
    from genomicode import parselib
    from Betsy import bie3

    outhandle = outhandle or sys.stdout

    ps = parselib.print_split

    x1 = ""
    if max_inputs is not None:
        x1 += " (up to %d)" % max_inputs
    x2 = "Acceptable combinations of DataTypes%s.  " % x1
    x3 = "The exact number of each DataType needed is not shown."
    x = x2 + x3
    ps(x, prefixn=0, outhandle=outhandle)
    outhandle.flush()

    datatype_combos = bie3.get_input_datatypes(
        network, user_attributes, skip_datatypes=excluded_datatypes,
        skip_private_datatypes=True, max_inputs=max_inputs)

    # Filter out the datatype_combos that do not contain all the
    # datatypes in required_datatypes.
    i = 0
    while i < len(datatype_combos):
        combo = datatype_combos[i]
        names = [x.name for x in combo]
        missing = False
        for n in required_datatypes:
            if n not in names:
                missing = True
                break
        if missing:
            del datatype_combos[i]
        else:
            i += 1
    # Sort by number of datatypes, then names.
    schwartz = []
    for combo in datatype_combos:
        x1 = len(combo)
        x2 = [x.name for x in combo]
        schwartz.append((x1, x2, combo))
    schwartz.sort()
    datatype_combos = [x[-1] for x in schwartz]

    if not datatype_combos:
        print >>outhandle, \
              "No combination of DataTypes can generate this network."
        return
    for i, combo in enumerate(datatype_combos):
        x = [x.name for x in combo]
        ps("%3d.  %s" % (i+1, ", ".join(x)), prefixn=6, outhandle=outhandle)
    print >>outhandle


def print_input_nodes(
    network, required_datatypes, excluded_datatypes, user_attributes,
    max_inputs=None, outhandle=None):
    # required_datatypes is a list of the names of datatypes that must
    # be in this combination.
    import sys
    from genomicode import parselib
    from Betsy import bie3

    outhandle = outhandle or sys.stdout

    ps = parselib.print_split

    print >>outhandle, "Possible Inputs"
    outhandle.flush()

    inputs = bie3.get_input_nodes(
        network, user_attributes, skip_datatypes=excluded_datatypes,
        skip_private_datatypes=True, max_inputs=max_inputs)
    dt2inputs = bie3.group_nodes_by_datatype(network, inputs)

    # Make a list of the datatypes to show.
    datatype_combos = dt2inputs.keys()  # list of datatypes

    ## # Remove any that contain "private" datatypes.
    ## i = 0
    ## while i < len(datatype_combos):
    ##     has_private = False
    ##     for x in datatype_combos[i]:
    ##         if x.name.startswith("_"):
    ##             has_private = True
    ##             break
    ##     if has_private:
    ##         del datatype_combos[i]
    ##     else:
    ##         i += 1


    # Remove any that do not contain all of the required_datatypes.
    # Also remove those that contain datatypes not in
    # required_datatypes.
    s_required_datatypes = sorted({}.fromkeys(required_datatypes))
    i = 0
    while i < len(datatype_combos):
        combo = datatype_combos[i]
        names = [x.name for x in combo]
        x = sorted({}.fromkeys(names))
        if x != s_required_datatypes:
            del datatype_combos[i]
        else:
            i += 1
    # Sort by number of datatypes, then names.
    schwartz = []
    for dt in datatype_combos:
        x1 = len(dt)
        x2 = [x.name for x in dt]
        schwartz.append((x1, x2, dt))
    schwartz.sort()
    datatype_combos = [x[-1] for x in schwartz]

    if not datatype_combos:
        print >>outhandle, "None"
        return
    for i, combo in enumerate(datatype_combos):
        x = [x.name for x in combo]
        ps("%d.  INPUTS: %s" % (i+1, ", ".join(x)), outhandle=outhandle)

        # node_id_combos a list of list of node_ids.  node_combos is
        # a list of list of nodes, all with the same combination of
        # types.
        # E.g. [
        #   [AgilentFiles, ClassLabelFile],
        #   [AgilentFiles, ClassLabelFile],
        # ]
        node_id_combos = dt2inputs[combo]
        node_combos = [None] * len(node_id_combos)
        for i, id_combo in enumerate(node_id_combos):
            # id_combo is a list of node_ids.  node_combo is a list of
            # DataNodes.
            node_combo = [network.nodes[x] for x in id_combo]
            # Sort the DataNodes by datatype name, then attributes.
            schwartz = []
            for id_, n in zip(id_combo, node_combo):
                x1 = n.datatype.name
                x2 = n.attributes
                schwartz.append((x1, x2, id_, n))
            schwartz.sort()
            node_id_combos[i] = [x[-2] for x in schwartz]
            node_combos[i] = [x[-1] for x in schwartz]

        # DEBUG: Make sure no IDs occur more than once
        #for x in node_id_combos:
        #    y = {}.fromkeys(x).keys()
        #    assert len(x) == len(y)

        # Merge the node_combos with very similar attributes.  Do this
        # with a greedy algorithm.  Take each node_combo, and make a
        # group of other node_combos to merge.  Every member of the
        # node_combo must vary by at most one attribute (and must be
        # the same attribute).
        j = 0
        while j < len(node_combos)-1:
            # Indexes of node_combos to merge, not including this
            # combo.
            group = []

            # Find other nodes that fit into this group.
            diff = None
            for k in range(j+1, len(node_combos)):
                combo1, combo2 = node_combos[j], node_combos[k]
                x = _list_differences_in_nodelists(combo1, combo2)
                if len(x) > 1:
                    continue
                if not x:
                    # Sometimes can have no differences.  Just merge
                    # them.
                    pass
                elif diff is None:
                    diff = x
                elif x != diff:
                    continue
                group.append(k)
            if not group:
                j += 1
                continue
            # Merge every combo in this group, and delete the unmerged
            # combinations.
            node_combo = node_combos[j]
            for k in group:
                node_combo = _merge_nodelists(node_combo, node_combos[k])
            node_combos[j] = node_combo
            node_combos = [
                node_combos[k] for k in range(len(node_combos))
                if k not in group]

        # Print out each combination.
        for j, node_combo in enumerate(node_combos):
            # Print one node.
            for k, node in enumerate(node_combo):
                assert isinstance(node, bie3.DataNode)
                #x = "- %s: %s" % (node.datatype.name, node.datatype.help)
                x = "- %s" % node.datatype.name
                ps(x, prefix1=4, prefixn=8, outhandle=outhandle)

                # Print the attributes.
                for name in sorted(node.attributes):
                    value = node.attributes[name]
                    # To save space, don't print the ones with default
                    # values.
                    attrdef = node.datatype.attribute_defs[name]
                    if attrdef.default_in == value:
                        continue
                    #if type(value) is type([]) and \
                    #       attrdef.default_in in value:
                    #    continue
                    x = "%s=%s" % (name, value)
                    ps(x, prefix1=8, prefixn=10, outhandle=outhandle)
            print >>outhandle  # space between each combination.
            outhandle.flush()
        print >>outhandle      # two spaces at end of this set of combos.


def print_diagnose(network, start_nodes, missing, outhandle=None):
    # start_nodes is a list of DataNode objects.  missing is a list of
    # the indexes of start_nodes that cannot be found in the network.
    from genomicode import parselib
    from Betsy import bie3

    x = [start_nodes[i] for i in missing]
    x = [str(x) for x in x]
    node_or_nodes = "node"
    if len(x) > 1:
        node_or_nodes = "nodes"
    x = "Input %s not found:\n%s" % (node_or_nodes, "\n".join(x))
    parselib.print_split(x, outhandle=outhandle)
    parselib.print_split(
        "Possible reasons are:", outhandle=outhandle)
    parselib.print_split(
        "- There is an incompatibility in the attributes somewhere.",
        prefix1=4, outhandle=outhandle)
    parselib.print_split(
        "- It is not needed for this pipeline.", prefix1=4,
        outhandle=outhandle)
    print >>outhandle
    # TODO: Move code to this module.
    bie3.diagnose_start_node(network, start_nodes, outhandle=outhandle)


def plot_network(filename, network, user_options=None, highlight_node_ids=None,
                 verbose=False):
    from Betsy import bie3

    highlight_node_ids = highlight_node_ids or []

    if filename is None:
        return
    print "Plotting network to %s." % filename
    bie3.plot_network_gv(
        filename, network, options=user_options,
        highlight_node_ids=highlight_node_ids, verbose=verbose)


def is_complete_input(network, user_attributes, start_node_ids):
    # Return a boolean indicating whether this is a valid set of start
    # nodes.  (SLOW).
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


def _delete_input_nodes(
    network, wanted_dt_names, unwanted_dt_names, user_attributes):
    # Only one of wanted_dt_names and unwanted_dt_names should be given.
    from Betsy import bie3

    assert not (wanted_dt_names and unwanted_dt_names)

    wanted_dt_names = wanted_dt_names or []
    unwanted_dt_names = unwanted_dt_names or []
    wanted_dict = {}.fromkeys(wanted_dt_names)
    unwanted_dict = {}.fromkeys(unwanted_dt_names)

    while True:
        # Keep all the internal IDs.  These are the nodes that someone
        # points to.
        keep_ids = {}
        for next_ids in network.transitions.values():
            x = {}.fromkeys(next_ids)
            keep_ids.update(x)
        # Everything else is a start ID.
        start_ids = [x for x in range(len(network.nodes)) if x not in keep_ids]
        if wanted_dict:
            # Keep all the start_ids that is in a wanted datatype.
            x = [
                x for x in start_ids
                if network.nodes[x].datatype.name in wanted_dict]
        elif unwanted_dict:
            # Keep all the start_ids that is not in the unwanted
            # datatype.
            x = [
                x for x in start_ids
                if network.nodes[x].datatype.name not in unwanted_dict]
        x = {}.fromkeys(x)
        keep_ids.update(x)

        # Remove every node not in keep_ids.
        remove_ids = [
            x for x in range(len(network.nodes)) if x not in keep_ids]
        if not remove_ids:
            break
        network = network.delete_nodes(remove_ids)
        network = bie3._OptimizeNoDanglingNodes().optimize(
            network, user_attributes)
        assert network.nodes, "Empty network"
    network = bie3.optimize_network(network, user_attributes)
    return network

def _remove_unwanted_input_nodes(network, unwanted_datatypes, user_attributes):
    # Remove all the input nodes that aren't one of these wanted
    # datatypes.
    # unwanted_datatypes is a list of names.
    return _delete_input_nodes(
        network, None, unwanted_datatypes, user_attributes)


def _keep_wanted_input_nodes(network, wanted_datatypes, user_attributes):
    # Keep only the input nodes that represents one of these wanted
    # datatypes.
    # wanted_datatypes is a list of names.
    return _delete_input_nodes(
        network, wanted_datatypes, None, user_attributes)


def _list_differences_in_nodelists(nodes1, nodes2):
    # nodes1 and nodes2 are lists of nodes.  Should have same
    # datatypes.  Return a list of tuples (node index, attribute name)
    # that are different.
    assert len(nodes1) == len(nodes2)
    diffs = []
    for i, (n1, n2) in enumerate(zip(nodes1, nodes2)):
        assert n1.datatype.name == n2.datatype.name
        for name in n1.attributes:
            v1, v2 = n1.attributes[name], n2.attributes[name]
            x = i, name
            if type(v1) is type("") and type(v2) is type(""):
                if v1 != v2:
                    diffs.append(x)
            elif type(v1) is type("") and type(v2) is not type(""):
                if v1 not in v2:
                    diffs.append(x)
            elif type(v1) is not type("") and type(v2) is type(""):
                if v2 not in v1:
                    diffs.append(x)
            elif type(v1) is not type("") and type(v2) is not type(""):
                if v1 != v2:
                    diffs.append(x)
            else:
                raise AssertionError
    return diffs


def _merge_nodelists(nodes1, nodes2):
    from Betsy import bie3

    assert len(nodes1) == len(nodes2)
    merged = []
    for (n1, n2) in zip(nodes1, nodes2):
        assert n1.datatype.name == n2.datatype.name
        attr = {}
        for name in n1.attributes:
            v1, v2 = n1.attributes[name], n2.attributes[name]
            if type(v1) is type("") and type(v2) is type(""):
                value = v1
                if v1 != v2:
                    value = [v1, v2]
            elif type(v1) is type("") and type(v2) is not type(""):
                value = v1
                if v1 not in v2:
                    value = v2 + [v1]
            elif type(v1) is not type("") and type(v2) is type(""):
                value = v1
                if v2 not in v1:
                    value = v1 + [v2]
            elif type(v1) is not type("") and type(v2) is not type(""):
                value = v1
                if v1 != v2:
                    value = sorted({}.fromkeys(v1 + v2))
            attr[name] = value
        n = bie3.DataNode(n1.datatype, **attr)
        merged.append(n)
    return merged


def _parse_args(args):
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
            # Make sure key not in these parameters already.
            x = [x[1] for x in in_parameters[index] if x[0] == key]
            x = {}.fromkeys(x)  # sorted list of unique values.
            assert not x, "Duplicate --dattr: %s" % key
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
            assert input_or_output == "--input", \
                   "--input_file must be after --input and before --output"
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


def main():
    import os
    import sys
    import argparse
    #import shutil
    import getpass
    import copy

    from genomicode import parselib
    from Betsy import config
    from Betsy import module_utils
    from Betsy import rule_engine_bie3
    from Betsy import userfile
    from Betsy import reportlib
    from Betsy import rulebase
    from Betsy import bie3

    WORKFLOW = (
        "1.  Look through the rulebase for the DataType of interest.\n"
        "2.  Enter the output DataType.\n"
        "FLAG: --output\n"
        "3.  Browse through the list of possible combinations of input\n"
        "    DataTypes.\n"
        "4.  Select one or more input DataTypes.\n"
        "    Now shows just the combinations that includes the input\n"
        "    DataTypes requested.\n"
        "FLAG: --inputs_complete\n"
        "5.  Work out the attributes of the inputs.\n"
        "    Shows the detailed information about each node.\n"
        "FLAG: --attrs_complete\n"
        "6.  System makes sure each node can be found in the network.\n"
        "    If not, try to diagnose differences.\n"
        "7.  System makes sure this is a complete set of input nodes.\n"
        "8.  System makes sure all required module attributes are given.\n"
        "FLAG:--run\n"
        "9.  System makes sure all input files are provided.\n"
        "10.  Actually run the analysis.\n"
        )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Hi!  I'm BETSY, and I'm here to do bioinformatics.\n\n"
        "Workflow:\n%s\n\n" % WORKFLOW)
    parser.add_argument(
        "--inputs_complete",  action="store_const", const=True,
        help="Have finished specifying --input DataTypes.")
    parser.add_argument(
        "--attrs_complete",  action="store_const", const=True,
        help="Have finished configuring --dattr attributes.")
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
    DEFAULT_MAX_INPUTS = 5
    parser.add_argument(
        '--max_inputs', default=DEFAULT_MAX_INPUTS, type=int,
        help="Maximum number of inputs to be shown (default %d).  "
        "(For optimization)" % DEFAULT_MAX_INPUTS)
    parser.add_argument(
        "--save_failed_data", action="store_true",
        help="If a module failed, do not clean up its working files.")

    group = parser.add_argument_group(title="Input/Output Nodes")
    group.add_argument(
        '--input', action='append', help='DataType of the input')
    group.add_argument(
        '--input_file', action='append',
        help='File corresponding to the previous --input.')
    group.add_argument(
        '--mattr', default=[], action='append',
        help='Set the option for a module.  Format should be: '
        '<key>=<value>.')
    group.add_argument(
        "--exclude_input", action="append", help="Experimental.")
    #group.add_argument(
    #    "--prune_other_inputs", action="store_true",
    #    help="There should be no other input data types in this network.  "
    #    "(For optimization)")

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
    #group.add_argument(
    #    '--unpruned_network_png', help='generate the output network png file')
    group.add_argument(
        '--sparse_network_png', action="store_true",
        help="Leave out details in network plot.")
    #group.add_argument(
    #    '--unpruned_network_text', help='')
    #group.add_argument(
    #    '--network_text', help='generate the output network text file')
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
    input_list, x = _parse_args(sys.argv)
    outtype, out_identifier, out_attributes = x

    args.clobber = True
    assert args.num_cores > 0,'num_cores should be greater than 0'
    assert args.max_inputs > 0 and args.max_inputs < 20
    if args.inputs_complete or args.attrs_complete or args.run:
        if args.run:
            x = "--run"
        elif args.attrs_complete:
            x = "--attrs_complete"
        else:
            x = "--inputs_complete"
        assert input_list, "%s given, but no --input." % x
    # TODO: Make sure args.exclude_input is valid.
    # args.exclude_input


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
            # Can be empty if no file given by user.
            assert os.path.exists(identifier), \
                   "File not found: %s" % identifier
    if args.output_file:
        args.output_file = os.path.realpath(args.output_file)
        if os.path.exists(args.output_file):
            assert args.clobber, "Output already exists: %s" % args.output_file

    # Making the IdentifiedDataNode objects.
    if input_list:
        print "Parsing input options."
    # List of identified data nodes.
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

    # Cache the files in the user's directory.  Don't depend on the
    # user keeping the file in the same place.  Needed for the first
    # module.
    print "Making a local copy of the input files."
    sys.stdout.flush()
    for i_data_node in in_data_nodes:
        filename = i_data_node.identifier
        if not filename or not os.path.exists(filename):
            continue
        x = userfile.set(getpass.getuser(), filename)
        i_data_node.identifier = x

    # Parse out the module attributes (or user_options).
    user_options = {}
    if args.mattr:
        print "Parsing module attributes."
        all_mattrs = get_all_option_names()
    for x in args.mattr:
        assert '=' in x, "--mattr should be in format: <option>=<value>"
        key, value = x.split('=', 1)
        assert key in all_mattrs, "Unknown module attribute: %s" % key
        user_options[key] = value

    # Possible cases.
    # 1.  No inputs and no outputs.
    #     Show the whole rulebase.
    # 2.  No inputs and has output.  (requires network)
    #     Show the input datatypes.
    # 3.  Input and no outputs.
    #     ERROR
    # 4.  Input and output.  (requires network)
    #     Show the input nodes.
    #     Run the analysis.

    # Case 1.  No inputs and no outputs.
    if not args.output and not args.input:
        # If there's no inputs and no outputs, just print the whole
        # rulebase.
        print
        print "Selecting --output DataType."
        print_rulebase(rulebase)
        return
    # Case 3.  Input and no outputs.
    if args.input and not args.output:
        raise AssertionError, "Missing --output"

    print "Generating a network that produces a %s..." % outtype
    sys.stdout.flush()
    # There's a bug in here somehow where impossible networks can be
    # created.  e.g. FastqFolder:orientation="unknown" -> merge_reads
    # -> FastqFolder:orientation="single", even though constraint
    # forces them to be the same.  Happens during complete_network or
    # optimize_network step.  Don't see anymore, because got rid of
    # orientation attribute in FastqFolder.
    network = bie3.backchain(
        rulebase.all_modules, out_datatype, user_attributes)
    network = bie3.complete_network(network, user_attributes)
    network = bie3.optimize_network(network, user_attributes)
    assert network, "Unexpected error: No network generated."
    print "Made a network with %d nodes." % len(network.nodes)

    #plot_network(
    #    args.network_png, network, user_options=user_options,
    #    verbose=(not args.sparse_network_png))
    #return

    network_orig = copy.deepcopy(network)

    # If the user specified some data types that should be excluded,
    # remove them from the network.
    if args.exclude_input:
        if len(args.exclude_input) <= 3:
            x = ", ".join(args.exclude_input)
            print "Remove as inputs: %s." % x
        else:
            print "Remove %d data types as inputs." % len(args.exclude_input)
        network = _remove_unwanted_input_nodes(
            network, args.exclude_input, user_attributes)


    ## if args.unpruned_network_png:
    ##     print "Plotting unpruned network: %s." % args.unpruned_network_png
    ##     bie3.plot_network_gv(
    ##         args.unpruned_network_png, network, options=user_options,
    ##         verbose=(not args.sparse_network_png))
    ## if args.unpruned_network_text:
    ##     print "Writing detailed network: %s." % args.unpruned_network_text
    ##     bie3.print_network(network, outhandle=args.unpruned_network_text)


    # Case 2.  No input and has output.
    if (not args.inputs_complete and not args.attrs_complete and not args.run)\
           or not in_data_nodes:
        # If there's an output and no inputs, then show the possible
        # data types.  The network is probably too big to show
        # everything.
        print
        print "Selecting --input DataTypes.  Add --inputs_complete when done."
        # TODO: Should not show the same datatype as the output.
        # e.g. if output is a FastQCFolder, don't show an input that
        # just consists of a FastQCFolder.  That's just lame.
        x = [x.data.datatype.name for x in in_data_nodes]
        print_input_datatypes(
            network, x, args.exclude_input, user_attributes,
            max_inputs=args.max_inputs)
        plot_network(
            args.network_png, network, user_options=user_options,
            verbose=(not args.sparse_network_png))
        return
    # Case 4.  Has input and output.
    assert args.output and in_data_nodes

    # If the user has already finalized the input data types, keep
    # only these as input nodes in the network.
    x = [x.data.datatype.name for x in in_data_nodes]
    x = {}.fromkeys(x).keys()
    assert x

    # This can lead to an empty network.  This can happen if user
    # specifies --inputs_complete, but inputs were not complete.
    network = _keep_wanted_input_nodes(network, x, user_attributes)

    # Configure the attributes if the user is not ready to run yet.
    if not args.attrs_complete and not args.run:
        print
        print "Configuring input data.  Add --attrs_complete when done."
        x = [x.data.datatype.name for x in in_data_nodes]
        print_input_nodes(
            network, x, args.exclude_input, user_attributes,
            max_inputs=args.max_inputs)
        plot_network(
            args.network_png, network, user_options=user_options,
            verbose=(not args.sparse_network_png))
        return

    # See if all the inputs can be found in the network.
    print "Searching for each input in the network..."
    sys.stdout.flush()
    start_node_ids = []
    missing = []  # list of indexes into in_data_nodes that are not found.
    for i, node in enumerate(in_data_nodes):
        x = bie3._find_start_nodes(network, node.data)
        start_node_ids.extend(x)
        if not x:
            missing.append(i)

    # If any nodes can't be found in the network, try to diagnose it.
    if missing:
        print
        print "Diagnosing problem with input attributes."
        x = [x.data for x in in_data_nodes]
        print_diagnose(network, x, missing)
        plot_network(
            args.network_png, network, user_options=user_options,
            highlight_node_ids=start_node_ids,
            verbose=(not args.sparse_network_png))
        return
    else:
        print "all input nodes found."
    sys.stdout.flush()

    ## if 0:  # This takes way too long.
    ##     # If this is taking too long, can speed it up by simplifying
    ##     # the network with more constraints.
    ##     if not is_complete_input(network, user_attributes, start_node_ids):
    ##         print_input_nodes(network, user_attributes)
    ##         return

    ## # Actually, don't prune the network.  It can be helpful for
    ## # diagnosing problems.
    ## # Prune network based on these inputs.
    ## x = [x.data for x in in_data_nodes]
    ## x = [x.datatype.name for x in x]
    ## x = sorted({}.fromkeys(x))
    ## x_str = ", ".join(x)
    ## input_or_inputs = "input"
    ## if len(x) > 1:
    ##     input_or_inputs = "inputs"
    ## x = "Pruning the network for paths that start from %s: %s." % (
    ##     input_or_inputs, x_str)
    ## parselib.print_split(x, prefixn=2)
    ## sys.stdout.flush()
    ## x = [x.data for x in in_data_nodes]
    ## p_network = bie3.select_start_node(network, x, user_attributes)
    ## if not p_network.nodes:
    ##     print "No network.  Problem with --inputs somehow."
    ##     plot_network(
    ##         args.network_png, network, user_options=user_options,
    ##         highlight_node_ids=start_node_ids,
    ##         verbose=(not args.sparse_network_png))
    ##     return
    ## network = p_network

    ## # Network has changed.  Need to find the new start_node_ids.
    ## start_node_ids = []
    ## for node in in_data_nodes:
    ##     x = bie3._find_start_nodes(network, node.data)
    ##     print "HERE 7", x
    ##     assert x
    ##     start_node_ids.extend(x)

    print "Making sure all required --mattr are provided..."
    # Make sure all required mattr are provided.  This can only be run
    # after the final network is generated.
    modules = [
        x for x in network.nodes if isinstance(x, bie3.ModuleNode)]
    opt2mods = get_required_option_names(modules)
    x = opt2mods.keys()
    x = [x for x in x if x not in user_options]
    for on in x:
        for mn in opt2mods[on]:
            print 'Missing --mattr: %s requires attribute "%s".' % (mn, on)
    if x:
        plot_network(
            args.network_png, network, user_options=user_options,
            highlight_node_ids=start_node_ids,
            verbose=(not args.sparse_network_png))
        return
    print "All --mattr are specified or have defaults."

    print "Success!  I have a network."
    print "The final network has %d nodes." % len(network.nodes)

    #if args.network_png:
    #    print "Plotting network: %s." % args.network_png
    #    bie3.plot_network_gv(
    #        args.network_png, network, options=user_options,
    #        verbose=(not args.sparse_network_png))
    #if args.network_text:
    #    print "Writing detailed network: %s." % args.network_text
    #    bie3.print_network(network, outhandle=args.network_text)
    #if args.network_json:
    #    print "Writing network in json format: %s." % args.network_json
    #    bie3.write_network(args.network_json, network)

    print "Looking for input files...."
    missing = False
    for x in in_data_nodes:
        if not x.identifier:
            print "no file given for: %s." % (x.data.datatype.name)
            missing = True
        elif not os.path.exists(x.identifier):
            print "File not found: %s." % x.identifier
            missing = True
    if missing:
        plot_network(
            args.network_png, network, user_options=user_options,
            highlight_node_ids=start_node_ids,
            verbose=(not args.sparse_network_png))
        return
    print "all input files found."

    if not args.run:
        print ("Please review the network to make sure the you agree with "
               "the analysis.")
        print "Add --run when ready to run network."
        plot_network(
            args.network_png, network, user_options=user_options,
            highlight_node_ids=start_node_ids,
            verbose=(not args.sparse_network_png))
        return

    print "Running the analysis."
    clean_up = not args.save_failed_data
    node_dict = output_file = None
    try:
        x = rule_engine_bie3.run_pipeline(
            network, in_data_nodes, user_attributes, user_options,
            user=args.user, job_name=args.job_name, clean_up=clean_up,
            num_cores=args.num_cores)
        if x:
            node_dict, output_file = x
    except AssertionError, x:
        if str(x).startswith("Inference error"):
            node_ids = rule_engine_bie3.DEBUG_POOL.keys()
            plot_network(
                args.network_png, network, user_options=user_options,
                highlight_node_ids=node_ids,
                verbose=(not args.sparse_network_png))
        raise 

    # Print the original network, showing the path that was taken.
    # Find the node_ids in the original network.
    if not node_dict:
        # Pipeline was not successful.  Print the filtered network to
        # facilitate diagnosing.
        plot_network(
            args.network_png, network, user_options=user_options,
            highlight_node_ids=start_node_ids,
            verbose=(not args.sparse_network_png))
    else:
        # Successful.  Show the full network.
        node_ids = []
        for inode in node_dict.itervalues():
            x = bie3.find_node(network_orig, inode.data)
            assert x is not None
            node_ids.append(x)
        plot_network(
            args.network_png, network_orig, user_options=user_options,
            highlight_node_ids=node_ids,
            verbose=(not args.sparse_network_png))

    if not output_file or not args.output_file:
        return


    print "Saving results at %s." % args.output_file
    reportlib.copy_file_or_path(output_file, args.output_file)
    ## # TODO: Use reportlib.copy_file_or_path.
    ## # Remove the previous files if they exist.
    ## if args.clobber and os.path.exists(args.output_file):
    ##     if os.path.isdir(args.output_file):
    ##         shutil.rmtree(args.output_file)
    ##     else:
    ##         os.unlink(args.output_file)
    ## # Copy the file or directory.
    ## if os.path.isdir(output_file):
    ##     shutil.copytree(output_file, args.output_file)
    ## else:
    ##     shutil.copy(output_file, args.output_file)


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
