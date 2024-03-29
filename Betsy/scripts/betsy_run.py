#! /usr/bin/env python

# Functions:
# check_output_provided      Step 1: Make sure there's an output.
#   _print_rulebase          To help select output DataType.
#   _pretty_print_datatype
#   _pretty_print_module
# generate_network           Step 2.
# check_more_than_one_node_network    Step 2.5.
# check_inputs_provided      Step 3: Make sure any --inputs are provided.
#   _print_input_datatypes   To help select input DataType(s).
# check_inputs_in_network    Step 4.
#   _print_node_score_table
# build_pipelines            Step 5.
# check_attributes_complete  Step 6.
# prune_pipelines            Step 7.
#   _print_attributes_pruned
# check_input_files          Step 8.
# check_output_file          Step 9.
# manually_verify_network    Step 10.
# 
# plot_network
# plot_network_show_pipelines
# plot_pipelines
# write_receipt
# 
# get_all_option_names       Get the names of all OptionDef in the rulebase.
# get_module_options         Get the options used in a list of modules.
#
# _list_differences_in_nodelists
# _merge_nodelists
# _find_lowest_datatype
# _parse_dattr
# _parse_args                 Parse the arguments from the user.
# _make_custom_attr


def check_output_provided(rulebase, output_file):
    if output_file:
        return True
    # If there is no output, just print the whole rulebase.
    #print "Selecting --output DataType."
    _print_rulebase(rulebase)
    return False


def _print_rulebase(rulebase):
    from Betsy import bie3

    # Make a list of the DataType objects.
    x = [getattr(rulebase, x) for x in dir(rulebase)]
    x = [x for x in x if isinstance(x, bie3.DataType)]
    datatypes = x

    # Count the attributes.
    num_attributes = 0
    for dt in datatypes:
        num_attributes += len(dt.attribute_defs)

    # Make a list of the modules
    modules = rulebase.all_modules

    print "Rulebase contains %d data types, %d attributes, and %d modules." % (
        len(datatypes), num_attributes, len(modules))

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


def generate_network(rulebase, outtype,
                     custom_attributes, out_custom_attribute):
    import sys
    from Betsy import bie3

    assert outtype

    # Get the out_datatype.
    # BUG: Should separate the custom attributes for the out_datatype.
    assert hasattr(rulebase, outtype), "Unknown datatype: %s" % outtype
    attrs = {}
    if out_custom_attribute:
        assert out_custom_attribute.datatype.name == outtype
        for x in out_custom_attribute.attributes:
            assert x.name not in attrs
            attrs[x.name] = x.value
    out_datatype = getattr(rulebase, outtype)
    out_data = out_datatype.output(**attrs)
    #for cattr in custom_attributes:
    #    if cattr.datatype.name != out_datatype.name:
    #        continue
    #    for x in cattr.attributes:
    #        assert x.name not in attrs
    #        attrs[x.name] = x.value

    a_or_an = "a"
    if outtype.lower()[0] in "aeiou":
        a_or_an = "an"
    print "Generating a network that produces %s:\n%s" % (
        a_or_an, out_data.datatype.name)
    for name in sorted(out_data.attributes):
        value = out_data.attributes[name]
        x = "  %s=%s" % (name, value)
        print x
    
    sys.stdout.flush()

    # There may be a bug in here somehow where impossible networks can
    # be created.  e.g. FastqFolder:orientation="unknown" ->
    # merge_reads -> FastqFolder:orientation="single", even though
    # constraint forces them to be the same.  Happens during
    # complete_network or optimize_network step.  Don't see anymore,
    # because got rid of orientation attribute in FastqFolder.
    network = bie3.make_network(
        rulebase.all_modules, out_data, custom_attributes)
    assert network, "Unexpected error: No network generated."

    ## # If the user specified some data types that should be excluded,
    ## # remove them from the network.
    ## if exclude_input:
    ##     if len(exclude_input) <= 3:
    ##         x = ", ".join(exclude_input)
    ##         print "Remove as inputs: %s." % x
    ##     else:
    ##         print "Remove %d data types as inputs." % len(exclude_input)
    ##     network = _remove_unwanted_input_nodes(
    ##         network, exclude_input, user_attributes)
    x = "nodes"
    if len(network.nodes) == 1:
        x = "node"
    print "Made a network with %d %s." % (len(network.nodes), x)
    return network


def check_more_than_one_node_network(network, rulebase):
    #from genomicode import parselib
    from Betsy import bie3
    
    if len(network.nodes) > 1:
        return True
    
    # Print out the node.
    node = network.nodes[0]
    print "I cannot generate this network."
    #print node.datatype.name
    #for x in node.attributes.iteritems():
    #    name, value = x
    #    x = "  %s=%s" % (name, value)
    #    print x
    print "This can be caused by:"
    print "1.  An incompatibility in the attributes (most likely)."
    print "    Please verify the attributes in the output data object."
    print "2.  The knowledge base is incomplete and missing a rule."

    # Look for rules that can almost generate this node, with at most
    # 1 mismatch.
    attr2values = {}  # attr -> list of values
    for module in rulebase.all_modules:
        assert not bie3._is_valid_output(module, node, save_conflicts=True)
        # Print out conflicts
        conflicts = bie3.DEBUG_IS_VALID_OUTPUT_CONFLICTS
        if not conflicts:  # Can happen if data type doesn't match.
            continue
        if len(conflicts) > 1:
            continue
        attr_name, desired_value, data_value = conflicts[0]
        #print module.name, attr_name, desired_value, data_value
        v1 = attr2values.get(attr_name, [])
        v2 = desired_value
        if type(v2) is type(""):
            v2 = [v2]
        attr2values[attr_name] = v1+v2
    if not attr2values:
        return False

    print "Is it possible that you meant:"
    i = 1
    for attr_name in sorted(attr2values):
        x = attr2values[attr_name]
        x = {}.fromkeys(x)
        values = sorted(x)
        #values_str = parselib.pretty_list(values, conjunction="or")
        dtype = node.datatype.name
        for value in values:
            print "%d.  %s.%s=%s" % (i, dtype, attr_name, value)
            i += 1
    return False


def check_inputs_provided(
    network, in_data_nodes, custom_attributes, user_options,
    max_inputs, network_png, verbose):
    if in_data_nodes:
        return True
    # If there's an output and no inputs, then show the possible
    # data types.  The network is probably too big to show
    # everything.
    x = [x.data.datatype.name for x in in_data_nodes]
    _print_input_datatypes(
        network, x, custom_attributes, max_inputs=max_inputs)
    plot_network(
        network_png, network, user_options=user_options, verbose=verbose)
    return False


def _print_input_datatypes(
    network, required_datatypes, custom_attributes,
    max_inputs=None, outhandle=None):
    # required_datatypes is a list of the names of datatypes that must
    # be in this combination.
    import sys
    from genomicode import parselib
    from Betsy import bie3

    outhandle = outhandle or sys.stdout

    ps = parselib.print_split

    print >>outhandle, \
          "Looking for --input DataTypes that can generate this output."
    outhandle.flush()

    excluded_datatypes = None
    datatype_combos = bie3.get_input_datatypes(
        network, custom_attributes, skip_datatypes=excluded_datatypes,
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

    # Don't do this, because sometimes you want to convert between the
    # same data type (e.g. to convert counts to cpm in an
    # UnprocessedSignalFile).
    ## If there are other data types, then exclude the same datatype as
    ## the output.  e.g. if output is a FastQCFolder, don't show an
    ## input that just consists of a FastQCFolder.  That's just lame.
    #out_dname = network.nodes[0].datatype.name
    #i = 0
    #while i < len(datatype_combos):
    #    combo = datatype_combos[i]
    #    if len(combo) == 1 and combo[0].name == out_dname:
    #        del datatype_combos[i]
    #    else:
    #        i += 1

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

    x1 = "These are"
    x2 = ""
    x3 = "combinations"
    if len(datatype_combos) == 1:
        x1 = "This is"
        x2 = "only "
        x3 = "combination"
    x4 = ""
    if max_inputs is not None:
        x4 += " (up to %d)" % max_inputs
    x5 = "%s the %sacceptable %s of --input DataTypes%s.  " % (x1, x2, x3, x4)
    x6 = "The exact number of each DataType needed is not shown."
    x = x5 + x6
    ps(x, prefixn=0, outhandle=outhandle)
    for i, combo in enumerate(datatype_combos):
        names = [x.name for x in combo]
        for j in range(len(names)):
            if names[j] not in required_datatypes:
                names[j] = "%s*" % names[j]
        ps("%3d.  %s" % (
            i+1, ", ".join(names)), prefixn=6, outhandle=outhandle)


def check_inputs_in_network(
    network, user_options, in_data_nodes, network_png, verbose):
    # Return a tuple of:
    # all_inputs_found (boolean)
    # list of (list of node IDs that match in_data_nodes).
    import sys
    import itertools
    from Betsy import bie3
    from genomicode import jmath

    print "Assigning --input's to nodes in the network."
    sys.stdout.flush()

    data_nodes = [x.data for x in in_data_nodes]
    
    # index of data_nodes -> list of compatible node_ids
    node2ids = [bie3._find_compat_data(network, x) for x in data_nodes]

    complete = []    # list of indexes into in_data_nodes
    incomplete = []
    for i, node_ids in enumerate(node2ids):
        if node_ids:
            complete.append(i)
        else:
            incomplete.append(i)

    # Every node can be assigned to the network.
    if len(complete) == len(in_data_nodes):
        # Make sure there is a unique assignment of in_data_nodes to
        # nodes in the network, e.g. in_data_nodes are not assigned to
        # the same node.  Test every possible combination of
        # assignment.
        x = [len(x) for x in node2ids]
        num_combos = jmath.prod(x)
        assert num_combos < 100000, "Too many possibilities"
        found = False
        for x in itertools.product(*node2ids):
            num_nodes = len({}.fromkeys(x))
            if num_nodes == len(in_data_nodes):
                found = True
                break
        if not found:
            print "Ambiguous or duplicated --inputs."
            return False, node2ids
        # Everything looks OK.
        print "All --input data found in network."
        return True, node2ids

    # For each in_data_node that cannot be found, see if we can figure
    # out the closest match.
    scores = []
    for data_index in incomplete:
        x = bie3._score_compat_data(network, data_nodes[data_index])
        # If x is empty, then there are no matching data types in the
        # network.
        print "--input %s is not in the network." % \
              data_nodes[data_index].datatype.name
        scores.extend(x)
    if scores:
        _print_node_score_table(network, scores)

    # Make a list of all possible start_ids.
    start_ids = bie3._uniq(bie3._flatten(node2ids))
    #for i, node_ids in enumerate(node2ids):
    #    start_ids.extend(node_ids)
    #start_ids = _uniq(start_ids)
    plot_network(
        network_png, network, user_options=user_options,
        highlight_green=start_ids, verbose=verbose)
    return False, node2ids


def _print_node_score_table(network, scores):
    # scores from bie3._score_compat_data

    # Figure out the score cutoff for each data type.
    dt2scores = {}
    for x in scores:
        score, node_id, user_data, attr_values = x
        dtname = network.nodes[node_id].datatype.name
        if dtname not in dt2scores:
            dt2scores[dtname] = []
        dt2scores[dtname].append(score)
    # Score cutoff is minimum score + 1
    score_cutoffs = {}  # dtname -> max score to print
    for dtname, score in dt2scores.iteritems():
        score_cutoffs[dtname] = min(score) + 1

    
    # Make an output table.
    table = []
    header = ["Node", "D", "Datatype", "Attribute", "Yours", "Network"]
    table.append(header)

    for x in scores:
        score, node_id, user_data, attr_values = x
        dt_name = network.nodes[node_id].datatype.name
        if not attr_values:
            x = node_id, score, dt_name, "", "", ""
            assert len(x) == len(header)
            table.append(x)
        for name, netw_value, user_value in attr_values:
            x = node_id, score, dt_name, name, user_value, netw_value
            assert len(x) == len(header)
            if score > score_cutoffs[dt_name]:
                continue
            table.append(x)

    # Figure out the maximum lengths of each column.
    num_cols = len(header)
    col_lengths = []
    for i in range(num_cols):
        x = [x[i] for x in table]     # Get values in column.
        x = [len(str(x)) for x in x]  # calculate lengths.
        x = max(x)
        col_lengths.append(x)

    # Set a maximum limit for Datatype and Attribute columns.
    # Should be max 79 columns long, including 5 for spaces.
    # Column      MIN  MAX  Notes
    # Node         4     4  As short as possible.
    # Delta        1     5  As short as possible.
    # Datatype     8    18  Datatype name.
    # Attribute    9    19  Attribute name.
    # User         4    10  Can be long.  But usually short.
    # Network      7    22  Might be long.  Lots of values.
    max_lengths = [4, 1, 18, 19, 10, 22]
    assert len(col_lengths) == len(max_lengths)
    # Just use the maximum lengths.
    col_lengths = max_lengths
    #for i in range(len(col_lengths)):
    #    col_lengths[i] = min(col_lengths[i], max_lengths[i])
    # Make sure the values aren't too long.
    for i in range(len(table)):
        row = list(table[i])
        for j in range(len(row)):
            x = row[j]
            x = str(x)
            x = x.rstrip()
            if len(x) > col_lengths[j]:
                x = x[:col_lengths[j]-3] + "..."
            row[j] = x
        table[i] = row

    fmt = "{!s:^%ds} {!s:^%ds} {:<%ds} {:<%ds} {:<%ds} {:<%ds}" % \
          tuple(col_lengths)
    for x in table:
        print fmt.format(*x)


def build_pipelines(
    network, user_options, in_data_nodes, data_node_ids, custom_attributes,
    max_inputs, network_png, verbose):
    # Return list of (path, start_ids).  start_ids is parallel to
    # in_data_nodes.  If no paths found, will return an empty list.
    import sys
    from genomicode import parselib
    from Betsy import bie3

    print "Constructing pipelines that use --input data types."
    sys.stdout.flush()

    # data_node_ids is parallel to in_data_nodes.  Each element is a
    # list of the node_ids that that data node can map onto.
    try:
        paths = bie3.find_paths_by_start_ids(
            network, custom_attributes, data_node_ids)
    except AssertionError, x:
        if str(x).startswith("Too many paths"):
            # Helpful for debugging.
            print "ERROR: Too many pipelines.  Could not generate network."
            plot_network(
                network_png, network, user_options=user_options,
                verbose=verbose)
        raise

    # Make sure all the --inputs are needed.  Any unnecessary ones may
    # indicate a problem.
    inputs_used = {}  # list of indexes of --inputs that are used
    for p in paths:
        used = [i for (i, x) in enumerate(p.start_ids) if x is not None]
        inputs_used.update({}.fromkeys(used))
    has_unused_inputs = len(inputs_used) != len(in_data_nodes)

    good_paths = [x for x in paths if not x.missing_ids]
    if good_paths and not has_unused_inputs:
        x = "pipelines"
        if len(good_paths) == 1:
            x = "pipeline"
        print "Found %d possible %s." % (len(good_paths), x)
        return good_paths

    # Print out --inputs that aren't used.
    if has_unused_inputs:
        for i in range(len(in_data_nodes)):
            if i in inputs_used:
                continue
            name = in_data_nodes[i].data.datatype.name
            x = (
                "%s is not used in any pipelines.  "
                "Please make sure the proposed pipelines are acceptable and "
                "remove this input." % name)
            parselib.print_split(x, prefixn=2)

    if not good_paths:
        print "No pipelines found.  Examine network to diagnose."
        print "Make sure that no --input is missing."
        x = [x.data.datatype.name for x in in_data_nodes]
        _print_input_datatypes(
            network, x, custom_attributes, max_inputs=max_inputs)
        print

        # DEBUG
        #for i, path in enumerate(paths):
        #    print "%2d.  %s %s" % (i, path.start_ids, path.missing_ids)

        # For each in_data_node, see if it might be a better match to
        # another node.
        scores = []
        for i in range(len(in_data_nodes)):
            x = bie3._score_compat_data(network, in_data_nodes[i].data)
            scores.extend(x)
        if scores:
            print ("Make sure that the attributes for the --input's are "
                   "correct.")
            _print_node_score_table(network, scores)
            print

    # Plot out the network.
    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        verbose=verbose)
    return []


def check_attributes_complete(
    network, user_options, paths, network_png, prune_network, verbose):
    # user_options is dict of name to value
    import sys
    from genomicode import parselib
    from Betsy import bie3

    # Make sure all required mattr are provided.  This can only be run
    # after the final network is generated.
    print "Making sure all required attributes (--mattr) are provided."
    sys.stdout.flush()

    assert paths

    all_missing = {}
    all_extra = []
    good_paths = []
    for i, p in enumerate(paths):
        module_ids = [
            x for x in p.node_ids if
            isinstance(network.nodes[x], bie3.ModuleNode)]
        modules = [network.nodes[x] for x in module_ids]

        missing = {}
        all_opt2mods = get_module_options(modules)
        required_opt2mods = get_module_options(modules, required_only=True)
        x = [x for x in required_opt2mods if x not in user_options]
        for on in x:
            for mn in all_opt2mods[on]:
                missing[(mn, on)] = 1
        extra = [x for x in user_options if x not in all_opt2mods]
        if not missing:
            good_paths.append(p)
        all_missing.update(missing)
        all_extra.extend(extra)
    all_extra = {}.fromkeys(all_extra)

    if all_extra:
        names = sorted(all_extra)
        x = parselib.pretty_list(names)
        x = ("The following --mattr options were provided, but may not "
             "be needed: %s" % x)
        parselib.print_split(x, prefixn=2)

    if good_paths:
        if len(good_paths) < len(paths):
            num_removed = len(paths)-len(good_paths)
            assert all_missing
            if len(all_missing) == 1:
                mn, on = all_missing.keys()[0]
                x = ("Removed %d pipelines because --mattr option %r was "
                     "not provided." % (num_removed, on))
            else:
                names = sorted([x[1] for x in all_missing.keys()])
                x = parselib.pretty_list(names)
                x = ("Removed %d pipelines because the following --mattr "
                     "options were not provided: %s" % (num_removed, x))
            parselib.print_split(x, prefixn=2)
        return good_paths
    for (mn, on) in all_missing:
        print 'Missing --mattr: %s requires attribute "%s".' % (mn, on)

    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        prune=prune_network, verbose=verbose)
    return []


def prune_pipelines(
    network, user_options, custom_attributes, paths, 
    network_png, prune_network, verbose):
    # Any any pipelines look weird, then remove them.
    import sys
    from Betsy import bie3

    print "Pruning redundant pipelines."
    sys.stdout.flush()
    paths_orig = paths[:]
    
    ## Optimizations:
    ## 1.  Convert node_ids to dictionaries.
    ## 2.  Convert the transitions to sets for easier comparison.
    #for i, p in enumerate(paths):
    #    # Not much more efficient.
    #    #p.node_ids_set = set(p.node_ids)
    #    p.node_ids = {}.fromkeys(p.node_ids)
    #    for key, value in p.transitions.iteritems():
    #        p.transitions[key] = set(value)

    paths = bie3.prune_paths(paths, network, custom_attributes)
    
    ## Convert node_ids back to lists.
    #for i, p in enumerate(paths):
    #    p.node_ids = p.node_ids.keys()
    #    for key, value in p.transitions.iteritems():
    #        p.transitions[key] = list(value)
            
    if not paths:
        print "All pipelines pruned.  This can happen if:"
        print "  1.  There is a conflicting data attribute in the pipeline."
        print "      Maybe an attribute is set for the wrong DataType?"
        print "  2.  A --mattr option is missing."
        print "  3.  There is a bug in the network generation."
        print "  4.  There is a bug in the pipeline pruning."
        print "Please review the network, attributes, and --mattr options."
        plot_network_show_pipelines(
            network_png, network, paths_orig, user_options=user_options,
            prune=prune_network, verbose=verbose)
        return paths
    
    num_pruned = len(paths_orig) - len(paths)
    if not num_pruned:
        x = "s"
        if len(paths) == 1:
            x = ""
        print "No redundant pipelines found.  %d final pipeline%s." % (
            len(paths), x)
    else:
        x1 = "s"
        if num_pruned == 1:
            x1 = ""
        x2 = "s"
        if len(paths) == 1:
            x2 = ""
        print "Pruned %d pipeline%s.  %d final pipeline%s." % (
            num_pruned, x1, len(paths), x2)

    return paths


def _print_attributes_pruned(pruned_attributes):
    # Print out a message about what was deleted.
    from Betsy import bie3
    from genomicode import parselib
    
    attr2values = {}
    attr2deleted = {}
    attr2user = {}  # whether this was deleted based on user request
    for x in pruned_attributes:
        datatype_name, attr_name, kept_value, del_value, user = x
        name = "%s.%s" % (datatype_name, attr_name)
        if name not in attr2values:
            attr2values[name] = []
        if name not in attr2deleted:
            attr2deleted[name] = []
        attr2values[name].append(kept_value)
        attr2values[name].append(del_value)
        attr2deleted[name].append(del_value)

        # Multiple values may be deleted for each attribute
        # (e.g. multiple aligners may be deleted.  If any one is
        # deleted based on a user attribute, then this should be True.
        if user or not name in attr2user:
            attr2user[name] = user
        
    for name in sorted(attr2values):
        all_values = sorted(bie3._uniq(attr2values[name]))
        kept_values = [x for x in all_values if x not in attr2deleted[name]]
        assert len(kept_values) == 1
        kept_value = kept_values[0]
        x1 = "%s can be %s." % (name, repr(all_values))
        x2 = "Arbitrarily using %s." % repr(kept_value)
        if attr2user[name]:
            x2 = "Using %s because it was specified by the user." % \
                 repr(kept_value)
        x = "%s  %s" % (x1, x2)
        parselib.print_split(x, prefixn=2)


def _dict_diff(dict1, dict2, _as_dict=False):
    # Only elements that are in dict1 and not in dict2.
    # Using set comparisons does not make this any faster than dicts.
    if _as_dict:
        dict3 = {}
        for x in dict1:
            if x not in dict2:
                dict3[x] = 1
    else:
        dict3 = [x for x in dict1 if x not in dict2]
    return dict3


def check_input_files(network, in_data_nodes, user_options, paths,
                      network_png, prune_network, verbose):
    import os
    import sys

    print "Making sure all input files provided."
    sys.stdout.flush()
    missing = False
    for x in in_data_nodes:
        if x.data.datatype.no_file:
            # This DataType does not require a file.
            continue
        elif not x.identifier:
            print "No file given: %s." % (x.data.datatype.name)
            missing = True
        elif not os.path.exists(x.identifier):
            print "File not found: %s." % x.identifier
            missing = True

    if not missing:
        #print "All input files found."
        #sys.stdout.flush()
        return True

    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        prune=prune_network, verbose=verbose)
    return []


def check_output_file(filename):
    if not filename:
        print "No --output_file specified.  Will not save results."

    return True


def manually_verify_network(
    network, user_options, paths, run, network_png, prune_network, verbose):
    import sys
    if run:
        return True

    print "Please review the network to make sure the analysis is correct."
    print "Add --run when ready to run the analysis."
    sys.stdout.flush()

    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        prune=prune_network, verbose=verbose)
    return False
    

def plot_network_show_pipelines(filename, network, paths, **keywds):
    # Pass keywds through to plot_network.
    from Betsy import bie3

    x1 = [x.node_ids for x in paths]
    # Hack: For optimization, sometimes node_ids is dict.
    if x1 and type(x1[0]) is type({}):
        x1 = [x.keys() for x in x1]
    # Hack: For optimization, sometimes node_ids is frozenset.
    if x1 and type(x1[0]) is frozenset:
        x1 = [list(x) for x in x1]
    x2 = [x.start_ids for x in paths]
    x1 = bie3._uniq(bie3._flatten(x1))
    x2 = bie3._uniq(bie3._flatten(x2))
    x2 = [x for x in x2 if x is not None]
    all_node_ids = x1
    all_start_ids = x2

    # Pull out the missing IDs from the pathways.
    x = [x.missing_ids for x in paths]
    if x and type(x[0]) is type({}):
        x = [x.keys() for x in x]
    elif x and type(x[0]) is frozenset:
        x = [list(x) for x in x]
    x = bie3._uniq(bie3._flatten(x))
    missing_ids = x

    if not missing_ids:
        # Find nodes with no parents in the network that aren't start_ids.
        nodeid2parents = bie3._make_parents_dict(network)
        no_parents = {}
        for node_id in all_node_ids:
            x = nodeid2parents.get(node_id, [])
            x = [x for x in x if x in all_node_ids]
            if not x:
                no_parents[node_id] = 1
        x = all_node_ids
        x = [x for x in x if x in no_parents]   # no parents
        x = [x for x in x if x not in all_start_ids]  # not start id
        missing_ids = x

    transitions = {}
    for path in paths:
        for node_id, next_ids in path.transitions.iteritems():
            for next_id in next_ids:
                transitions[(node_id, next_id)] = 1
                
    highlight_yellow = None
    if not keywds.get("prune"):
        x = all_node_ids
        x = [x for x in x if x not in all_start_ids]  # not start_ids
        x = [x for x in x if x not in missing_ids]
        highlight_yellow = x

    plot_network(
        filename, network, bold=all_node_ids,
        bold_transitions=transitions, highlight_green=all_start_ids,
        highlight_orange=missing_ids, highlight_yellow=highlight_yellow,
        **keywds)


def plot_network(
    filename, network, user_options=None,
    bold=[], bold_transitions=[], 
    highlight_green=[], highlight_orange=[], highlight_purple=[],
    highlight_yellow=[], prune=False, verbose=False):
    # If show_node_ids is not None, should be a list of dict of the
    # node_ids to include.
    import sys
    from Betsy import bie3

    if filename is None:
        return
    print "Plotting (%d node) network to %s." % (len(network.nodes), filename)
    sys.stdout.flush()

    show_node_ids = None
    if prune:
        x = set()
        if bold:
            x.update(bold)
        if highlight_green:
            x.update(highlight_green)
        if highlight_orange:
            x.update(highlight_orange)
        if highlight_purple:
            x.update(highlight_purple)
        if highlight_yellow:
            x.update(highlight_yellow)
        for (y, z) in bold_transitions:
            x.update([y, z])
        show_node_ids = x
    
    bie3.plot_network_gv(
        filename, network, options=user_options, bold=bold,
        bold_transitions=bold_transitions,
        highlight_green=highlight_green, highlight_orange=highlight_orange,
        highlight_purple=highlight_purple, highlight_yellow=highlight_yellow,
        show_node_ids=show_node_ids, verbose=verbose)


def write_network(filename, network):
    from Betsy import bie3

    if not filename:
        return
    print "Writing network in json format: %s." % filename
    bie3.write_network(filename, network)


def plot_pipelines(filestem, network, paths, user_options, max_pipelines=None,
                   prune=False, verbose=False):
    if max_pipelines is not None:
        paths = paths[:max_pipelines]
    for i, path in enumerate(paths):
        filename = "%s-%02d.png" % (filestem, i)
        plot_network_show_pipelines(
            filename, network, [path], user_options=user_options,
            prune=prune, verbose=verbose)


def write_receipt(outfilename, network, start_ids, transitions, node_dict):
    import os
    import sys
    from genomicode import parselib
    from Betsy import bie3
    from Betsy import rule_engine
    from Betsy import module_utils as mlib

    print "Writing receipt to %s." % outfilename
    
    # Figure out the order to write the nodes.
    node_order = []
    stack = start_ids
    # Do a breadth first search.
    while stack:
        nid = stack.pop(0)
        if nid in node_order:
            continue
        node_order.append(nid)
        for (n1, n2) in transitions:
            if nid == n1:
                stack.append(n2)

    # No module nodes.
    node_order = [
        x for x in node_order
        if not isinstance(network.nodes[x], bie3.ModuleNode)]

    handle = open(outfilename, 'w')

    # Write the command line.
    print >>handle, "COMMAND"
    print >>handle, "-------"
    cmd = " ".join(sys.argv)
    parselib.print_split(cmd, outhandle=handle)
    print >>handle
    
    # Write out each module.
    for nid in node_order:
        inode = node_dict[nid]
        #node = inode.data

        # Input data type.
        if not hasattr(inode, "out_path"):
            #print >>handle, "Input: %s" % node.datatype.name
            #if inode.identifier:
            #    print >>handle, inode.identifier
            #print >>handle
            continue

        # Module.
        # module name.
        x = os.path.join(inode.out_path, rule_engine.BETSY_PARAMETER_FILE)
        params = rule_engine._read_parameter_file(x)
        metadata = params.get("metadata", {})
        
        module_name = params.get("module_name")
        x = "%s [Node %d]" % (module_name, nid)
        print >>handle, x
        print >>handle, "-"*len(x)
        # run time.
        assert "start_time" in params, "Missing: start_time"
        start_time = params["start_time"]
        run_time = params.get("elapsed_pretty")
        if run_time == "instant":
            x = "ran instantly"
        else:
            x = "took %s" % run_time
        print >>handle, "RUN: %s (%s)." % (start_time, x)
        #time_ = time.strptime(start_time, rule_engine.TIME_FMT)

        # input files.
        antecedents = params.get("antecedents")
        if antecedents:
            print >>handle, "INPUTS:"
            for i, (name, filename) in enumerate(antecedents):
                x = "%d.  %s" % (i+1, name)
                parselib.print_split(
                    x, prefix1=4, prefixn=8, outhandle=handle)
                if filename:
                    parselib.print_split(
                        filename, prefix1=8, prefixn=8, outhandle=handle)

        # output files.
        outfile = params["outfile"]
        outfilename = os.path.join(inode.out_path, outfile)
        size = parselib.pretty_filesize(mlib.get_dirsize(outfilename))
        x = "OUTPUT: %s (%s)" % (outfile, size)
        parselib.print_split(x, prefix1=0, prefixn=4, outhandle=handle)

        # working path.
        x = "WORKING PATH: %s" % inode.out_path
        parselib.print_split(x, outhandle=handle)

        # user options (--mattr)
        user_options = params.get("user_options")
        if user_options:
            print >>handle, "MODULE ATTRIBUTES:"
            for name in sorted(user_options):
                value = user_options[name]
                print >>handle, "    %s=%s" % (name, value)

        # commands
        commands = metadata.get("commands")
        if commands:
            if type(commands) is type(""):
                commands = [commands]
            for x in metadata["commands"]:
                print >>handle, x

        # miscellaneous metadata
        meta_lines = []
        for key, value in metadata.iteritems():
            if key == "commands":
                continue
            x = "%s=%s" % (key.upper(), value)
            meta_lines.append(x)
        if meta_lines:
            print >>handle, "METADATA:"
        for x in meta_lines:
            parselib.print_split(x, prefix1=4, prefixn=8, outhandle=handle)
        print >>handle


def get_all_option_names():
    # Return a list of the names for all OptionDef objects in the
    # rulebase.
    from Betsy import rulebase

    names = {}
    for module in rulebase.all_modules:
        for option in module.option_defs:
            names[option.name] = 1
    return sorted(names)


def get_module_options(modules, required_only=False):
    # From a list of modules, return a dictionary where the keys are
    # option name and the values are a list of the modules that
    # use it.  If required_only is True, then will only return the
    # options that are required.
    option2modules = {}
    for module in modules:
        for option in module.option_defs:
            is_required = option.default is None
            if required_only and not is_required:
                continue
            on, mn = option.name, module.name
            if on not in option2modules:
                option2modules[on] = []
            option2modules[on].append(mn)
    return option2modules


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


def _find_lowest_datatype(network, datatype_name, allowed_node_ids):
    # Return a node_id or None if not found.
    from Betsy import bie3

    x = allowed_node_ids
    x = [x for x in x if isinstance(network.nodes[x], bie3.DataNode)]
    x = [x for x in x if network.nodes[x].datatype.name == datatype_name]
    node_ids = x

    if not node_ids:
        return None
    if len(node_ids) == 1:
        return node_ids[0]
    descendents = bie3._make_descendent_dict(network)
    good_node_ids = []
    for nid in node_ids:
        x = descendents.get(nid, [])
        x = [x for x in x if isinstance(network.nodes[x], bie3.DataNode)]
        x = [x for x in x if network.nodes[x].datatype.name == datatype_name]
        x = [x for x in x if x in allowed_node_ids]
        if not x:
            good_node_ids.append(nid)
    if len(good_node_ids) == 1:
        return good_node_ids[0]
    return None


def _find_highest_datatype(network, datatype_name, allowed_node_ids):
    # Return a node_id or None if not found.
    from Betsy import bie3

    x = allowed_node_ids
    x = [x for x in x if isinstance(network.nodes[x], bie3.DataNode)]
    x = [x for x in x if network.nodes[x].datatype.name == datatype_name]
    node_ids = x

    if not node_ids:
        return None
    if len(node_ids) == 1:
        return node_ids[0]
    ancestors = bie3._make_ancestor_dict(network)
    good_node_ids = []
    for nid in node_ids:
        x = ancestors.get(nid, [])
        x = [x for x in x if isinstance(network.nodes[x], bie3.DataNode)]
        x = [x for x in x if network.nodes[x].datatype.name == datatype_name]
        x = [x for x in x if x in allowed_node_ids]
        if not x:
            good_node_ids.append(nid)
    if len(good_node_ids) == 1:
        return good_node_ids[0]
    return None


def _parse_dattr(dattr_str):
    # Format: <datatype>[*].<key>=<value>
    # Return <datatype>, <key>, <value>, <all_nodes>.

    #err_msg = "--dattr should be <datatype>.<key>=<value>[.<next_id>]"
    err_msg = "--dattr should be <datatype>[*].<key>=<value>"
    x = dattr_str.split(".", 1)
    assert len(x) == 2, err_msg
    #assert len(x) in [2, 3], err_msg
    datatype, kv = x[:2]
    #next_id = None
    #if len(x) == 3:
    #    next_id = int(x[2])
    x = kv.split("=", 1)
    assert len(x) == 2, err_msg
    key, value = x
    all_nodes = False
    if datatype.endswith("*"):
        datatype = datatype[:-1]
        all_nodes = True
    return datatype, key, value, all_nodes


def _parse_args(args):
    inputs = []            # list of names of datatypes
    in_identifiers = {}    # index (into inputs) -> identifier (e.g. filename)
    in_parameters = {}     # index (into inputs) -> list of (key, value)

    output = None          # name of datatype
    out_identifier = None  # identifier
    out_parameters = []    # list of (datatype, key, value, all_nodes)
    # out_parameters can refer to output, or other internal nodes in
    # the network.

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
            assert inputs
            assert len(args) >= i+1
            dattr = args[i+1]
            i += 2
            datatype, key, value, all_nodes = _parse_dattr(dattr)
            #assert next_id is None
            assert not all_nodes
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
            in_parameters[index].append((key, value))
        elif arg == "--dattr" and input_or_output == "--output":
            assert len(args) >= i+1
            dattr = args[i+1]
            i += 2
            datatype, key, value, all_nodes = _parse_dattr(dattr)
            assert output
            # This check won't work.  out_parameters also includes
            # attributes for internal nodes in the network that might
            # have different data types.
            #assert datatype == output, \
            #       "Datatype mismatch: --dattr %s and --output %s." % (
            #    datatype, output)
            out_parameters.append((datatype, key, value, all_nodes))
        elif arg == "--dattr":
            raise AssertionError, "--dattr before --input or --output"
        elif arg.startswith("--input_file"):
            # Possible formats:
            # 1.  --input_file fastq01
            # 2.  --input_file=fastq01
            assert input_or_output == "--input", \
                   "--input_file must be after --input and before --output"
            if arg == "--input_file":
                assert len(args) >= i+1
                filename = args[i+1]
                i += 2
            else:
                x = arg.split("=")
                assert len(x) == 2, "Invalid arg: %s" % arg
                assert x[0] == "--input_file"
                filename = x[1]
                i += 1
            index = len(inputs) - 1
            assert index >= 0
            assert index not in in_identifiers, \
                   "Multiple --input_file provided for %s" % inputs[-1]
            in_identifiers[index] = filename
        elif arg == '--output_file':
            assert input_or_output == "--output", \
                   "--output must precede --output_file"
            assert len(args) >= i+1
            filename = args[i+1]
            i += 2
            assert not out_identifier
            out_identifier = filename
        elif arg.startswith("--output_file"):
            x = arg.split("=")
            assert len(x) == 2, "Invalid arg: %s" % arg
            assert x[0] == "--output_file"
            filename = x[1]
            i += 1
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


def _make_custom_attr(rulebase, datatype, attributes, all_nodes):
    from Betsy import bie3
    
    assert attributes
    assert hasattr(rulebase, datatype), "Unknown datatype: %s" % datatype
    fn = getattr(rulebase, datatype)
    x = [bie3.Attribute(fn, name, value) for (name, value) in attributes]
    x = bie3.CustomAttributes(x, all_nodes)
    return x


def main():
    import os
    import sys
    import argparse
    import getpass
    #import itertools

    from Betsy import config
    from Betsy import rule_engine
    from Betsy import userfile
    from Betsy import reportlib
    from Betsy import rulebase
    from Betsy import bie3

    WORKFLOW = (
        " 1.  Look through the rulebase for the DataType of interest.\n"
        "     Run betsy_run.py (by itself)\n"
        " 2.  Specify the output DataType and its attributes.\n"
        "     Add argument: --output\n"
        "     Add argument(s): --dattr\n"
        " 3.  Browse through the list of possible combinations of input\n"
        "     DataTypes.\n"
        " 4.  Select one or more input DataTypes.\n"
        "     Now shows just the combinations that includes the input\n"
        "     DataTypes requested.\n"
        "     Add argument(s): --input\n"
        " 5.  Configure the attributes of the inputs.\n"
        "     Shows the detailed information about each node.\n"
        "     Add argument(s): --dattr\n"
        " 6.  System makes sure each node can be found in the network.\n"
        "     If not, try to diagnose differences.\n"
        " 7.  System makes sure this is a complete set of input nodes.\n"
        " 8.  System makes sure all required module attributes are given.\n"
        " 9.  System makes sure all input files are provided.\n"
        "10.  Visualize the pipeline to make sure this is what you want.\n"
        "     Plot the network: --network_png\n"
        "11.  Configure the attributes to alter the pipeline, if necessary.\n"
        "     Add argument(s): --dattr\n"
        "     Add argument(s): --mattr\n"
        "12.  Actually run the analysis.\n"
        "     Add argument: --run\n"
        )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Hi!  I'm BETSY, and I like to do bioinformatics.\n\n"
        "Workflow:\n%s\n\n" % WORKFLOW)
    #parser.add_argument(
    #    "--inputs_complete", action="store_const", const=True,
    #    help="Have finished specifying --input DataTypes.")
    #parser.add_argument(
    #    "--attrs_complete", action="store_const", const=True,
    #    help="Have finished configuring --dattr attributes.")
    parser.add_argument(
        "--run",  action="store_true", 
        help="Perform the analysis.")
    parser.add_argument(
        '--user', default=getpass.getuser(),
        help='the username who run the command')
    parser.add_argument(
        '--job_name', default='', help='the name of this job')
    parser.add_argument(
        '--num_cores', default=4, type=int,
        help='number of cores used in the processing')
    DEFAULT_MAX_INPUTS = 6
    parser.add_argument(
        '--max_inputs', default=DEFAULT_MAX_INPUTS, type=int,
        help="Maximum number of inputs to be shown (default %d).  "
        "(For optimization)" % DEFAULT_MAX_INPUTS)
    #parser.add_argument(
    #    "--save_failed_data", action="store_true",
    #    help="If a module failed, do not clean up its working files.")
    parser.add_argument(
        "--cache_input_files", action="store_true",
        help="Save a copy of the input files.  "
        "Helpful for very big files.")
    parser.add_argument(
        "-v", "--verbose", default=0, action="count",
        help="Give more output.")

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
    #group.add_argument(
    #    "--exclude_input", action="append", help="Experimental.")
    group.add_argument(
        "--dont_find_mattr_files", action="store_true",
        help="By default, if an --mattr option looks like a filename, "
        "will convert to full path so that modules can find them.  "
        "This will turn off that feature.")

    group.add_argument(
        '--output',  help='Desired DataType for the output.')
    group.add_argument(
        '--dattr', default=[], action='append',
        help="Specify a Datatype's attribute.  Since there may be multiple "
        "inputs or outputs with the same Datatype, the position of "
        "this argument is important.  --dattr should be given "
        "immediately following the Datatype that it refers to.  "
        "No --dattr should be given before the first --input.  "
        "--dattr that refer to internal nodes in the network should be "
        "given after the --output.  Format: <datatype>.<key>=<value>.")
    group.add_argument(
        '--output_file', help='file or folder of output result')
    group.add_argument(
        '--also_save_lowest', default=[], action="append",
        help="Will save the contents of other nodes.  "
        "Format: <datatype>,<filename>.  Will save the bottom-most "
        "node with this datatype.")
    group.add_argument(
        '--also_save_highest', default=[], action="append",
        help="Will save the contents of other nodes.  "
        "Format: <datatype>,<filename>.  Will save the top-most "
        "node with this datatype.")

    group = parser.add_argument_group(title="Output")
    group.add_argument(
        '--network_png',
        help='Generate a PNG that shows the data flow graph.')
    group.add_argument(
        '--sparse_network_png', action="store_true",
        help="Leave out details in network plot.")
    group.add_argument(
        '--prune_network', action="store_true",
        help="Prune nodes that are not included in any pipeline.  "
        "Mostly for debugging.")
    #group.add_argument(
    #    '--network_text', help='generate the output network text file')
    group.add_argument(
        '--network_json', help='generate the output network json file')
    group.add_argument(
        '--restart_from_network', action="store_true",
        help="If the --network_json file already exists, "
        "then will use this network rather than recreating a new one.")
    group.add_argument(
        "--receipt",
        help="Name of file to write a receipt (as a text file) "
        "describing the analysis.  ")
    #parser.add_argument(
    #    '--clobber', action='store_const', const=True, default=False,
    #    help='overwrite the output_data if it already exists')

    print "Starting rule engine."
    sys.stdout.flush()

    # Parse the arguments.
    args = parser.parse_args()
    input_list, x = _parse_args(sys.argv)
    outtype, out_identifier, out_attributes = x
    verbose_network = (not args.sparse_network_png)

    print "Checking parameters."
    sys.stdout.flush()

    args.clobber = True
    assert args.num_cores > 0, "num_cores should be greater than 0"
    assert args.max_inputs > 0 and args.max_inputs <= 20
    args.save_failed_data = True   # Always save now.
    #if args.inputs_complete or args.attrs_complete or args.run:
    #    if args.run:
    #        x = "--run"
    #    elif args.attrs_complete:
    #        x = "--attrs_complete"
    #    else:
    #        x = "--inputs_complete"
    #    assert input_list, "%s given, but no --input." % x

    #if args.run:
    #    x = "--run"
    #    assert input_list, "%s given, but no --input." % x
    ## TODO: Make sure args.exclude_input is valid.
    ## args.exclude_input

    if args.restart_from_network:
        assert args.network_json, \
               "Cannot restart_from_network: no --network_json specified."

    # Make sure configuration directory exists.
    assert hasattr(config, "CACHE_PATH"), \
           "Not defined in .genomicoderc: CACHE_PATH"
    if not os.path.exists(config.CACHE_PATH):
        print "Making BETSY working path: %s." % config.CACHE_PATH
        os.mkdir(config.CACHE_PATH)

    # Make sure files exist.
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
    # List of identified data nodes.
    in_data_nodes = []
    for x in input_list:
        intype, identifier, attributes = x
        assert hasattr(rulebase, intype), "Unknown datatype: %s" % intype
        fn = getattr(rulebase, intype)
        params = {}
        for key, value in attributes:
            params[key] = value
        # Make sure the identifier is a full path since we'll be
        # changing directories.
        if identifier:
            identifier = os.path.realpath(identifier)
        in_data = fn.input(**params)
        x = bie3.IdentifiedDataNode(in_data, identifier)
        in_data_nodes.append(x)

    # test outtype and build the list of custom_attributes.
    custom_attributes = []  # List of bie3.CustomAttributes objects.
    out_custom_attribute = None  # bie3.CustomAttributes
    if outtype:
        # Pull out the custom_attributes.  These attributes can refer
        # to:
        # 1.  The input data (to guide the inferencing engine).
        #     e.g. Input has indel_realigned=yes, so downstream nodes
        #     should also have indel_realigned=yes.
        # 2.  The output data (out_attributes).
        # 3.  Internal nodes (out_attributes, different data type).

        # Case 1.  Get the custom attributes for the input nodes.
        attrs = []  # list of (datatype, list of (key, value), all_nodes)
        for x in input_list:
            intype, identifier, attributes = x
            if not attributes:
                continue
            x = intype, attributes, False
            attrs.append(x)
        for x in attrs:
            datatype, attributes, all_nodes = x
            x = _make_custom_attr(rulebase, datatype, attributes, all_nodes)
            custom_attributes.append(x)

        # Case 2 and 3.  out_attributes contains attributes for both
        # output and internal nodes.  Group them by datatype.
        # out_attributes is a list of (datatype, key, value, all_nodes)
        x = [x[0] for x in out_attributes]
        x = sorted({}.fromkeys(x))
        datatypes = x
        
        # datatypes can be empty if no attributes specified.
        for dt in datatypes:
            # Get all_nodes for this data type.  Should all be the same.
            x = [x[3] for x in out_attributes if x[0] == dt]  # all_nodes
            x = sorted({}.fromkeys(x))
            assert len(x) == 1, "Inconsistent all_nodes: %s" % dt
            all_nodes = x[0]
            out_attrs = [x for x in out_attributes
                         if x[0] == dt and x[3] == all_nodes]
            assert out_attrs
            attributes = [(x[1], x[2]) for x in out_attrs]
            # Since there's only one output object, all attributes for
            # that data type should be put into the same
            # CustomAttributes object.
            if dt == outtype:
                x = _make_custom_attr(rulebase, dt, attributes, all_nodes)
                out_custom_attribute = x
            # Otherwise, make separate CustomAttribute objects.
            else:
                for x in attributes:
                    x = _make_custom_attr(rulebase, dt, [x], all_nodes)
                    custom_attributes.append(x)
                #x = _make_custom_attr(
                #    rulebase, dt, attributes, all_nodes)
                #custom_attributes.append(x)
                
        
        ## # out_attributes contains attributes for the output node, as
        ## # well as other nodes in the network. Just group them by
        ## # datatype.
        ## # Since there's only one output, it should be put into the
        ## # same CustomAttributes object.
        ## # out_attributes is a list of (datatype, key, value, all_nodes)
        ## if out_attributes:
        ##     x1 = [x[0] for x in out_attributes]  # datatype
        ##     x2 = [x[3] for x in out_attributes]  # all_nodes
        ##     x1 = sorted({}.fromkeys(x1))
        ##     x2 = sorted({}.fromkeys(x2))
        ##     for x in itertools.product(x1, x2):
        ##         datatype, all_nodes = x
        ##         out_attrs = [x for x in out_attributes
        ##                  if x[0] == datatype and x[3] == all_nodes]
        ##         if not out_attrs:
        ##             continue
        ##         attributes = [(x[1], x[2]) for x in out_attrs]
        ##         x = datatype, attributes, all_nodes
        ##         attrs.append(x)
        
    # Cache the files in the user's directory.  Don't depend on the
    # user keeping the file in the same place.  Needed for the first
    # module.
    if args.cache_input_files:
        print "Making a local copy of the input files."
        sys.stdout.flush()
        for i_data_node in in_data_nodes:
            filename = i_data_node.identifier
            if not filename or not os.path.exists(filename):
                continue
            x = userfile.set(getpass.getuser(), filename)
            i_data_node.identifier = x

    # Parse out the module attributes (AKA user_options).
    user_options = {}
    if args.mattr:
        #print "Parsing module attributes."
        all_mattrs = get_all_option_names()
    for x in args.mattr:
        assert '=' in x, "--mattr should be in format: <option>=<value>\n%s" %x
        key, value = x.split('=', 1)
        assert key in all_mattrs, "Unknown module attribute: %s" % key
        user_options[key] = value

    # Since the modules will be run from their own directory, any
    # module attributes that point to files (as relative paths) will
    # be messed up.  See if any of these look like files and convert
    # them into full path names.
    if not args.dont_find_mattr_files:
        for key, value in user_options.iteritems():
            if not os.path.exists(value):
                continue
            x = os.path.realpath(value)
            user_options[key] = x
    
    # Step 1: Make sure there's an output provided.
    #         If not, show the rulebase.
    # Step 2: Generate network that produces this output.
    # Step 3: Make sure there are inputs provided.
    #         If not, show the input datatypes.
    # Step 4: Make sure the inputs can be found in the network.
    # Step 5: Create the pipelines.
    # Step 6: Make sure required attributes are given.
    # Step 7: Prune redundant pipelines.
    # Step 8: Make sure input files exist.
    # Step 9: Manual verification of the network by the user.

    # Step 1: Make sure an output is provided.
    if not check_output_provided(rulebase, args.output):
        return
    # Step 2: Generate network.
    if args.network_json and os.path.exists(args.network_json) and \
       args.restart_from_network:
        network = bie3.read_network(args.network_json)
    else:
        network = generate_network(
            rulebase, outtype, custom_attributes, out_custom_attribute)
        if args.network_json:
            write_network(args.network_json, network)
    #plot_network(
    #    "network.png", network, user_options=user_options,
    #    verbose=verbose_network)

    # Step 2.5: Make sure network has more than one node.
    if not check_more_than_one_node_network(network, rulebase):
        return
    
    # Step 3: Make sure some inputs are provided.
    if not check_inputs_provided(
        network, in_data_nodes, custom_attributes, user_options,
        args.max_inputs, args.network_png, verbose_network):
        return
    # Step 4: Make sure each of the input nodes match a node in the
    # network.
    x = check_inputs_in_network(
        network, user_options, in_data_nodes, 
        args.network_png, verbose_network)
    inputs_ok, data_node_ids = x
    if not inputs_ok:
        return
    # Step 5: Search for pipelines that can be run given the INPUT
    # nodes.
    paths = build_pipelines(
        network, user_options, in_data_nodes, data_node_ids, custom_attributes,
        args.max_inputs, args.network_png, verbose_network)
    #plot_pipelines(
    #    "pipeline", network, paths, user_options, max_pipelines=16,
    #    prune=True, verbose=True)

    if not paths:
        return
    # Step 6: Make sure required attributes are given.
    paths = check_attributes_complete(
        network, user_options, paths, 
        args.network_png, args.prune_network, verbose_network)
    if not paths:
        return
    # DEBUG: To debug the pruning, save the network and paths so we
    # don't have to re-generate them each time we run.
    #import pickle
    #open("network.txt", 'w').write(pickle.dumps(network))
    #open("paths.txt", 'w').write(pickle.dumps(paths))
    #import sys; sys.exit(0)
    #network = pickle.loads(open("network.txt").read())
    #paths = pickle.loads(open("paths.txt").read())
    # DEBUG: Print out each of the pipelines.
    #plot_pipelines(
    #    "pipeline", network, paths, user_options, max_pipelines=16,
    #    verbose=True)
    #plot_network_show_pipelines(
    #    args.network_png, network, paths, user_options=user_options,
    #    verbose=verbose_network)
    
    # Step 7: Prune undesired pipelines.
    paths = prune_pipelines(
        network, user_options, custom_attributes, paths, 
        args.network_png, args.prune_network, verbose_network)
    if not paths:
        return
    # DEBUG: Print out each of the pipelines.
    #plot_pipelines(
    #    "pipeline", network, paths, user_options, max_pipelines=16,
    #    verbose=True)
        
    # Step 8: Look for input files.
    if not check_input_files(
        network, in_data_nodes, user_options, paths, 
        args.network_png, args.prune_network, verbose_network):
        return

    # Step 9: Look for output file.
    if not check_output_file(args.output_file):
        return

    #print "The network (%d nodes) is complete." % len(network.nodes)
    #print "There are %d possible pipelines." % len(paths)
    x = "core"
    if args.num_cores > 1:
        x = "cores"
    print "Ready to go!  Will run the analysis using a maximum of %d %s." % (
        args.num_cores, x)
    # Step 10: Manual verification of the network.
    if not manually_verify_network(
        network, user_options, paths, args.run, 
        args.network_png, args.prune_network, verbose_network):
        return

    #plot_network_show_pipelines(
    #    args.network_png, network, paths, user_options=user_options,
    #    verbose=args.verbose)

    print "Running the analysis."
    sys.stdout.flush()
    clean_up = not args.save_failed_data
    node_dict = transitions = None
    try:
        x = rule_engine.run_pipeline(
            network, in_data_nodes, custom_attributes, user_options, paths,
            user=args.user, job_name=args.job_name, clean_up=clean_up,
            num_cores=args.num_cores, verbosity=args.verbose)
        if x:
            node_dict, transitions = x
    except AssertionError, x:
        if str(x).startswith("Inference error"):
            node_ids = rule_engine.DEBUG_POOL.keys()
            plot_network(
                args.network_png, network, user_options=user_options,
                highlight_green=node_ids, verbose=verbose_network)
            #write_network(args.network_json, network) 
        raise

    if args.output_file and node_dict and 0 in node_dict:
        print "Saving output %s to %s." % (outtype, args.output_file)
        sys.stdout.flush()
        output_file = node_dict[0].identifier
        reportlib.copy_file_or_path(output_file, args.output_file)

    # See what else to save.
    also_save = []  # list of (arg, "lowest" or "highest")
    for arg in args.also_save_lowest:
        also_save.append((arg, "lowest"))
    for arg in args.also_save_highest:
        also_save.append((arg, "highest"))
    for arg, which_one in also_save:
        # Format: <datatype>,<file_or_path>
        x = arg.split(",", 1)
        assert len(x) == 2, "Invalid also_save: %s" % arg
        dname, out_filename = x
        if which_one == "lowest":
            node_id = _find_lowest_datatype(network, dname, node_dict.keys())
        else:
            node_id = _find_highest_datatype(network, dname, node_dict.keys())
        assert node_id, "Unable to find: %s" % dname
        in_filename = node_dict[node_id].identifier
        print "Saving %s to %s." % (dname, out_filename)
        reportlib.copy_file_or_path(in_filename, out_filename)
        
    # Draw the final network.
    node_dict = node_dict or {}
    transitions = transitions or {}
    start_ids = []
    for x in data_node_ids:
        x = [x for x in x if x in node_dict]
        start_ids.extend(x)
    node_ids = [x for x in node_dict if x not in start_ids]
    # node_dict only contains DataNodes.  Add the ModuleNodes where
    # there is a transition into and out of.
    for nid in range(len(network.nodes)):
        if not isinstance(network.nodes[nid], bie3.ModuleNode):
            continue
        if nid in node_ids:
            continue
        found1 = found2 = False
        for (n1, n2) in transitions:
            if nid == n1:
                found1 = True
            elif nid == n2:
                found2 = True
            if found1 and found2:
                break
        if found1 and found2:
            node_ids.append(nid)
    plot_network(
        args.network_png, network,
        user_options=user_options,
        bold=start_ids+node_ids,
        bold_transitions=transitions,
        highlight_green=start_ids, highlight_yellow=node_ids,
        verbose=verbose_network)
    if args.receipt:
        write_receipt(args.receipt, network, start_ids, transitions, node_dict)
    #if args.network_text:
    #    print "Writing detailed network: %s." % args.network_text
    #    bie3.print_network(network, outhandle=args.network_text)

    print "Done."


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
