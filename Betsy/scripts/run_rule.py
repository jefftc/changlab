#! /usr/bin/env python

# Functions:
# check_output_provided      Step 1: Make sure there's an output.
#   _print_rulebase          To help select output DataType.
#   _pretty_print_datatype
#   _pretty_print_module
# generate_network           Step 2.
# check_inputs_provided      Step 3: Make sure any --inputs are provided.
#   _print_input_datatypes   To help select input DataType(s).
# check_inputs_in_network    Step 4.
#   _print_node_score_table
# build_pipelines            Step 5.
# check_attributes_complete  Step 6.
# prune_pipelines            Step 7.
#   _prune_alternate_attributes
#   _find_alternate_attributes
#   _prune_superset_pipelines
#   _is_superset_pipeline
#   _does_path_go_through
#   _prune_parallel_pipelines
#   _is_parallel_pipeline
# check_input_files          Step 8.
# manually_verify_network    Step 9.
# 
# plot_network
# plot_network_show_pipelines
# plot_pipelines             
# 
# get_all_option_names       Get the names of all OptionDef in the rulebase.
# get_required_option_names  Get the required options in a list of modules.
#
# _list_differences_in_nodelists
# _merge_nodelists
# _parse_dattr
# _parse_args                 Parse the arguments from the user.


def check_output_provided(rulebase, output_file):
    if output_file:
        return True
    # If there is no output, just print the whole rulebase.
    print "Selecting --output DataType."
    _print_rulebase(rulebase)
    return False


def _print_rulebase(rulebase):
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


def generate_network(rulebase, outtype, custom_attributes):
    import sys
    from Betsy import bie3

    assert outtype
    a_or_an = "a"
    if outtype.lower()[0] in "aeiou":
        a_or_an = "an"
    print "Generating a network that produces %s %s." % (a_or_an, outtype)
    sys.stdout.flush()

    # Get the out_datatype.
    assert hasattr(rulebase, outtype), "Unknown datatype: %s" % outtype
    out_datatype = getattr(rulebase, outtype)
    attrs = {}
    for attr in custom_attributes:
        if attr.datatype.name != out_datatype.name:
            continue
        attrs[attr.name] = attr.value
    out_data = out_datatype.output(**attrs)

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
    print "Made a network with %d nodes." % len(network.nodes)
    return network


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
    x2 = "combinations"
    if len(datatype_combos) == 1:
        x1 = "This is"
        x2 = "combination"
    x3 = ""
    if max_inputs is not None:
        x3 += " (up to %d)" % max_inputs
    x4 = "%s the acceptable %s of --input DataTypes%s.  " % (x1, x2, x3)
    x5 = "The exact number of each DataType needed is not shown."
    x = x4 + x5
    ps(x, prefixn=0, outhandle=outhandle)
    for i, combo in enumerate(datatype_combos):
        x = [x.name for x in combo]
        ps("%3d.  %s" % (i+1, ", ".join(x)), prefixn=6, outhandle=outhandle)


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
    # Score cutoff is minimum score + 2
    score_cutoffs = {}  # dtname -> max score to print
    for dtname, score in dt2scores.iteritems():
        score_cutoffs[dtname] = min(score) + 2

    
    # Make an output table.
    table = []
    header = [
        "Node", "S", "Datatype", "Attribute", "Yours", "Network"]
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
    from Betsy import bie3

    print "Building pipelines that use --input data types."
    sys.stdout.flush()

    # data_node_ids is parallel to in_data_nodes.  Each element is a
    # list of the node_ids that that data node can map onto.

    paths = bie3.find_paths_by_start_ids(
        network, custom_attributes, data_node_ids)
    
    # Make sure all the --inputs are needed.  Any unnecessary ones may
    # indicate a problem.
    inputs_used = {}  # list of indexes of --inputs that are used
    for p in paths:
        used = [i for (i, x) in enumerate(p.start_ids) if x is not None]
        inputs_used.update({}.fromkeys(used))
    has_unused_inputs = len(inputs_used) != len(in_data_nodes)

    good_paths = [x for x in paths if not x.missing_ids]
    if good_paths and not has_unused_inputs:
        print "Found %d possible pipelines." % len(good_paths)
        return good_paths

    # Print out --inputs that aren't used.
    if has_unused_inputs:
        for i in range(len(in_data_nodes)):
            if i in inputs_used:
                continue
            name = in_data_nodes[i].data.datatype.name
            print (
                "%s is not used in any pipelines.  "
                "Please make sure the proposed pipelines are acceptable and "
                "remove this input." % name)
            
    if not good_paths:
        print "No pipelines found.  Examine network to diagnose."
        print "Make sure that no --input is missing."
        x = [x.data.datatype.name for x in in_data_nodes]
        _print_input_datatypes(
            network, x, custom_attributes, max_inputs=max_inputs)
        print

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
    network, user_options, paths, network_png, verbose):
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
    for p in paths:
        module_ids = [
            x for x in p.node_ids if
            isinstance(network.nodes[x], bie3.ModuleNode)]
        modules = [network.nodes[x] for x in module_ids]

        missing = {}
        opt2mods = get_required_option_names(modules)
        x = [x for x in opt2mods if x not in user_options]
        for on in x:
            for mn in opt2mods[on]:
                missing[(mn, on)] = 1
        extra = [x for x in user_options if x not in opt2mods]
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
                x = ("Removed %d pathways because --mattr option %r was "
                     "not provided." % (num_removed, on))
            else:
                names = sorted([x[1] for x in all_missing.keys()])
                x = parselib.pretty_list(names)
                x = ("Removed %d pathways because the following --mattr "
                     "options were not provided: %s" % (num_removed, x))
            parselib.print_split(x, prefixn=2)
        return good_paths
    for (mn, on) in all_missing:
        print 'Missing --mattr: %s requires attribute "%s".' % (mn, on)

    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        verbose=verbose)
    return []


def prune_pipelines(
    network, user_options, custom_attributes, paths, network_png, verbose):
    # Any any pipelines look weird, then remove them.
    import sys
    from Betsy import bie3

    print "Pruning redundant pipelines."
    sys.stdout.flush()
    paths_orig = paths[:]
    

    nodeid2parents = bie3._make_parents_dict(network)
    
    # Optimizations:
    # 1.  Convert node_ids to dictionaries.
    # 2.  Convert the transitions to sets for easier comparison.
    for i, p in enumerate(paths):
        # Not much more efficient.
        #p.node_ids_set = set(p.node_ids)
        p.node_ids = {}.fromkeys(p.node_ids)
        for key, value in p.transitions.iteritems():
            p.transitions[key] = set(value)

    # Do the O(N) pruning, then fast O(NN) pruning, then slow O(NN)
    # pruning.
    # Just try different orders until I find the fastest.
    paths = _prune_by_custom_attributes(
        network, custom_attributes, paths, nodeid2parents)
    paths = _prune_alternate_attributes2(
        network, custom_attributes, paths, nodeid2parents)
    paths = _prune_alternate_attributes1(
        network, custom_attributes, paths, nodeid2parents)
    paths = _prune_superset_pipelines(network, paths)
    paths = _prune_parallel_pipelines(network, paths, nodeid2parents)
    
    # Convert node_ids back to lists.
    for i, p in enumerate(paths):
        p.node_ids = p.node_ids.keys()
        for key, value in p.transitions.iteritems():
            p.transitions[key] = list(value)
            
    if not paths:
        print "All pipelines pruned.  This can happen if:"
        print "  1.  A --mattr option is missing."
        print "  2.  There is a bug in the network generation."
        print "  3.  There is a bug in the pipeline pruning."
        print "Please review the network and --mattr options."
        plot_network_show_pipelines(
            network_png, network, paths_orig, user_options=user_options,
            verbose=verbose)
        return paths
    
    num_pruned = len(paths_orig) - len(paths)
    if not num_pruned:
        print "No redundant pipelines found.  %d left." % len(paths)
    else:
        x = ""
        if num_pruned >= 2:
            x = "s"
        print "Pruned %d pipeline%s.  %d left." % (num_pruned, x, len(paths))
    return paths


def _prune_by_custom_attributes(network, custom_attributes, paths,
                                nodeid2parents):
    # Keep only the paths that match the attributes desired by the
    # user.  custom_attributes is a list of Attribute objects.
    from Betsy import bie3
    if not custom_attributes:
        return paths
    
    dname2attrs = {}  # datatype name -> attr name -> value
    for x in custom_attributes:
        if x.datatype.name not in dname2attrs:
            dname2attrs[x.datatype.name] = {}
        dname2attrs[x.datatype.name][x.name] = x.value

    # If a module's input datatype is different from the output
    # datatype, and it has has no descendents with the same datatype,
    # then the attributes of that input datatype should match the user
    # attributes.


    # Search through the network for data nodes that might be subject
    # to custom_attributes.
    all_node_ids = {}
    for x in paths:
        all_node_ids.update(x.node_ids)

    descendents = bie3._make_descendent_dict(network)
    ## Only care about the descendents that are in this network.
    #desc = {}
    #for (node_id, next_ids) in descendents.iteritems():
    #    if node_id not in all_node_ids:
    #        continue
    #    next_ids = [x for x in next_ids if x in all_node_ids]
    #    desc[node_id] = next_ids
    #descendents = desc
    
    module_ids = [
        x for x in all_node_ids
        if isinstance(network.nodes[x], bie3.ModuleNode)]
    data_node_ids = []
    for module_id in module_ids:
        module = network.nodes[module_id]
        
        # If the module takes only one kind of datatype and produces
        # the same kind, then this won't work.
        if len(module.in_datatypes) == 1 and \
           module.out_datatype == module.in_datatypes[0]:
            continue

        # Make sure at least one input datatype (that is not
        # DefaultAttributesFrom) has a user attribute.
        daf = [x.input_index for x in module.default_attributes_from]
        in_dtype_names = [
            x.name for (i, x) in enumerate(module.in_datatypes)
            if i not in daf]
        x = [x for x in in_dtype_names if x in dname2attrs]
        if not x:
            continue
        # This should be a converting module.

        # Find the node_ids that match a user attribute.
        x = nodeid2parents.get(module_id, [])
        x = [x for x in x if network.nodes[x].datatype.name in dname2attrs]
        node_ids = x

        # Make sure these node_ids don't have any descendents of the
        # same type.
        good_ids = []
        for node_id in node_ids:
            x = descendents.get(node_id, [])
            x = [x for x in x if x in all_node_ids]
            x = [x for x in x if isinstance(network.nodes[x], bie3.DataNode)]
            x = [network.nodes[x].datatype.name for x in x]
            if network.nodes[node_id].datatype.name not in x:
                good_ids.append(node_id)
        node_ids = good_ids

        data_node_ids.extend(node_ids)
    data_node_ids = {}.fromkeys(data_node_ids)

    # List the node_ids that either don't match the user attributes,
    # or are ambiguous, e.g.
    # Fastq.adapters = [no, yes]; user wants no
    mismatch_node_ids = {}
    ambiguous_node_ids = {}
    for node_id in data_node_ids:
        node = network.nodes[node_id]
        attrs = dname2attrs[node.datatype.name]
        for name, uvalue in attrs.iteritems():
            dvalue = node.attributes[name]
            utype = bie3._get_attribute_type(uvalue)
            dtype = bie3._get_attribute_type(dvalue)
            assert utype == bie3.TYPE_ATOM
            if dtype == bie3.TYPE_ENUM:
                ambiguous_node_ids[node_id] = 1
                continue
            if dvalue != uvalue:
                mismatch_node_ids[node_id] = 1

    if not mismatch_node_ids and not ambiguous_node_ids:
        return paths

    # For each pathway, see if it contains a delete or ambiguous node.
    bc_cache = {}
    fc_cache = {}
    path_cache = {}
    
    delete = {}
    for (i, p) in enumerate(paths):
        #node_ids, start_ids, data_indexes = x
        # If there's a mismatch, delete this path.
        x = [x for x in mismatch_node_ids if x in p.node_ids]
        if x:
            delete[i] = 1
            continue

        # If there's no ambiguity, then keep this path.
        x = [x for x in ambiguous_node_ids if x in p.node_ids]
        if not x:
            continue
        check_ids = x

        # See if this pathway can produce a node that has a proper
        # attribute.
        check_paths = []
        for out_data_id in check_ids:
            x = _bc_to_input_and_module_ids(
                network, out_data_id, custom_attributes, p.node_ids,
                nodeid2parents, bc_cache)
            check_paths.extend(x)
        assert check_paths

        # Each node ID in check_ids needs to match.
        match = {}  # node_id -> 1
        for x in check_paths:
            in_data_ids, module_id, out_data_id = x
            if out_data_id in match:
                continue
            name = network.nodes[out_data_id].datatype.name
            attrs = dname2attrs[name]
            key = module_id, in_data_ids, out_data_id
            if key not in path_cache:
                x = _is_valid_output_from_input_and_module_ids(
                    network, module_id, in_data_ids, out_data_id,
                    attrs, fc_cache)
                path_cache[key] = x
            if path_cache[key]:
                match[out_data_id] = 1
            if len(match) == len(check_ids):
                break
        if len(match) != len(check_ids):
            delete[i] = 1

    paths = [x for (i, x) in enumerate(paths) if i not in delete]
    return paths


def _bc_to_input_and_module_ids(
    network, out_id, custom_attributes, allowed_ids, nodeid2parents, cache):
    from Betsy import bie3
    
    paths = []
    module_ids = nodeid2parents.get(out_id, [])
    module_ids = [x for x in module_ids if x in allowed_ids]
    for module_id in module_ids:
        key = module_id, out_id
        if key not in cache:
            x = bie3._bc_to_input_ids(
                network, module_id, custom_attributes, all_output_ids=[out_id],
                nodeid2parents=nodeid2parents)
            cache[key] = x
        combos = cache[key]
        for in_ids in combos:
            assert module_id in allowed_ids
            assert out_id in allowed_ids
            x = [x for x in in_ids if x in allowed_ids]
            if len(x) != len(in_ids):
                continue
            x = in_ids, module_id, out_id
            paths.append(x)
    return paths


def _is_valid_output_from_input_and_module_ids(
    network, module_id, in_data_ids, out_data_id, user_attrs, cache):
    from Betsy import bie3

    in_datas = [network.nodes[x] for x in in_data_ids]
    module = network.nodes[module_id]

    key = module_id, in_data_ids
    if key not in cache:
        x = bie3._fc_to_outputs(module, in_datas)
        cache[key] = x
    out_datas = cache[key]

    name = network.nodes[out_data_id].datatype.name
    for out_data in out_datas:
        assert out_data.datatype.name == name
        for name, uvalue in user_attrs.iteritems():
            dvalue = out_data.attributes[name]
            utype = bie3._get_attribute_type(uvalue)
            dtype = bie3._get_attribute_type(dvalue)
            assert utype == bie3.TYPE_ATOM
            assert dtype == bie3.TYPE_ATOM
            if uvalue != dvalue:
                return False
    return True


def _prune_alternate_attributes1(
    network, custom_attributes, paths, nodeid2parents):
    # If there is an object with an attribute that can take several
    # different values, then pipelines may be generated that can
    # create each value.  This causes unnecessary computation, because
    # only one value needs to be used.  Arbitrarily choose one, and
    # delete the other pipelines.
    # 
    # align_with_bowtie1 -> SamFolder (bowtie1)
    # align_with_bowtie2 -> SamFolder (bowtie2)
    # align_with_bwa_aln -> SamFolder (bwa_backtrack)
    # align_with_bwa_mem -> SamFolder (bwa_mem)
    from Betsy import bie3
    #from genomicode import parselib

    # 1.  Look in two pipelines for:
    #     - Same DataObject.
    #     - Different Modules pointing to that DataObject.
    #     - Modules have Consequence SET_TO to different values in the
    #       same attribute.
    #     - That attribute is TYPE_ENUM in the DataObject.
    # 2.  Keep the pipeline with the Module that generates the value
    #     that comes first alphabetically.
    # 3.  Print a message showing which value was chosen.

    # List of (datatype name, attr name, kept value, deleted value)
    deleted = []
    paths = paths[:]

    # For each data node, list the alternate attributes.
    dataid2alts = {}
    for i in range(len(network.nodes)):
        if not isinstance(network.nodes[i], bie3.DataNode):
            continue
        x = _list_alternate_attributes1(network, i, nodeid2parents)
        dataid2alts[i] = x

    # Only check pathways that has a data node in dataid2alts.
    has_alts = []
    for i, p in enumerate(paths):
        #node_ids, start_ids, data_indexes = x
        x = [x for x in p.node_ids if dataid2alts.get(x, [])]
        if x:
            has_alts.append(i)
    has_alts = {}.fromkeys(has_alts)

    ## Convert node_ids to bitwise numbers for faster superset
    ## comparison.
    #paths_node_ids_bit = [_intlist2bits(x[0]) for x in paths]

    seen = {}  # (path_i, path_j) -> 1
    while True:
        found = None
        # Actually more efficient to restart this loop after deleting
        # paths than to run through it without deleting anything.
        # After deleting paths, leads to overall fewer calls to
        # _find_alternate_attributes.  O(N*N) growth.
        for i in range(len(paths)-1):
            if not has_alts.get(i):
                continue
            for j in range(i+1, len(paths)):
                if not has_alts.get(j):
                    continue
                if (i, j) in seen:
                    continue
                # Almost all time spent in this function.
                x = _find_alternate_attributes(
                    network, paths[i], paths[j], dataid2alts)
                seen[(i, j)] = 1
                if x:
                    found = i, j, x
                    break
            if found:
                break
        if not found:
            break

        path_1, path_2, x = found
        attr_name, module_id_1, attr_value_1, \
                   module_id_2, attr_value_2, data_id = x
        assert attr_value_1 != attr_value_2
        datatype_name = network.nodes[data_id].datatype.name

        # Figure out which path to delete.
        path_to_delete = None
        # Look in custom_attributes to see if one of them is requested
        # by the user.
        user_request = None  # value requested by user
        for x in custom_attributes:
            if x.datatype.name != datatype_name:
                continue
            if x.name != attr_name:
                continue
            assert user_request is None, "Multiple custom_attributes."
            user_request = x.value
        deleted_based_on_user = True
        if user_request:
            if user_request == attr_value_2:
                path_to_delete = 1
            elif user_request == attr_value_1:
                path_to_delete = 2
        if path_to_delete is None:
            # Keep the one that comes first alphabetically.
            path_to_delete = 1
            if attr_value_1 < attr_value_2:
                path_to_delete = 2
            deleted_based_on_user = False

        # Delete the appropriate path.
        assert path_to_delete in [1, 2]
        if path_to_delete == 1:
            del paths[path_1]
            #del paths_node_ids_bit[path_1]
            p1, p2 = path_2, path_1
            v1, v2 = attr_value_2, attr_value_1
        else:
            del paths[path_2]
            #del paths_node_ids_bit[path_2]
            p1, p2 = path_1, path_2
            v1, v2 = attr_value_1, attr_value_2

        # Fix the indexes of the seen variable.
        x = seen.keys()
        x = bie3._fix_node_id_pairs_after_merge(x, p1, p2)
        seen = {}.fromkeys(x)
            
        x = datatype_name, attr_name, v1, v2, deleted_based_on_user
        deleted.append(x)

    _print_attributes_pruned(deleted)

    return paths


def _list_alternate_attributes1(network, data_id, nodeid2parents):
    # Return list of (data_id, parent_ids, attr_name, attr_values).
    # align_with_bowtie1 -> SamFolder (bowtie1)
    # align_with_bowtie2 -> SamFolder (bowtie2)
    # align_with_bwa_aln -> SamFolder (bwa_backtrack)
    # align_with_bwa_mem -> SamFolder (bwa_mem)
    from Betsy import bie3

    if data_id not in nodeid2parents:
        return []
    prev_ids = nodeid2parents[data_id]

    # 2 parents must have Consequences SET_TO to the same
    # attribute.
    attr2values = {}     # name -> list of values to set to
    attr2moduleids = {}  # name -> list of module_ids
    for prev_id in prev_ids:
        module_node = network.nodes[prev_id]
        for cons in module_node.consequences:
            if cons.behavior != bie3.SET_TO:
                continue
            n, v = cons.name, cons.arg1
            assert v is not None
            if n not in attr2values:
                attr2values[n] = []
                attr2moduleids[n] = []
            if v not in attr2values[n]:
                attr2values[n].append(v)
                attr2moduleids[n].append(prev_id)
    # List of attributes with at least 2 values.
    attrs = [n for (n, v) in attr2values.iteritems() if len(v) >= 2]
    retvals = []
    for attr_name in attrs:
        x = data_id, attr2moduleids, attr_name, attr2values[attr_name]
        retvals.append(x)
    return retvals


def _find_alternate_attributes(network, path_1, path_2, dataid2alts):
    from Betsy import bie3
    
    shared_ids = [x for x in path_1.node_ids if x in path_2.node_ids]
    unique_ids_1 = [x for x in path_1.node_ids if x not in path_2.node_ids]
    unique_ids_2 = [x for x in path_2.node_ids if x not in path_1.node_ids]

    # Look for a DataNode with alternates that is shared in both pathways.
    data_ids = [x for x in shared_ids if dataid2alts.get(x, [])]
    if not data_ids:
        return None

    # Look for ModuleNodes that are unique for each pathway.
    module_ids_1 = [x for x in unique_ids_1
                    if isinstance(network.nodes[x], bie3.ModuleNode)]
    module_ids_2 = [x for x in unique_ids_2
                    if isinstance(network.nodes[x], bie3.ModuleNode)]
    # This is slower.
    #module_ids_1 = [x for x in node_ids_1 if x not in shared_ids and 
    #                isinstance(network.nodes[x], bie3.ModuleNode)]
    #module_ids_2 = [x for x in node_ids_2 if x not in shared_ids and 
    #                isinstance(network.nodes[x], bie3.ModuleNode)]
    if not module_ids_1 or not module_ids_2:
        return None

    # Check each data_id carefully.
    for data_id in data_ids:
        for x in dataid2alts[data_id]:
            x, parent_ids, attr_name, attr_values = x

            # Look for parent_ids that are in different pathways.
            good_1 = [x for x in parent_ids if x in module_ids_1]
            good_2 = [x for x in parent_ids if x in module_ids_2]
            if not good_1 or not good_2:
                continue
            
            # Found one!
            module_id_1 = good_1[0]
            module_id_2 = good_2[0]
            i = parent_ids.index(module_id_1)
            j = parent_ids.index(module_id_2)
            attr_value_1 = attr_values[i]
            attr_value_2 = attr_values[j]
            x = attr_name, module_id_1, attr_value_1, \
                module_id_2, attr_value_2, data_id
            return x
    return None


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


def _intlist2bits(int_list):
    bits = 0
    for i in int_list:
        bits = bits | (1<<i)
    return bits


def _prune_alternate_attributes2(
    network, custom_attributes, paths, nodeid2parents):
    # If a module takes DataNodes with several values, then we only
    # need to calculate one value.
    from Betsy import bie3

    # Fastq.trimmed=no                              -> align
    # Fastq.trimmed=no -> trim -> Fastq.trimmed=yes -> align
    # 
    # 1.  Look in two pipelines for:
    #     - Same ModuleNode.
    #     - Different DataNodes can go into that Node.
    #     - One DataNode is upstream of the other DataNode.
    # 2.  Keep the pipeline that is shorter (or found in
    #     custom_attributes).
    # 3.  Print a message showing which value was chosen.

    path_ids = {}
    transitions = {}
    for path in paths:
        path_ids.update(path.node_ids)
        for node_id, next_ids in path.transitions.iteritems():
            x = transitions.get(node_id, set()).union(next_ids)
            transitions[node_id] = x
    ancestors = bie3._make_ancestor_dict(network)

    # Find module nodes that can take DataNodes with different
    # attributes.
    alternates = []
    for i in range(len(network.nodes)):
        if not isinstance(network.nodes[i], bie3.ModuleNode):
            continue
        x = _list_alternate_attributes2(
            network, i, custom_attributes, path_ids, transitions,
            nodeid2parents)
        alternates.extend(x)

    # For each of the alternates, select the desired transition and
    # rule out the others.
    desired_alternates = [None] * len(alternates)
    for i, x in enumerate(alternates):
        module_id, parent_ids, attr_name, attr_values = x
        # If there is a value specified in the custom attribute, use
        # that one.
        for x in custom_attributes:
            if x.datatype != network.nodes[parent_ids[0]].datatype:
                continue
            if x.name != attr_name:
                continue
            if x.value in attr_values:
                desired_alternates[i] = attr_values.index(x.value)
                break
        if desired_alternates[i] is not None:
            continue
        # Hack: If the options are "no" and "yes", choose "no".
        if sorted(attr_values) == ["no", "yes"]:
            desired_alternates[i] = attr_values.index("no")
        if desired_alternates[i] is not None:
            continue
        # If one is an ancestor of the other, choose the ancestor.
        if len(parent_ids) == 2:
            if parent_ids[0] in ancestors.get(parent_ids[1], []):
                desired_alternates[i] = 0
            elif parent_ids[1] in ancestors.get(parent_ids[0], []):
                desired_alternates[i] = 1
        if desired_alternates[i] is not None:
            continue
        # Choose the one that comes first in the alphabet.
        x = sorted(attr_values)[0]
        desired_alternates[i] = attr_values.index(x)

    # Make a list of the transitions to avoid in the network.
    to_prune = {}  # transitions to avoid
    for i, x in enumerate(alternates):
        module_id, parent_ids, attr_name, attr_values = x
        for j, id_ in enumerate(parent_ids):
            if j == desired_alternates[i]:
                continue
            to_prune[(parent_ids[j], module_id)] = 1

    # Find the paths to prune.
    delete = {}
    for i, path in enumerate(paths):
        p = False
        for parent_id, child_id in to_prune:
            if child_id in path.transitions.get(parent_id, []):
                p = True
                break
        if p:
            delete[i] = 1
    
    paths = [x for (i, x) in enumerate(paths) if i not in delete]
    return paths


def _list_alternate_attributes2(
    network, module_id, custom_attributes, path_ids, transitions,
    nodeid2parents):
    # Return list of (module_id, parent_ids, attr_name, attr_values).
    import itertools
    from Betsy import bie3
    
    # Fastq.trimmed=no                              -> align
    # Fastq.trimmed=no -> trim -> Fastq.trimmed=yes -> align

    # At least 2 DataNodes of the same DataType must transition
    # into this ModuleNode.
    x = nodeid2parents[module_id]
    if len(x) < 2:
        return []
    x = [x for x in x if x in path_ids]
    if len(x) < 2:
        return []
    x = [x for x in x if module_id in transitions.get(x, [])]
    if len(x) < 2:
        return []
    parent_ids = x

    # Keep only if there are at least two nodes with the same datatype.
    datatype_names = [network.nodes[x].datatype.name for x in parent_ids]
    counts = {}
    for n in datatype_names:
        counts[n] = counts.get(n, 0) + 1
    x = [x for (i, x) in enumerate(parent_ids)
         if counts[datatype_names[i]] >= 2]
    parent_ids = x
    if len(parent_ids) < 2:
        return []

    # Previous DataNodes must be in different combinations.
    inputnum2parentids = {}  # which input to module -> list of parent IDs
    combos = bie3._bc_to_input_ids(
        network, module_id, custom_attributes, nodeid2parents=nodeid2parents)
    for x in itertools.product(combos, parent_ids):
        combo, parentid = x
        if parentid not in combo:
            continue
        input_num = combo.index(parentid)
        if input_num not in inputnum2parentids:
            inputnum2parentids[input_num] = []
        inputnum2parentids[input_num].append(parentid)

    # Previous DataNodes must have different attribute values.
    alt_attributes = []
    for input_num, parent_ids in inputnum2parentids.iteritems():
        if len(parent_ids) < 2:
            continue
        id_1 = parent_ids[0]
        node_1 = network.nodes[id_1]
        attr_names = []
        for id_2 in parent_ids[1:]:
            node_2 = network.nodes[id_2]
            x, x, diff_attrs = bie3._score_same_data(node_1, node_2)
            # How to handle multiple different attributes?
            if len(diff_attrs) != 1:
                continue
            attr_names.append(diff_attrs[0])
        # Must have all the same different attributes.
        if len(attr_names) != len(parent_ids)-1:
            continue
        all_same = True
        for i in range(1, len(attr_names)):
            if attr_names[i] != attr_names[0]:
                all_same = False
                break
        if not all_same:
            continue
        name = attr_names[0]

        # Make sure the attribute is not SAME_AS_CONSTRAINT.  Module
        # should not be passing the value of this attribute along.
        # Not sure whether BASED_ON_DATA should be ignored too?
        module_passes_attr = False
        for cons in network.nodes[module_id].consequences:
            if cons.name == name and \
                   cons.behavior in [bie3.SAME_AS_CONSTRAINT]:
                module_passes_attr = True
                break
        if module_passes_attr:
            continue
        
        values = [network.nodes[x].attributes[name] for x in parent_ids]
        x = module_id, parent_ids, name, values
        alt_attributes.append(x)
    return alt_attributes


def _prune_superset_pipelines(network, paths):
    # Remove pipelines that are just supersets of another pipeline.
    import itertools

    path2length = [len(x.node_ids) for x in paths]

    # Convert node_ids to bitwise numbers for faster superset
    # comparison.
    paths_node_ids_bit = [_intlist2bits(x.node_ids) for x in paths]

    superset = []
    for (i, j) in itertools.product(range(len(paths)), range(len(paths))):
        if i == j:
            continue
        # Optimization: Check for length here to avoid function call.
        if path2length[i] <= path2length[j]:
            continue
        # Optimization: Check for superset here.
        #if not set(node_ids_1).issuperset(node_ids_2):
        b1, b2 = paths_node_ids_bit[i], paths_node_ids_bit[j]
        if (b1 | b2) != b1:
            continue
        if _is_superset_pipeline(network, paths[i], paths[j]):
            superset.append(i)
    superset = {}.fromkeys(superset)
    paths = [x for (i, x) in enumerate(paths) if i not in superset]
    return paths


def _is_superset_pipeline(network, path_1, path_2):
    # Test if path_1 is superset, given path_2.  I.e. path_1 has more
    # processing steps that are not necessary in path_2.
    from Betsy import bie3
    
    # If they're superset, then the start nodes must be different.
    # Actually, this is not true.  They can have the same start nodes,
    # but different internal paths.
    #if sorted(start_ids_1) == sorted(start_ids_2):
    #    return False

    # Comparisons now done in calling function.
    ## If they're superset, path_1 must be a superset of path_2.
    #if len(node_ids_1) <= len(node_ids_2):
    #    # Returns from here ~60% of the time.
    #    return False
    ## Takes about 50% of the time in this function.
    #if not set(node_ids_1).issuperset(node_ids_2):
    #    return False

    # At least one of the extra nodes must be a start node.
    # Actually, this is not true.  They can have the same start nodes,
    # but different internal paths.
    #x = list(set(node_ids_1).difference(node_ids_2))
    #x = [x for x in x if x in start_ids_1]
    #super_start_ids_1 = x
    #if not super_start_ids_1:
    #    return False

    # All paths starting from start_ids_1 must go through one of
    # start_ids_2.
    start_ids_1 = [x for x in path_1.start_ids if x is not None]
    start_ids_2 = [x for x in path_2.start_ids if x is not None]
    for start_id in start_ids_1:
        if not _does_path_go_through(
            network, start_id, path_1.node_ids, start_ids_2):
            return False

    # Looks like a superset, but it's not.
    # is_compressed -> Fastq (yes) -> uncompress -> Fastq (no)
    # is_compressed                              -> Fastq (no)
    super_ids = [x for x in path_1.node_ids if x not in path_2.node_ids]
    assert super_ids
    # If any of the nodes in the super_ids list is downstream of a
    # Module with a BASED_ON_DATA Consequence, then this is not a
    # superset.
    for (node_id, node) in enumerate(network.nodes):
        if node_id in super_ids: # ignore if this is part of the superset
            continue
        if not isinstance(node, bie3.ModuleNode):
            continue
        x = [x for x in node.consequences if x.behavior == bie3.BASED_ON_DATA]
        if not x:
            continue
        next_ids = network.transitions.get(node_id, [])
        x = [x for x in next_ids if x in super_ids]
        if x:
            return False
    return True


def _does_path_go_through(network, start_id, good_ids, intermediate_ids):
    # Do all paths starting from start_id go through intermediate_ids?
    stack = [start_id]
    while stack:
        node_id = stack.pop(0)
        if node_id in intermediate_ids:
            continue
        x = network.transitions.get(node_id, [])
        next_ids = [x for x in x if x in good_ids]
        if not next_ids:
            # Goes to end of network.  Did not encounter intermediate_ids.
            return False
        stack.extend(next_ids)
    return True


def _prune_parallel_pipelines(network, paths, nodeid2parents):
    # Remove pipelines that are parallel to another pipeline.
    global SUBPATH_CACHE   # for optimization

    SUBPATH_CACHE = {}

    #plot_pipelines(
    #    "pipeline", network, paths, {}, max_pipelines=16, verbose=True)

    path_lengths = [len(x.node_ids) for x in paths]

    # Prune parallel pipelines with the same lengths first to get a
    # smaller set of pipelines.  Then, do the more complex pruning.
    prune = {}
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            if i in prune and j in prune:
                continue
            
            # Optimization: Don't check if the lengths are not the
            # same.
            if path_lengths[i] != path_lengths[j]:
                continue
            p = _is_parallel_pipeline1(
                network, paths[i], paths[j], nodeid2parents)
            if not p:
                p = _is_parallel_pipeline2(
                    network, paths[i], paths[j], nodeid2parents)
            #if not p:
            #    p = _is_parallel_pipeline3(
            #        network, paths, i, j, nodeid2parents)
            if not p:
                continue
            if p == 1:
                prune[i] = 1
            else:
                prune[j] = 1
    #orig_path_ids = [
    #    x for (i, x) in enumerate(range(len(paths))) if i not in prune]
    paths = [x for (i, x) in enumerate(paths) if i not in prune]
    path_lengths = [x for (i, x) in enumerate(path_lengths) if i not in prune]

    # Now do the more computationally expensive pruning.
    prune = {}
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            if i in prune and j in prune:
                continue
            p = _is_parallel_pipeline3(
                network, paths, i, j, nodeid2parents)
            if not p:
                continue
            if p == 1:
                prune[i] = 1
            else:
                prune[j] = 1
                
    paths = [x for (i, x) in enumerate(paths) if i not in prune]
    return paths


def _is_parallel_pipeline1(network, path_1, path_2, nodeid2parents):
    # Test if path_1 is parallel to path_2.  If not parallel, returns
    # False.  Otherwise, returns a 1 or 2 indicating which one should
    # be pruned.
    # 
    # Parallel:
    # BAM (no) -> sort by name  -> BAM (name)  -> count_htseq
    # BAM (no) -> sort by coord -> BAM (coord) -> count_htseq
    # Different from:
    # FASTQ (unknown) -> is_compress                              -> FASTQ (no)
    # FASTQ (unknown) -> is_compress -> FASTQ (yes) -> uncompress -> FASTQ (no)
    from Betsy import bie3

    #node_ids_1, start_ids_1, data_indexes_1 = path_1
    #node_ids_2, start_ids_2, data_indexes_2 = path_2

    # If they're parallel, then they have the same number of nodes.
    if len(path_1.node_ids) != len(path_2.node_ids):
        return False
    # They differ by only two nodes.
    unique_ids_1 = [x for x in path_1.node_ids if x not in path_2.node_ids]
    if len(unique_ids_1) != 2:
        return False
    unique_ids_2 = [x for x in path_2.node_ids if x not in path_1.node_ids]
    if len(unique_ids_2) != 2:
        return False
    # The module node is upstream of the data node.
    module_id_1, data_id_1 = unique_ids_1
    module_id_2, data_id_2 = unique_ids_2
    if isinstance(network.nodes[module_id_1], bie3.DataNode):
        module_id_1, data_id_1 = data_id_1, module_id_1
    if isinstance(network.nodes[module_id_2], bie3.DataNode):
        module_id_2, data_id_2 = data_id_2, module_id_2
    if not isinstance(network.nodes[module_id_1], bie3.ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_1], bie3.DataNode):
        return False
    if not isinstance(network.nodes[module_id_2], bie3.ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_2], bie3.DataNode):
        return False
    if data_id_1 not in network.transitions.get(module_id_1, []):
        return False
    if data_id_2 not in network.transitions.get(module_id_2, []):
        return False
    # The module nodes have the same parent.
    parent_ids_1 = nodeid2parents.get(module_id_1, [])
    parent_ids_2 = nodeid2parents.get(module_id_2, [])
    common_ids = [x for x in parent_ids_1 if x in parent_ids_2]
    if not common_ids:
        return False
    # The data nodes have the same child.
    child_ids_1 = network.transitions.get(data_id_1, [])
    child_ids_2 = network.transitions.get(data_id_2, [])
    common_ids = [x for x in child_ids_1 if x in child_ids_2]
    if not common_ids:
        return False
    # The data nodes have the same type.
    data_1 = network.nodes[data_id_1]
    data_2 = network.nodes[data_id_2]
    if data_1.datatype != data_2.datatype:
        return False
    # The data nodes differ by only one attribute.
    assert sorted(data_1.attributes) == sorted(data_2.attributes)
    diff = [x for x in data_1.attributes
            if data_1.attributes[x] != data_2.attributes[x]]
    if len(diff) != 1:
        return False
    # Prune the one whose value comes later in the alphabet.
    name = diff[0]
    value_1 = data_1.attributes[name]
    value_2 = data_2.attributes[name]
    if value_1 < value_2:
        return 2
    return 1


def _is_parallel_pipeline2(network, path_1, path_2, nodeid2parents):
    # Test if path_1 is parallel to path_2.  If not parallel, returns
    # False.  Otherwise, returns a 1 or 2 indicating which one should
    # be pruned.
    # 
    # Parallel:
    # BAM (no) -> sort -> BAM (sort=y) -> addgroup ->  BAM (sort=y, group=y)
    #          -> addgroup -> BAM (group=y) -> sort ->
    from Betsy import bie3

    # If they're parallel, then they have the same number of nodes.
    if len(path_1.node_ids) != len(path_2.node_ids):
        return False
    # They differ by only three nodes.
    unique_ids_1 = [x for x in path_1.node_ids if x not in path_2.node_ids]
    if len(unique_ids_1) != 3:
        return False
    unique_ids_2 = [x for x in path_2.node_ids if x not in path_1.node_ids]
    if len(unique_ids_2) != 3:
        return False
    
    # The module node is upstream of the data node.
    top_id_1, data_id_1, bottom_id_1 = unique_ids_1
    top_id_2, data_id_2, bottom_id_2 = unique_ids_2
    if isinstance(network.nodes[top_id_1], bie3.DataNode):
        top_id_1, data_id_1 = data_id_1, top_id_1
    if isinstance(network.nodes[bottom_id_1], bie3.DataNode):
        bottom_id_1, data_id_1 = data_id_1, bottom_id_1
    if data_id_1 in network.transitions[bottom_id_1]:
        top_id_1, bottom_id_1 = bottom_id_1, top_id_1
    if isinstance(network.nodes[top_id_2], bie3.DataNode):
        top_id_2, data_id_2 = data_id_2, top_id_2
    if isinstance(network.nodes[bottom_id_2], bie3.DataNode):
        bottom_id_2, data_id_2 = data_id_2, bottom_id_2
    if data_id_2 in network.transitions[bottom_id_2]:
        top_id_2, bottom_id_2 = bottom_id_2, top_id_2

    # Make sure the node types are correct.
    if not isinstance(network.nodes[top_id_1], bie3.ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_1], bie3.DataNode):
        return False
    if not isinstance(network.nodes[bottom_id_1], bie3.ModuleNode):
        return False
    if not isinstance(network.nodes[top_id_2], bie3.ModuleNode):
        return False
    if not isinstance(network.nodes[data_id_2], bie3.DataNode):
        return False
    if not isinstance(network.nodes[bottom_id_2], bie3.ModuleNode):
        return False

    # Make sure the connections are correct.
    if data_id_1 not in path_1.transitions.get(top_id_1, []):
        return False
    if bottom_id_1 not in path_1.transitions.get(data_id_1, []):
        return False
    if data_id_2 not in path_2.transitions.get(top_id_2, []):
        return False
    if bottom_id_2 not in path_2.transitions.get(data_id_2, []):
        return False

    # The top nodes have the same parent.
    parent_ids_1 = nodeid2parents.get(top_id_1, [])
    parent_ids_2 = nodeid2parents.get(top_id_2, [])
    x = [x for x in parent_ids_1 if x in parent_ids_2]
    x = [x for x in x if x in path_1.node_ids and x in path_2.node_ids]
    common_ids = x
    if not common_ids:
        return False
    
    # The bottom nodes have the same child.
    child_ids_1 = network.transitions.get(bottom_id_1, [])
    child_ids_2 = network.transitions.get(bottom_id_2, [])
    x = [x for x in child_ids_1 if x in child_ids_2]
    x = [x for x in x if x in path_1.node_ids and x in path_2.node_ids]
    common_ids = x
    if not common_ids:
        return False
    
    # The data nodes have the same type.
    data_1 = network.nodes[data_id_1]
    data_2 = network.nodes[data_id_2]
    if data_1.datatype != data_2.datatype:
        return False
    
    # The data nodes differ by two attributes.
    assert sorted(data_1.attributes) == sorted(data_2.attributes)
    diff = [x for x in data_1.attributes
            if data_1.attributes[x] != data_2.attributes[x]]
    if len(diff) != 2:
        return False
    
    # Prune the one whose value comes later in the alphabet.
    name = sorted(diff)[0]
    value_1 = data_1.attributes[name]
    value_2 = data_2.attributes[name]
    if value_1 < value_2:
        return 2
    return 1


def _is_parallel_pipeline3(
    network, paths, path_id_1, path_id_2, nodeid2parents):
    # Test if path_1 is parallel to path_2.  If not parallel, returns
    # False.  Otherwise, returns a 1 or 2 indicating which one should
    # be pruned.
    # More careful and slower version of _is_parallel_pipeline3.
    # 
    # Parallel:
    # DataNode -> sort_coord -> mark_dup -> add_read_group -> DataNode
    #          -> add_read_group -> sort_coord -> mark_dup ->
    from Betsy import bie3
    
    path_1, path_2 = paths[path_id_1], paths[path_id_2]
    x = _compare_paths(network, path_1, path_2)
    shared_ids, unique_ids_1, unique_ids_2 = x
    # shared_ids is almost always a very long list.
    
    # For downstream tests, calculate some useful variables describing
    # the networks.
    x1 = _build_subpath(
        network, paths, path_id_1, unique_ids_1, nodeid2parents)
    x2 = _build_subpath(
        network, paths, path_id_2, unique_ids_2, nodeid2parents)
    subpath_1, subpath_2 = x1, x2

    # There should only be 1 top and 1 bottom node.  This test removes
    # 95% of the possibilities.  The others don't filter as much.
    if len(subpath_1.top_node_ids) != 1 or len(subpath_1.bottom_node_ids) != 1:
        return False
    if len(subpath_2.top_node_ids) != 1 or len(subpath_2.bottom_node_ids) != 1:
        return False

    # The top should be modules.  The bottom might be a DataNode if
    # the last step is the same.
    # path1 -> BamFolder.sorted=name  -> sort_by_contig
    # path2 -> BamFolder.sorted=coord -> sort_by_contig
    x = [subpath_1.top_node_ids[0], subpath_2.top_node_ids[0]]
    for x in x:
        if not isinstance(network.nodes[x], bie3.ModuleNode):
            return False

    # The parents of top and children of bottom should be the same.
    if sorted(subpath_1.parents_of_top) != sorted(subpath_2.parents_of_top):
        return False
    if sorted(subpath_1.children_of_bottom) != \
           sorted(subpath_2.children_of_bottom):
        return False

    # Any attributes based on data should be the same.
    bod1 = subpath_1.based_on_data
    bod2 = subpath_2.based_on_data
    if sorted(bod1) != sorted(bod2):
        return False
    for name in bod1:
        if bod1[name] != bod2[name]:
            return False

    # Be sure to always take the same parallel pipeline.  Use the one
    # whose module name comes first in the alphabet.
    id_1, id_2 = subpath_1.top_node_ids[0], subpath_2.top_node_ids[0]
    node_1, node_2 = network.nodes[id_1], network.nodes[id_2]
    value_1, value_2 = node_1.name, node_2.name
    # Prune the one whose value comes later in the alphabet.
    if value_1 < value_2:
        return 2
    return 1


def _compare_paths(network, path_1, path_2):
    from Betsy import bie3
    
    # If a node occurs in only one pathway, then it is unique for
    # sure.
    unique_ids_1 = _dict_diff(path_1.node_ids, path_2.node_ids, _as_dict=True)
    unique_ids_2 = _dict_diff(path_2.node_ids, path_1.node_ids, _as_dict=True)
    shared_ids = _dict_diff(path_1.node_ids, unique_ids_1, _as_dict=True)

    # If a Module occurs in both pathways, it might be unique if it
    # has different transitions.
    # sort_by_contig -> BamFolder.readgroups (no)
    # sort_by_contig -> BamFolder.readgroups (yes)
    # i.e same module, but different uses
    
    for node_id in list(shared_ids):
        if not isinstance(network.nodes[node_id], bie3.ModuleNode):
            continue
        #assert node_id in path_1.transitions
        #assert node_id in path_2.transitions
        children_1 = path_1.transitions[node_id]
        children_2 = path_2.transitions[node_id]
        # Actually takes more time to compare lengths first.  Sorted
        # pretty fast.
        #if len(children_1) != len(children_2):
        #    continue
        #if sorted(children_1) == sorted(children_2):
        if children_1 == children_2:
            continue
        del shared_ids[node_id]
        unique_ids_1[node_id] = 1
        unique_ids_2[node_id] = 1
    return shared_ids, unique_ids_1, unique_ids_2


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


class PathSubset:
    # Have a network that includes all nodes.
    # A subset of those nodes comprises a path through the network.
    # A subset of the path comprises a subpath
    #
    # path_node_ids
    # sub_node_ids
    # path_transitions
    # 
    # top_node_ids         list of sub node IDs (no parents in subpath)
    # bottom_node_ids      list of sub node IDs (no children in subpath)
    # parents_of_top       list of path node IDs (above top_node_ids)
    # children_of_bottom   list of path node IDs (below top_node_ids)
    #
    # sub_parents      dict of sub_node_id -> ids of parents in subpath
    # sub_children     dict of sub_node_id -> ids of children in subpath
    # path_parents     dict of sub_node_id -> ids of parents not in subpath
    # path_children    dict of sub_node_id -> ids of children not in subpath
    #
    # based_on_data    dict of attrname -> value; from subpath
    # 
    # Since this is very expensive computationally, calculate these
    # variables on demand.
    def __init__(
        self, network, path_node_ids, path_transitions, sub_node_ids,
        nodeid2parents):
        # Make sure these are dicts, or will be really slow.
        assert type(path_node_ids) is type({})
        assert type(sub_node_ids) is type({})

        self._network = network
        self._nodeid2parents = nodeid2parents
        self.path_node_ids = path_node_ids
        self.path_transitions = path_transitions
        self.sub_node_ids = sub_node_ids
    def __getattr__(self, name):
        var2info = {
            "sub_parents" : (
                "_sub_parents", self._calc_sub_parents_and_top_nodes),
            "sub_children" : (
                "_sub_children", self._calc_sub_children_and_bottom_nodes),
            "path_parents" : ("_path_parents", self._calc_path_parents),
            "path_children" : ("_path_children", self._calc_path_children),
            "top_node_ids" : (
                "_top_node_ids", self._calc_sub_parents_and_top_nodes),
            "bottom_node_ids" : (
                "_bottom_node_ids", self._calc_sub_children_and_bottom_nodes),
            "parents_of_top" : ("_parents_of_top", self._calc_parents_of_top),
            "children_of_bottom" : (
                "_children_of_bottom", self._calc_children_of_bottom),
            "based_on_data" : ("_based_on_data", self._calc_based_on_data),
            }
        if name not in var2info:
            raise AttributeError, name
        x = var2info[name]
        member_name, fn = x
        if not hasattr(self, member_name):
            fn()
            assert hasattr(self, member_name)
            #setattr(self, member_name, fn())
        return getattr(self, member_name)
    def _calc_sub_parents_and_top_nodes(self):
        pars = {}
        top_node_ids = []
        for node_id in self.sub_node_ids:
            x = self._nodeid2parents.get(node_id, [])
            x = [x for x in x if x in self.sub_node_ids]
            if x:
                pars[node_id] = x
            else:
                top_node_ids.append(node_id)
        self._sub_parents = pars
        self._top_node_ids = top_node_ids
    def _calc_sub_children_and_bottom_nodes(self):
        children = {}
        bottom_node_ids = []
        for node_id in self.sub_node_ids:
            x = self.path_transitions.get(node_id, [])
            x = [x for x in x if x in self.sub_node_ids]
            if x:
                children[node_id] = x
            else:
                bottom_node_ids.append(node_id)
        self._sub_children = children
        self._bottom_node_ids = bottom_node_ids
    def _calc_path_parents(self):
        pars = {}
        for node_id in self.sub_node_ids:
            x = self._nodeid2parents.get(node_id, [])
            x = [x for x in x if x not in self.sub_node_ids]
            if x:
                pars[node_id] = x
        self._path_parents = pars
    def _calc_path_children(self):
        pars = {}
        for node_id in self.sub_node_ids:
            x = self.path_transitions.get(node_id, [])
            x = [x for x in x if x not in self.sub_node_ids]
            if x:
                pars[node_id] = x
        self._path_children = pars
    #def _calc_top_node_ids(self):
    #    sub_parents = self.sub_parents
    #    x = [x for x in self.sub_node_ids if not sub_parents.get(x)]
    #    self._top_node_ids = x
    #def _calc_bottom_node_ids(self):
    #    sub_children = self.sub_children
    #    x = [x for x in self.sub_node_ids if not sub_children.get(x)]
    #    self._bottom_node_ids = x
    def _calc_parents_of_top(self):
        parents_of_top = {}
        for node_id in self.top_node_ids:
            x = [x for x in self._nodeid2parents.get(node_id, [])]
            # Is this right?
            x = [x for x in x if x in self.path_node_ids]
            for x in x:
                parents_of_top[x] = 1
        self._parents_of_top = parents_of_top.keys()
    def _calc_children_of_bottom(self):
        children_of_bottom = {}
        for node_id in self.bottom_node_ids:
            x = [x for x in self.path_transitions.get(node_id, [])]
            x = [x for x in x if x in self.path_node_ids]
            for x in x:
                children_of_bottom[x] = 1
        self._children_of_bottom = children_of_bottom.keys()
    def _calc_based_on_data(self):
        # For all the nodes in the subpath, get the values of the
        # attributes that are set based on data.
        from Betsy import bie3

        based_on_data = {}
        for node_id in self.sub_node_ids:
            node = self._network.nodes[node_id]
            if not isinstance(node, bie3.ModuleNode):
                continue
            # Make a list of the attributes that are BASED_ON_DATA.
            x = [
                x for x in node.consequences if
                x.behavior == bie3.BASED_ON_DATA]
            attr_names = [x.name for x in x]
            if not attr_names:
                continue

            # Find the values of the attributes in the children.
            x = self.path_transitions.get(node_id, [])
            # Is this necessary?
            child_ids = [x for x in x if x in self.path_node_ids]
            for name in attr_names:
                for node_id in child_ids:
                    node = self._network.nodes[node_id]
                    if name not in node.attributes:
                        continue
                    value = node.attributes[name]
                    assert name not in based_on_data
                    based_on_data[name] = value
                    #if name not in based_on_data:
                    #    based_on_data[name] = []
                    #based_on_data[name].append(value)
        self._based_on_data = based_on_data



SUBPATH_CACHE = {}  # (path_id, sub_node_ids) -> PathSubset
def _build_subpath(network, paths, path_id, sub_node_ids, nodeid2parents):
    global SUBPATH_CACHE
    key = path_id, tuple(sorted(sub_node_ids))
    if key not in SUBPATH_CACHE:
        x = _build_subpath_h(
            network, paths, path_id, sub_node_ids, nodeid2parents)
        SUBPATH_CACHE[key] = x
    return SUBPATH_CACHE[key]

    
def _build_subpath_h(network, paths, path_id, sub_node_ids, nodeid2parents):
    path = paths[path_id]
    return PathSubset(
        network, path.node_ids, path.transitions, sub_node_ids, nodeid2parents)


def check_input_files(network, in_data_nodes, user_options, paths,
                      network_png, verbose):
    import os
    import sys

    print "Looking for input files."
    sys.stdout.flush()
    missing = False
    for x in in_data_nodes:
        if not x.identifier:
            print "no file given for: %s." % (x.data.datatype.name)
            missing = True
        elif not os.path.exists(x.identifier):
            print "File not found: %s." % x.identifier
            missing = True

    if not missing:
        #print "all input files found."
        #sys.stdout.flush()
        return True

    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        verbose=verbose)
    return []


def manually_verify_network(
    network, user_options, paths, run, network_png, verbose):
    import sys
    if run:
        return True

    print "Please review the network to make sure the analysis is correct."
    print "Add --run when ready to run the analysis."
    sys.stdout.flush()

    plot_network_show_pipelines(
        network_png, network, paths, user_options=user_options,
        verbose=verbose)
    return False
    

def plot_network_show_pipelines(filename, network, paths, **keywds):
    from Betsy import bie3

    x1 = [x.node_ids for x in paths]
    # Hack: For optimization, sometimes node_ids is dict.
    if x1 and type(x1[0]) is type({}):
        x1 = [x.keys() for x in x1]
    x2 = [x.start_ids for x in paths]
    x1 = bie3._uniq(bie3._flatten(x1))
    x2 = bie3._uniq(bie3._flatten(x2))
    x2 = [x for x in x2 if x is not None]
    all_node_ids = x1
    all_start_ids = x2

    transitions = {}
    for path in paths:
        for node_id, next_ids in path.transitions.iteritems():
            for next_id in next_ids:
                transitions[(node_id, next_id)] = 1
    plot_network(
        filename, network, bold=all_node_ids,
        bold_transitions=transitions, highlight_green=all_start_ids,
        highlight_purple=all_node_ids, **keywds)


def plot_network(
    filename, network, user_options=None, bold=[], bold_transitions=[], 
    highlight_green=[], highlight_orange=[], highlight_purple=[],
    highlight_yellow=[], verbose=False):
    import sys
    from Betsy import bie3

    if filename is None:
        return
    print "Plotting network to %s." % filename
    sys.stdout.flush()
    bie3.plot_network_gv(
        filename, network, options=user_options, bold=bold,
        bold_transitions=bold_transitions, 
        highlight_green=highlight_green, highlight_orange=highlight_orange,
        highlight_purple=highlight_purple, highlight_yellow=highlight_yellow,
        verbose=verbose)


def plot_pipelines(filestem, network, paths, user_options, max_pipelines=None,
                   verbose=False):
    if max_pipelines is not None:
        paths = paths[:max_pipelines]
    for i, path in enumerate(paths):
        filename = "%s-%02d.png" % (filestem, i)
        plot_network_show_pipelines(
            filename, network, [path], user_options=user_options,
            verbose=verbose)


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


def _parse_dattr(dattr_str):
    # Format: <datatype>.<key>=<value>
    # Return <datatype>, <key>, <value>.

    err_msg = "--dattr should be <datatype>.<key>=<value>"
    x = dattr_str.split(".", 1)
    assert len(x) == 2, err_msg
    datatype, x = x
    x = x.split("=", 1)
    assert len(x) == 2, err_msg
    key, value = x
    return datatype, key, value


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
            assert inputs
            assert len(args) >= i+1
            dattr = args[i+1]
            i += 2
            datatype, key, value = _parse_dattr(dattr)
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
            datatype, key, value = _parse_dattr(dattr)
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
    import getpass

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
        "FLAG: --input\n"
        "5.  Work out the attributes of the inputs.\n"
        "    Shows the detailed information about each node.\n"
        "6.  System makes sure each node can be found in the network.\n"
        "    If not, try to diagnose differences.\n"
        "7.  System makes sure this is a complete set of input nodes.\n"
        "8.  System makes sure all required module attributes are given.\n"
        "9.  System makes sure all input files are provided.\n"
        "10.  Actually run the analysis.\n"
        "FLAG:--run\n"
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
        '--num_cores', default=4, type=int,
        help='number of cores used in the processing')
    DEFAULT_MAX_INPUTS = 5
    parser.add_argument(
        '--max_inputs', default=DEFAULT_MAX_INPUTS, type=int,
        help="Maximum number of inputs to be shown (default %d).  "
        "(For optimization)" % DEFAULT_MAX_INPUTS)
    parser.add_argument(
        "--save_failed_data", action="store_true",
        help="If a module failed, do not clean up its working files.")
    parser.add_argument(
        "--cache_input_files", action="store_true",
        help="Save a copy of the input files.  "
        "Helpful for very big files.")

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
        "--dattr that refers to nodes not in the network should be "
        "given after the --output.  Format: <datatype>.<key>=<value>.")
    group.add_argument(
        '--output_file', help='file or folder of output result')

    group = parser.add_argument_group(title="Outfiles")
    group.add_argument(
        '--network_png', help='generate the output network png file')
    group.add_argument(
        '--sparse_network_png', action="store_true",
        help="Leave out details in network plot.")
    #group.add_argument(
    #    '--network_text', help='generate the output network text file')
    group.add_argument(
        '--network_json', help='generate the output network json file')
    #parser.add_argument(
    #    '--clobber', action='store_const', const=True, default=False,
    #    help='overwrite the output_data if it already exists')

    print "Starting rule engine."
    sys.stdout.flush()

    # Parse the arguments.
    args = parser.parse_args()
    input_list, x = _parse_args(sys.argv)
    outtype, out_identifier, out_attributes = x
    verbose = (not args.sparse_network_png)    

    print "Checking parameters."
    sys.stdout.flush()

    args.clobber = True
    assert args.num_cores > 0, "num_cores should be greater than 0"
    assert args.max_inputs > 0 and args.max_inputs < 20
    if args.inputs_complete or args.attrs_complete or args.run:
        if args.run:
            x = "--run"
        elif args.attrs_complete:
            x = "--attrs_complete"
        else:
            x = "--inputs_complete"
        assert input_list, "%s given, but no --input." % x
    ## TODO: Make sure args.exclude_input is valid.
    ## args.exclude_input

    # Make sure configuration directory exists.
    if not os.path.exists(config.OUTPUTPATH):
        print "Making BETSY working path: %s." % config.OUTPUTPATH
        os.mkdir(config.OUTPUTPATH)

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
        x = module_utils.IdentifiedDataNode(in_data, identifier)
        in_data_nodes.append(x)

    # test outtype and build the list of custom_attributes.
    custom_attributes = []  # List of bie3.Attribute objects.
    if outtype:
        # Get the custom_attributes.  These are the attributes from the
        # output data (out_attributes), plus the attributes from the
        # input data (to guide the inferencing engine).
        attrs = []  # list of (datatype, key, value)
        for x in input_list:
            intype, identifier, attributes = x
            for k, v in attributes:
                attrs.append((intype, k, v))
        attrs.extend(out_attributes)
        #for x in out_attributes:
        for x in attrs:
            datatype, key, value = x
            assert hasattr(rulebase, datatype), \
                   "Unknown datatype: %s" % datatype
            fn = getattr(rulebase, datatype)
            custom_attributes.append(bie3.Attribute(fn, key, value))

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
    network = generate_network(rulebase, outtype, custom_attributes)
    # Step 3: Make sure some inputs are provided.
    if not check_inputs_provided(
        network, in_data_nodes, custom_attributes, user_options,
        args.max_inputs, args.network_png, verbose):
        return
    # Step 4: Make sure each of the input nodes match a node in the
    # network.
    x = check_inputs_in_network(
        network, user_options, in_data_nodes, args.network_png, verbose)
    inputs_ok, data_node_ids = x
    if not inputs_ok:
        return
    # Step 5: Search for pipelines that can be run given the INPUT
    # nodes.
    paths = build_pipelines(
        network, user_options, in_data_nodes, data_node_ids, custom_attributes,
        args.max_inputs, args.network_png, verbose)
    #plot_pipelines(
    #    "pipeline", network, paths[:1], user_options, max_pipelines=16,
    #    verbose=True)
    if not paths:
        return
    # Step 6: Make sure required attributes are given.
    paths = check_attributes_complete(
        network, user_options, paths, args.network_png, verbose)
    if not paths:
        return
    # DEBUG: Print out each of the pipelines.
    #plot_pipelines(
    #    "pipeline", network, paths[:1], user_options, max_pipelines=16,
    #    verbose=True)
    # Step 7: Prune undesired pipelines.
    paths = prune_pipelines(
        network, user_options, custom_attributes, paths, args.network_png,
        verbose)
    if not paths:
        return
    # DEBUG: Print out each of the pipelines.
    #plot_pipelines(
    #    "pipeline", network, paths, user_options, max_pipelines=32,
    #    verbose=True)
        
    # Step 8: Look for input files.
    if not check_input_files(
        network, in_data_nodes, user_options, paths, args.network_png,
        verbose):
        return

    #print "The network (%d nodes) is complete." % len(network.nodes)
    #print "There are %d possible pipelines." % len(paths)
    x = "core"
    if args.num_cores > 1:
        x = "cores"
    print "Ready to go!  Will run the analysis using a maximum of %d %s." % (
        args.num_cores, x)
    # Step 9: Manual verification of the network.
    if not manually_verify_network(
        network, user_options, paths, args.run, args.network_png, verbose):
        return

    print "Running the analysis."
    sys.stdout.flush()
    clean_up = not args.save_failed_data
    node_dict = output_file = None
    try:
        x = rule_engine_bie3.run_pipeline(
            network, in_data_nodes, custom_attributes, user_options, paths,
            user=args.user, job_name=args.job_name, clean_up=clean_up,
            num_cores=args.num_cores)
        if x:
            node_dict, transitions, output_file = x
    except AssertionError, x:
        if str(x).startswith("Inference error"):
            node_ids = rule_engine_bie3.DEBUG_POOL.keys()
            plot_network(
                args.network_png, network, user_options=user_options,
                highlight_green=node_ids, verbose=verbose)
        raise

    # Draw out the network.
    node_ids = []
    if node_dict:
        node_ids = node_dict.keys()
    plot_network(
        args.network_png, network, user_options=user_options,
        bold=node_ids, bold_transitions=transitions,
        highlight_green=node_ids, verbose=verbose)
    if args.network_json:
        print "Writing network in json format: %s." % args.network_json
        bie3.write_network(args.network_json, network)
    #if args.network_text:
    #    print "Writing detailed network: %s." % args.network_text
    #    bie3.print_network(network, outhandle=args.network_text)

    if output_file and args.output_file:
        print "Saving results at %s." % args.output_file
        sys.stdout.flush()
        reportlib.copy_file_or_path(output_file, args.output_file)
    print "Done."


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
