#! /usr/bin/env python

# Functions:
# check_output_provided      Step 1: Make sure there's an output.
#   _print_rulebase          To help select output DataType.
#   _pretty_print_datatype
#   _pretty_print_module
# generate_network           Step 2.
# check_inputs_provided      Step 3.
#   _print_input_datatypes   To help select input DataType(s).
# check_datatypes_complete   Step 4.
# check_inputs_complete      Step 5.
#   _assign_data_nodes_to_start_ids
#   _print_global_incompatibilities
#   _print_node_score_table
# check_attributes_complete  Step 6.
# prune_pipelines            Step 7.
#   _prune_redundant_pipelines
#   _is_redumndant_pipeline
#   _does_path_go_through
# check_input_files          Step 8.
# manually_verify_network    Step 9.
# 
# plot_network
# 
# is_complete_input
# get_all_option_names        Get the names of all OptionDef in the rulebase.
# get_required_option_names   Get the required options in a list of modules.
#
# _list_differences_in_nodelists
# _merge_nodelists
# _uniq
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


def generate_network(rulebase, outtype, user_attributes, exclude_input):
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

    # There's a bug in here somehow where impossible networks can be
    # created.  e.g. FastqFolder:orientation="unknown" -> merge_reads
    # -> FastqFolder:orientation="single", even though constraint
    # forces them to be the same.  Happens during complete_network or
    # optimize_network step.  Don't see anymore, because got rid of
    # orientation attribute in FastqFolder.
    network = bie3.make_network(
        rulebase.all_modules, out_datatype, user_attributes)
    network = bie3.complete_network(network, user_attributes)
    network = bie3.optimize_network(network, user_attributes)
    assert network, "Unexpected error: No network generated."
    
    # If the user specified some data types that should be excluded,
    # remove them from the network.
    if exclude_input:
        if len(exclude_input) <= 3:
            x = ", ".join(exclude_input)
            print "Remove as inputs: %s." % x
        else:
            print "Remove %d data types as inputs." % len(exclude_input)
        network = _remove_unwanted_input_nodes(
            network, exclude_input, user_attributes)
    print "Made a network with %d nodes." % len(network.nodes)
    return network


def check_inputs_provided(
    network, exclude_input, user_attributes, user_options, in_data_nodes, 
    max_inputs, network_png, verbose):
    if in_data_nodes:
        return True
    # If there's an output and no inputs, then show the possible
    # data types.  The network is probably too big to show
    # everything.
    print
    print "Selecting --input DataTypes.  Add --inputs_complete when done."
    # TODO: Should not show the same datatype as the output.
    # e.g. if output is a FastQCFolder, don't show an input that
    # just consists of a FastQCFolder.  That's just lame.
    x = [x.data.datatype.name for x in in_data_nodes]
    _print_input_datatypes(
        network, x, exclude_input, user_attributes, max_inputs=max_inputs)
    plot_network(
        network_png, network, user_options=user_options, verbose=verbose)
    return False


def _print_input_datatypes(
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


def check_datatypes_complete(
    network, user_attributes, user_options, in_data_nodes, network_png,
    verbose):
    # Return list of (path, start_ids).  Empty list means no paths
    # found.
    import sys
    from Betsy import bie3
    
    print "Searching the network for pipelines that use --input data types."
    sys.stdout.flush()
    # Make a list of the names of data types for start nodes.
    assert in_data_nodes
    start_dtype_names = [x.data.datatype.name for x in in_data_nodes]

    complete = []    # list of (path, start, missing)
    incomplete = []  # list of (path, start, missing)
    for i, x in enumerate(bie3.find_paths_by_datatypes(
            network, user_attributes, start_dtype_names)):
        path, start_ids, missing_ids = x

        # Ignore paths that don't contain any of the requested data types.
        if not start_ids:
            continue
        # Ignore paths with too many missing data types.
        if len(missing_ids) > 5:
            continue
        
        x = path, start_ids, missing_ids
        if not missing_ids:
            assert len(start_ids) <= len(in_data_nodes)
            complete.append(x)
        else:
            incomplete.append(x)

        ## # DEBUG: print out interesting paths.
        ## if 0 and len(missing_ids) == 0:
        ##     filename = "path-%d.png" % i
        ##     plot_network(
        ##         filename, network, user_options=user_options,
        ##         bold=path, highlight_green=start_ids,
        ##         highlight_purple=missing_ids, verbose=True)

    if complete:
        print "Found %d possible pipelines." % len(complete)
        x = [tuple(x[:2]) for x in complete]
        return x

    # Print a diagnostic message.
    print "No paths found."
    if not incomplete:
        print "Please check network and provide more --inputs."
    else:
        print "Possible --inputs to add:"
        seen = {}
        index = 1
        for x in incomplete:
            path, start_ids, missing_ids = x

            # Convert IDs to datatype names.
            #x = start_ids
            #x = [network.nodes[x].datatype.name for x in x]
            #x = tuple(sorted(x))
            #start_names = x
            x = missing_ids
            x = [network.nodes[x].datatype.name for x in x]
            x = tuple(sorted(x))
            missing_names = x

            if missing_names in seen:
                continue
            seen[missing_names] = 1

            #x1 = ", ".join(start_names)
            x = ", ".join(missing_names)
            print "%d.  %s" % (index, x)
            index += 1
        print
    sys.stdout.flush()

    # Make a list of all possible start_ids.
    start_ids = []
    for i, x in enumerate(incomplete):
        path, sids, mids = x
        start_ids.extend(sids)
    start_ids = _uniq(start_ids)
    plot_network(
        network_png, network, user_options=user_options,
        highlight_green=start_ids, verbose=verbose)
    return []


def check_inputs_complete(
    network, user_options, in_data_nodes, paths, network_png, verbose):
    # Return list of 
    # paths is a list of (path, start_ids).
    import sys
    from Betsy import bie3

    print "Assigning --input's to nodes in the network."
    sys.stdout.flush()

    data_nodes = [x.data for x in in_data_nodes]
    # index of data_nodes -> list of compatible node_ids
    #node2ids = [bie3._find_start_nodes(network, x) for x in data_nodes]

    # For each path, map data_nodes to start_ids.
    complete = []    # list of path, start_ids, data_indexes
    incomplete = []  # list of path, start_ids, data_indexes
    for x in paths:
        path, start_ids = x
        assert len(start_ids) <= len(in_data_nodes)
        x = _assign_data_nodes_to_start_ids(network, data_nodes, start_ids)
        for x in x:
            score, data_indexes = x
            x = path, start_ids, data_indexes
            if score == 0:
                # May have duplicates if user gives redundant inputs.
                if x not in complete:
                    complete.append(x)
            else:
                if x not in incomplete:
                    incomplete.append(x)

    if complete:
        print "%d pipelines can use these inputs." % len(complete)
        return complete
    print "I cannot find a pipeline given the current --input's."

    if not incomplete:
        # No good assignments.
        print "Cannot find inputs in the network with the right data types."
        print "Please verify the data types carefully."
    else:
        pipeline2scores = []
        for x in incomplete:
            path, start_ids, data_indexes = x
            assert len(start_ids) == len(data_indexes)
            scores = []
            for j in range(len(start_ids)):
                data_index = data_indexes[j]
                start_id = start_ids[j]
                x = bie3._score_start_nodes(
                    network, data_nodes[data_index], good_ids=[start_id])
                scores.extend(x)
            pipeline2scores.append(scores)

        # Print out incompatibilities seen in all pipelines.
        _print_global_incompatibilities(network, pipeline2scores)
        
        # Print out incompatibilities between the data nodes and
        # network nodes for each pipeline.
        print
        for i, scores in enumerate(pipeline2scores):
            print "Pipeline %d." % (i+1)
            _print_node_score_table(network, scores)
            print

    # Make a list of all possible start_ids.
    start_ids = []
    for i, x in enumerate(incomplete):
        path, sids, inds = x
        start_ids.extend(sids)
    start_ids = _uniq(start_ids)
    plot_network(
        network_png, network, user_options=user_options,
        highlight_green=start_ids, verbose=verbose)
    return []


def _assign_data_nodes_to_start_ids(network, data_nodes, start_ids):
    # Try every combination and return the one with the lowest score.
    # Return list of (score, list of indexes of data_nodes parallel to
    # start_ids).  Low scores is better match.  0 is exact match.
    import itertools
    from Betsy import bie3

    assert len(start_ids) <= len(data_nodes)
    data_indexes = range(len(data_nodes))

    # Score each of the data nodes against each of the start_ids.
    score_cache = {}   # (data_index, start_id) -> score
    for index in data_indexes:
        scores = bie3._score_start_nodes(
            network, data_nodes[index], good_ids=start_ids)
        for start_id in start_ids:
            x = [x for x in scores if x[1] == start_id]
            x = [x[0] for x in x]
            if not x:  # wrong data type
                continue
            score_cache[(index, start_id)] = min(x)

    scored_combos = []
    for x in itertools.permutations(data_indexes):
        indexes = x[:len(start_ids)]  # assignment

        compatible = True
        total = 0
        for i in range(len(indexes)):
            data_index = indexes[i]
            start_id = start_ids[i]
            if (data_index, start_id) not in score_cache:
                compatible = False
                break
            score = score_cache[(data_index, start_id)]
            total += score
        if not compatible:  # mismatched data types
            continue
        x = total, indexes
        scored_combos.append(x)
    return scored_combos


def _print_global_incompatibilities(network, pipeline2scores):
    # If there is an incompatibility seen across all pipelines, print
    # it out.
    # scores from bie3._score_start_nodes

    incompatibility = {}  # (datatype_name, attr_name, network_value) -> count
    for scores in pipeline2scores:
        seen = {}  # only count once for each pipeline
        for x in scores:
            score, node_id, user_data, attr_values = x
            dt_name = network.nodes[node_id].datatype.name
            for name, netw_value, user_value in attr_values:
                if type(netw_value) is type([]):
                    netw_value = tuple(netw_value)
                x = dt_name, name, netw_value
                if x in seen:
                    continue
                seen[x] = 1
                incompatibility[x] = incompatibility.get(x, 0) + 1
    printed = False
    for (x, count) in incompatibility.iteritems():
        if count < len(pipeline2scores):
            continue
        dt_name, attr_name, netw_value = x
        prefix = ""
        if not printed:
            prefix = "Given these inputs, "
            printed = True
        print "%s%s:%s needs to be %s." % (
            prefix, dt_name, attr_name, netw_value)


def _print_node_score_table(network, scores):
    # scores from bie3._score_start_nodes
    
    # Make an output table.
    table = []
    header = [
        "Node", "D", "Datatype", "Attribute", "Yours", "Network"]
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


def check_attributes_complete(
    network, user_options, paths, network_png, verbose):
    import sys
    from Betsy import bie3

    # Make sure all required mattr are provided.  This can only be run
    # after the final network is generated.
    print "Making sure all required attributes (--mattr) are provided."
    sys.stdout.flush()

    assert paths
    
    all_missing = {}
    good_paths = []
    for x in paths:
        path, start_ids, data_indexes = x
        assert len(start_ids) == len(data_indexes)

        module_ids = [
            x for x in path if isinstance(network.nodes[x], bie3.ModuleNode)]
        modules = [network.nodes[x] for x in module_ids]

        missing = {}
        opt2mods = get_required_option_names(modules)
        x = opt2mods.keys()
        x = [x for x in x if x not in user_options]
        for on in x:
            for mn in opt2mods[on]:
                missing[(mn, on)] = 1
        if not missing:
            x = path, start_ids, data_indexes
            good_paths.append(x)
        all_missing.update(missing)

    if good_paths:
        return good_paths
    for (mn, on) in all_missing:
        print 'Missing --mattr: %s requires attribute "%s".' % (mn, on)

    start_ids = []
    for x in paths:
        path, sids, inds = x
        start_ids.extend(sids)
    start_ids = _uniq(start_ids)
    plot_network(
        network_png, network, user_options=user_options,
        highlight_green=start_ids, verbose=verbose)
    return []


def prune_pipelines(network, user_options, paths, verbose):
    # Any any pipelines look weird, then remove them.
    import sys

    print "Pruning redundant pipelines."
    sys.stdout.flush()
    num_paths_orig = len(paths)

    if 0:  # DEBUG
        for i, x in enumerate(paths):
            path, start_ids, data_indexes = x
            filename = "pipeline-%02d.png" % i
            plot_network(
                filename, network, user_options=user_options,
                bold=path, highlight_green=start_ids, verbose=verbose)
    paths = _prune_redundant_pipelines(network, paths)
    assert paths, "All pipelines pruned.  This should not happen."
    num_pruned = num_paths_orig - len(paths)
    if not num_pruned:
        print "No redundant pipelines found."
    else:
        print "Pruned %d pipelines.  %d left." % (num_pruned, len(paths))
    return paths


def _prune_redundant_pipelines(network, paths):
    # Remove pipelines that are just supersets of another pipeline.
    import itertools

    redundant = []
    for (i, j) in itertools.product(range(len(paths)), range(len(paths))):
        if i == j:
            continue
        if _is_redundant_pipeline(network, paths[i], paths[j]):
            redundant.append(i)
    redundant = sorted({}.fromkeys(redundant).keys())
    redundant = reversed(redundant)
    for i in redundant:
        del paths[i]
    return paths


def _is_redundant_pipeline(network, path_1, path_2):
    # Test if path_1 is redundant, given path_2.  I.e. path_1 has more
    # processing steps that are not necessary in path_2.
    node_ids_1, start_ids_1, data_indexes_1 = path_1
    node_ids_2, start_ids_2, data_indexes_2 = path_2
    
    # If they're redundant, then the start nodes must be different.
    if sorted(start_ids_1) == sorted(start_ids_2):
        return False
    # If they're redundant, path_1 must be a superset of path_2.
    if len(node_ids_1) <= len(node_ids_2):
        return False
    if not set(node_ids_1).issuperset(node_ids_2):
        return False
    # At least one of the extra nodes must be a start node.
    x = list(set(node_ids_1).difference(node_ids_2))
    x = [x for x in x if x in start_ids_1]
    if not x:
        return False
    super_start_ids_1 = x

    # All paths starting from super_start_ids_1 must go through one of
    # start_ids_2.
    for start_id in super_start_ids_1:
        if not _does_path_go_through(network, start_id, start_ids_2):
            return False
    return True


def _does_path_go_through(network, start_id, intermediate_ids):
    # Do all paths starting from start_id go through intermediate_ids?
    stack = [start_id]
    while stack:
        node_id = stack.pop(0)
        if node_id in intermediate_ids:
            continue
        next_ids = network.transitions.get(node_id, [])
        if not next_ids:
            # Goes to end of network.  Did not encounter intermediate_ids.
            return False
        stack.extend(next_ids)
    return True


def check_input_files(in_data_nodes):
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

    if missing:
        return False
    #print "all input files found."
    #sys.stdout.flush()
    return True


def manually_verify_network(network, user_options, paths, run, network_png,
                            verbose):
    import sys
    if run:
        return True

    print "Please review the network to make sure the analysis is correct."
    print "Add --run when ready to run the analysis."
    sys.stdout.flush()

    path_ids = []
    start_ids = []
    for x in paths:
        path, sids, inds = x
        path_ids.extend(path)
        start_ids.extend(sids)
    path_ids = _uniq(path_ids)
    start_ids = _uniq(start_ids)
    plot_network(
        network_png, network, user_options=user_options,
        bold=path_ids, highlight_green=start_ids, verbose=verbose)
    return False
    

## def print_input_nodes(
##     network, required_datatypes, excluded_datatypes, user_attributes,
##     max_inputs=None, outhandle=None):
##     # required_datatypes is a list of the names of datatypes that must
##     # be in this combination.
##     import sys
##     from genomicode import parselib
##     from Betsy import bie3

##     outhandle = outhandle or sys.stdout

##     ps = parselib.print_split

##     print >>outhandle, "Possible Inputs"
##     outhandle.flush()

##     inputs = bie3.get_input_nodes(
##         network, user_attributes, skip_datatypes=excluded_datatypes,
##         skip_private_datatypes=True, max_inputs=max_inputs)
##     dt2inputs = bie3.group_nodes_by_datatype(network, inputs)

##     # Make a list of the datatypes to show.
##     datatype_combos = dt2inputs.keys()  # list of datatypes

##     ## # Remove any that contain "private" datatypes.
##     ## i = 0
##     ## while i < len(datatype_combos):
##     ##     has_private = False
##     ##     for x in datatype_combos[i]:
##     ##         if x.name.startswith("_"):
##     ##             has_private = True
##     ##             break
##     ##     if has_private:
##     ##         del datatype_combos[i]
##     ##     else:
##     ##         i += 1


##     # Remove any that do not contain all of the required_datatypes.
##     # Also remove those that contain datatypes not in
##     # required_datatypes.
##     s_required_datatypes = sorted(_uniq(required_datatypes))
##     i = 0
##     while i < len(datatype_combos):
##         combo = datatype_combos[i]
##         names = [x.name for x in combo]
##         x = sorted(_uniq(names))
##         if x != s_required_datatypes:
##             del datatype_combos[i]
##         else:
##             i += 1
##     # Sort by number of datatypes, then names.
##     schwartz = []
##     for dt in datatype_combos:
##         x1 = len(dt)
##         x2 = [x.name for x in dt]
##         schwartz.append((x1, x2, dt))
##     schwartz.sort()
##     datatype_combos = [x[-1] for x in schwartz]

##     if not datatype_combos:
##         print >>outhandle, "None"
##         return
##     for i, combo in enumerate(datatype_combos):
##         x = [x.name for x in combo]
##         ps("%d.  INPUTS: %s" % (i+1, ", ".join(x)), outhandle=outhandle)

##         # node_id_combos a list of list of node_ids.  node_combos is
##         # a list of list of nodes, all with the same combination of
##         # types.
##         # E.g. [
##         #   [AgilentFiles, ClassLabelFile],
##         #   [AgilentFiles, ClassLabelFile],
##         # ]
##         node_id_combos = dt2inputs[combo]
##         node_combos = [None] * len(node_id_combos)
##         for i, id_combo in enumerate(node_id_combos):
##             # id_combo is a list of node_ids.  node_combo is a list of
##             # DataNodes.
##             node_combo = [network.nodes[x] for x in id_combo]
##             # Sort the DataNodes by datatype name, then attributes.
##             schwartz = []
##             for id_, n in zip(id_combo, node_combo):
##                 x1 = n.datatype.name
##                 x2 = n.attributes
##                 schwartz.append((x1, x2, id_, n))
##             schwartz.sort()
##             node_id_combos[i] = [x[-2] for x in schwartz]
##             node_combos[i] = [x[-1] for x in schwartz]

##         # DEBUG: Make sure no IDs occur more than once
##         #for x in node_id_combos:
##         #    y = {}.fromkeys(x).keys()
##         #    assert len(x) == len(y)

##         # Merge the node_combos with very similar attributes.  Do this
##         # with a greedy algorithm.  Take each node_combo, and make a
##         # group of other node_combos to merge.  Every member of the
##         # node_combo must vary by at most one attribute (and must be
##         # the same attribute).
##         j = 0
##         while j < len(node_combos)-1:
##             # Indexes of node_combos to merge, not including this
##             # combo.
##             group = []

##             # Find other nodes that fit into this group.
##             diff = None
##             for k in range(j+1, len(node_combos)):
##                 combo1, combo2 = node_combos[j], node_combos[k]
##                 x = _list_differences_in_nodelists(combo1, combo2)
##                 if len(x) > 1:
##                     continue
##                 if not x:
##                     # Sometimes can have no differences.  Just merge
##                     # them.
##                     pass
##                 elif diff is None:
##                     diff = x
##                 elif x != diff:
##                     continue
##                 group.append(k)
##             if not group:
##                 j += 1
##                 continue
##             # Merge every combo in this group, and delete the unmerged
##             # combinations.
##             node_combo = node_combos[j]
##             for k in group:
##                 node_combo = _merge_nodelists(node_combo, node_combos[k])
##             node_combos[j] = node_combo
##             node_combos = [
##                 node_combos[k] for k in range(len(node_combos))
##                 if k not in group]

##         # Print out each combination.
##         for j, node_combo in enumerate(node_combos):
##             # Print one node.
##             for k, node in enumerate(node_combo):
##                 assert isinstance(node, bie3.DataNode)
##                 #x = "- %s: %s" % (node.datatype.name, node.datatype.help)
##                 x = "- %s" % node.datatype.name
##                 ps(x, prefix1=4, prefixn=8, outhandle=outhandle)

##                 # Print the attributes.
##                 for name in sorted(node.attributes):
##                     value = node.attributes[name]
##                     # To save space, don't print the ones with default
##                     # values.
##                     attrdef = node.datatype.attribute_defs[name]
##                     if attrdef.default_in == value:
##                         continue
##                     #if type(value) is type([]) and \
##                     #       attrdef.default_in in value:
##                     #    continue
##                     x = "%s=%s" % (name, value)
##                     ps(x, prefix1=8, prefixn=10, outhandle=outhandle)
##             print >>outhandle  # space between each combination.
##             outhandle.flush()
##         print >>outhandle      # two spaces at end of this set of combos.


## def print_diagnose(network, start_nodes, missing, outhandle=None):
##     # start_nodes is a list of DataNode objects.  missing is a list of
##     # the indexes of start_nodes that cannot be found in the network.
##     from genomicode import parselib
##     from Betsy import bie3

##     x = [start_nodes[i] for i in missing]
##     x = [str(x) for x in x]
##     node_or_nodes = "node"
##     if len(x) > 1:
##         node_or_nodes = "nodes"
##     x = "Input %s not found:\n%s" % (node_or_nodes, "\n".join(x))
##     parselib.print_split(x, outhandle=outhandle)
##     parselib.print_split(
##         "Possible reasons are:", outhandle=outhandle)
##     parselib.print_split(
##         "- There is an incompatibility in the attributes somewhere.",
##         prefix1=4, outhandle=outhandle)
##     parselib.print_split(
##         "- It is not needed for this pipeline.", prefix1=4,
##         outhandle=outhandle)
##     print >>outhandle
    
##     #bie3.diagnose_start_node(network, start_nodes, outhandle=outhandle)
##     results = bie3._score_start_nodes(network, start_nodes)

##     # This can happen if the network doesn't contain any of the same
##     # datatypes as start_nodes.  This can be if the network is really
##     # small, e.g. if there's no way at all to generate the desired out
##     # data type.
##     if not results:
##         node_or_nodes = "node"
##         if len(network.nodes) > 1:
##             node_or_nodes = "nodes"
##         print >>outhandle, \
##               "The network (%d %s) contains no matching data types." % (
##             len(network.nodes), node_or_nodes)
##         return

##     # Make an output table.
##     table = []
##     header = [
##         "Node", "D", "Datatype", "Attribute", "Network", "User"]
##     table.append(header)
##     for x in results:
##         score, node_id, user_data, attr_values = x
##         dt_name = network.nodes[node_id].datatype.name
##         if not attr_values:
##             x = node_id, score, dt_name, "", "", ""
##             assert len(x) == len(header)
##             table.append(x)
##         for name, netw_value, user_value in attr_values:
##             x = node_id, score, dt_name, name, netw_value, user_value
##             assert len(x) == len(header)
##             table.append(x)
    
##     # Figure out the maximum lengths of each column.
##     num_cols = len(header)
##     col_lengths = []
##     for i in range(num_cols):
##         x = [x[i] for x in table]     # Get values in column.
##         x = [len(str(x)) for x in x]  # calculate lengths.
##         x = max(x)
##         col_lengths.append(x)

##     # Set a maximum limit for Datatype and Attribute columns.
##     # Should be max 79 columns long, including 5 for spaces.
##     # Column      MIN  MAX  Notes
##     # Node         4     4  As short as possible.
##     # Delta        1     5  As short as possible.
##     # Datatype     8    18  Datatype name.
##     # Attribute    9    19  Attribute name.
##     # Network      7    22  Might be long.  Lots of values.
##     # User         4    10  Can be long.  But usually short.
##     max_lengths = [4, 1, 18, 19, 22, 10]
##     assert len(col_lengths) == len(max_lengths)
##     col_lengths = max_lengths
##     #for i in range(len(col_lengths)):
##     #    col_lengths[i] = min(col_lengths[i], max_lengths[i])
##     # Just use the maximum lengths.
##     # Make sure the values aren't too long.
##     for i in range(len(table)):
##         row = list(table[i])
##         for j in range(len(row)):
##             x = row[j]
##             x = str(x)
##             x = x.rstrip()
##             if len(x) > col_lengths[j]:
##                 x = x[:col_lengths[j]-3] + "..."
##             row[j] = x
##         table[i] = row

##     fmt = "{!s:^%ds} {!s:^%ds} {:<%ds} {:<%ds} {:<%ds} {:<%ds}" % \
##           tuple(col_lengths)
##     for x in table:
##         print >>outhandle, fmt.format(*x)


def plot_network(filename, network, user_options=None, bold=[],
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
        highlight_green=highlight_green, highlight_orange=highlight_orange,
        highlight_purple=highlight_purple, highlight_yellow=highlight_yellow,
        verbose=verbose)


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


## def _delete_input_nodes(
##     network, wanted_dt_names, unwanted_dt_names, user_attributes):
##     # Only one of wanted_dt_names and unwanted_dt_names should be given.
##     from Betsy import bie3

##     assert not (wanted_dt_names and unwanted_dt_names)

##     wanted_dt_names = wanted_dt_names or []
##     unwanted_dt_names = unwanted_dt_names or []
##     wanted_dict = {}.fromkeys(wanted_dt_names)
##     unwanted_dict = {}.fromkeys(unwanted_dt_names)

##     while True:
##         # Keep all the internal IDs.  These are the nodes that someone
##         # points to.
##         keep_ids = {}
##         for next_ids in network.transitions.values():
##             x = {}.fromkeys(next_ids)
##             keep_ids.update(x)
##         # Everything else is a start ID.
##         start_ids = [x for x in range(len(network.nodes)) if x not in keep_ids]
##         if wanted_dict:
##             # Keep all the start_ids that is in a wanted datatype.
##             x = [
##                 x for x in start_ids
##                 if network.nodes[x].datatype.name in wanted_dict]
##         elif unwanted_dict:
##             # Keep all the start_ids that is not in the unwanted
##             # datatype.
##             x = [
##                 x for x in start_ids
##                 if network.nodes[x].datatype.name not in unwanted_dict]
##         x = {}.fromkeys(x)
##         keep_ids.update(x)

##         # Remove every node not in keep_ids.
##         remove_ids = [
##             x for x in range(len(network.nodes)) if x not in keep_ids]
##         if not remove_ids:
##             break
##         network = network.delete_nodes(remove_ids)
##         network = bie3._OptimizeNoDanglingNodes().optimize(
##             network, user_attributes)
##         assert network.nodes, "Empty network"
##     network = bie3.optimize_network(network, user_attributes)
##     return network


## def _remove_unwanted_input_nodes(network, unwanted_datatypes, user_attributes):
##     # Remove all the input nodes that aren't one of these wanted
##     # datatypes.
##     # unwanted_datatypes is a list of names.
##     return _delete_input_nodes(
##         network, None, unwanted_datatypes, user_attributes)


## def _keep_wanted_input_nodes(network, wanted_datatypes, user_attributes):
##     # Keep only the input nodes that represents one of these wanted
##     # datatypes.
##     # wanted_datatypes is a list of names.
##     return _delete_input_nodes(
##         network, wanted_datatypes, None, user_attributes)


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


def _uniq(seq):
    return {}.fromkeys(seq).keys()


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
    group.add_argument(
        "--exclude_input", action="append", help="Experimental.")
    group.add_argument(
        "--dont_find_mattr_files", action="store_true",
        help="By default, if an --mattr option looks like a filename, "
        "will convert to full path so that modules can find them.  "
        "This will turn off that feature.")
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
        "given after the --output.  Format: <datatype>.<key>=<value>.")
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
    # TODO: Make sure args.exclude_input is valid.
    # args.exclude_input


    # Make sure configuration directory exists.
    if not os.path.exists(config.OUTPUTPATH):
        print "Making BETSY working path: %s." % config.OUTPUTPATH
        os.mkdir(config.OUTPUTPATH)

    # Make sure files exist.
    #if input_list or args.output_file:
    #    print "Looking for files."
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
    #if input_list:
    #    print "Parsing input options."
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

    # test outtype and build the list of user_attributes.
    user_attributes = []  # List of bie3.Attribute objects.
    if outtype:
        #print "Parsing output options."

        # Get the user_attributes.  These are the attributes from the
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
            user_attributes.append(bie3.Attribute(fn, key, value))

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
        assert '=' in x, "--mattr should be in format: <option>=<value>"
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
    

    # Cases:   OBSOLETE
    # 1.  No inputs and no outputs.
    #     Show the whole rulebase.
    # 2.  No inputs and has output.  (requires network)
    #     Show the input datatypes.
    # 3.  Input and no outputs.
    #     ERROR
    # 4.  Input and output.  (requires network)
    #     Show the input nodes.
    #     Run the analysis.

    # Step 1: Make sure there's an output provided.
    #         If not, show the rulebase.
    # Step 2: Generate network that produces this output.
    # Step 3: Make sure there are inputs provided.
    #         If not, show the input datatypes.
    # Step 4: Make sure the proper data types are provided.
    # Step 5: Make sure the inputs can generate a path through the
    #         network.
    # Step 6: Make sure required attributes are given.
    # Step 7: Prune bad pipelines.
    # Step 8: Make sure input files exist.
    # Step 9: Manual verification of the network by the user.

    # Step 1: Make sure an output is provided.
    if not check_output_provided(rulebase, args.output):
        return
    # Step 2: Generate network.
    network = generate_network(
        rulebase, outtype, user_attributes, args.exclude_input)
    #network_orig = copy.deepcopy(network)
    # Step 3: Make sure inputs are provided.
    if not check_inputs_provided(
        network, args.exclude_input, user_attributes, user_options,
        in_data_nodes, args.max_inputs, args.network_png, verbose):
        return
    # Step 4: Make sure the data types provided will generate a viable
    # path through the network.
    paths = check_datatypes_complete(
        network, user_attributes, user_options, in_data_nodes,
        args.network_png, verbose)
    if not paths:
        return
    # Step 5: Make sure the inputs can generate a path through the
    # network.
    paths = check_inputs_complete(
        network, user_options, in_data_nodes, paths, args.network_png, verbose)
    if not paths:
        return
    # Step 6: Make sure required attributes are given.
    paths = check_attributes_complete(
        network, user_options, paths, args.network_png, verbose)
    if not paths:
        return
    # Step 7: Prune bad pipelines.
    paths = prune_pipelines(network, user_options, paths, verbose)
    if not paths:
        return
    # Step 8: Look for input files.
    if not check_input_files(in_data_nodes):
        return
    #print "The network (%d nodes) is complete." % len(network.nodes)
    #print "There are %d possible pipelines." % len(paths)
    print "Ready to go!"
    sys.stdout.flush()
    # Step 9: Manual verification of the network.
    if not manually_verify_network(
        network, user_options, paths, args.run, args.network_png, verbose):
        return

    # DEBUG: Print out each of the paths.
    if 0:
        for i, x in enumerate(paths):
            path, start_ids, data_indexes = x
            filename = "path%03d.png" % i
            plot_network(
                filename, network, user_options=user_options,
                bold=path, highlight_purple=path, highlight_green=start_ids,
                verbose=verbose)
        

    print "Running the analysis."
    sys.stdout.flush()
    clean_up = not args.save_failed_data
    node_dict = output_file = None
    try:
        x = rule_engine_bie3.run_pipeline(
            network, in_data_nodes, user_attributes, user_options, paths,
            user=args.user, job_name=args.job_name, clean_up=clean_up,
            num_cores=args.num_cores)
        if x:
            node_dict, output_file = x
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
        bold=node_ids, highlight_green=node_ids, verbose=verbose)
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


if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
