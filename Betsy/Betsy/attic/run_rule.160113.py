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
        #path, start_ids, missing_ids = x
        used = [i for (i, x) in enumerate(p.start_ids) if x is not None]
        inputs_used.update({}.fromkeys(used))
    has_unused_inputs = len(inputs_used) != len(in_data_nodes)

    ## # Filter for no missing IDs.  Also, convert the start_ids to
    ## # start_ids and data_indexes.
    ## good_paths = []
    ## for p in paths:
    ##     #path, start_ids, missing_ids = x
    ##     if p.missing_ids:
    ##         continue
    ##     sids = []
    ##     indexes = []
    ##     for i, sid in enumerate(p.start_ids):
    ##         if sid is None:
    ##             continue
    ##         sids.append(sid)
    ##         indexes.append(i)
    ##     x = path, sids, indexes
    ##     good_paths.append(x)
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
        
        print
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


## def check_datatypes_complete(
##     network, user_attributes, user_options, in_data_nodes, network_png,
##     verbose):
##     # Return list of (path, start_ids).  Empty list means no paths
##     # found.
##     import sys
##     from Betsy import bie3

##     print "Searching the network for pipelines that use --input data types."
##     sys.stdout.flush()
##     # Make a list of the names of data types for start nodes.
##     assert in_data_nodes
##     start_dtype_names = [x.data.datatype.name for x in in_data_nodes]

##     complete = []    # list of (path, start, missing)
##     incomplete = []  # list of (path, start, missing)
##     for i, x in enumerate(bie3.find_paths_by_datatypes(
##             network, user_attributes, start_dtype_names)):
##         path, start_ids, missing_ids = x

##         # Ignore paths that don't contain any of the requested data types.
##         if not start_ids:
##             continue
##         # Ignore paths with too many missing data types.
##         if len(missing_ids) > 5:
##             continue
        
##         x = path, start_ids, missing_ids
##         if not missing_ids:
##             assert len(start_ids) <= len(in_data_nodes)
##             complete.append(x)
##         else:
##             incomplete.append(x)

##         ## # DEBUG: print out interesting paths.
##         ## if 0 and len(missing_ids) == 0:
##         ##     filename = "path-%d.png" % i
##         ##     plot_network(
##         ##         filename, network, user_options=user_options,
##         ##         bold=path, highlight_green=start_ids,
##         ##         highlight_purple=missing_ids, verbose=True)

##     if complete:
##         print "Found %d possible pipelines." % len(complete)
##         x = [tuple(x[:2]) for x in complete]
##         return x

##     # Print a diagnostic message.
##     print "No pipelines found."
##     if not incomplete:
##         print "Please check network and provide more --inputs."
##     else:
##         print "Possible --inputs to add:"
##         seen = {}
##         index = 1
##         for x in incomplete:
##             path, start_ids, missing_ids = x

##             # Convert IDs to datatype names.
##             #x = start_ids
##             #x = [network.nodes[x].datatype.name for x in x]
##             #x = tuple(sorted(x))
##             #start_names = x
##             x = missing_ids
##             x = [network.nodes[x].datatype.name for x in x]
##             x = tuple(sorted(x))
##             missing_names = x

##             if missing_names in seen:
##                 continue
##             seen[missing_names] = 1

##             #x1 = ", ".join(start_names)
##             x = ", ".join(missing_names)
##             print "%d.  %s" % (index, x)
##             index += 1
##         print
##     sys.stdout.flush()

##     # Make a list of all possible start_ids.
##     start_ids = []
##     for i, x in enumerate(incomplete):
##         path, sids, mids = x
##         start_ids.extend(sids)
##     start_ids = _uniq(start_ids)
##     plot_network(
##         network_png, network, user_options=user_options,
##         highlight_green=start_ids, verbose=verbose)
##     return []


## def check_inputs_complete(
##     network, user_options, in_data_nodes, paths, network_png, verbose):
##     # Return list of 
##     # paths is a list of (path, start_ids).
##     import sys
##     from Betsy import bie3

##     print "Assigning --input's to nodes in the network."
##     sys.stdout.flush()

##     data_nodes = [x.data for x in in_data_nodes]
##     # index of data_nodes -> list of compatible node_ids
##     #node2ids = [bie3._find_start_nodes(network, x) for x in data_nodes]

##     # For each path, map data_nodes to start_ids.
##     complete = []    # list of path, start_ids, data_indexes
##     incomplete = []  # list of path, start_ids, data_indexes
##     for x in paths:
##         path, start_ids = x
##         assert len(start_ids) <= len(in_data_nodes)
##         x = _assign_data_nodes_to_start_ids(network, data_nodes, start_ids)
##         for x in x:
##             score, data_indexes = x
##             x = path, start_ids, data_indexes
##             if score == 0:
##                 # May have duplicates if user gives redundant inputs.
##                 if x not in complete:
##                     complete.append(x)
##             else:
##                 if x not in incomplete:
##                     incomplete.append(x)

##     if complete:
##         print "%d pipelines can use these inputs." % len(complete)
##         return complete
##     print "I cannot find a pipeline given the current --input's."

##     if not incomplete:
##         # No good assignments.
##         print "Cannot find inputs in the network with the right data types."
##         print "Please verify the data types carefully."
##     else:
##         pipeline2scores = []
##         for x in incomplete:
##             path, start_ids, data_indexes = x
##             assert len(start_ids) == len(data_indexes)
##             scores = []
##             for j in range(len(start_ids)):
##                 data_index = data_indexes[j]
##                 start_id = start_ids[j]
##                 x = bie3._score_start_nodes(
##                     network, data_nodes[data_index], good_ids=[start_id])
##                 scores.extend(x)
##             pipeline2scores.append(scores)

##         # Print out incompatibilities seen in all pipelines.
##         _print_global_incompatibilities(network, pipeline2scores)
        
##         # Print out incompatibilities between the data nodes and
##         # network nodes for each pipeline.
##         print
##         for i, scores in enumerate(pipeline2scores):
##             print "Pipeline %d." % (i+1)
##             _print_node_score_table(network, scores)
##             print

##     # Make a list of all possible start_ids.
##     start_ids = []
##     for i, x in enumerate(incomplete):
##         path, sids, inds = x
##         start_ids.extend(sids)
##     start_ids = _uniq(start_ids)
##     plot_network(
##         network_png, network, user_options=user_options,
##         highlight_green=start_ids, verbose=verbose)
##     return []


## def _assign_data_nodes_to_start_ids(network, data_nodes, start_ids):
##     # Try every combination and return the one with the lowest score.
##     # Return list of (score, list of indexes of data_nodes parallel to
##     # start_ids).  Low scores is better match.  0 is exact match.
##     import itertools
##     from Betsy import bie3

##     assert len(start_ids) <= len(data_nodes)
##     data_indexes = range(len(data_nodes))

##     # Score each of the data nodes against each of the start_ids.
##     score_cache = {}   # (data_index, start_id) -> score
##     for index in data_indexes:
##         scores = bie3._score_start_nodes(
##             network, data_nodes[index], good_ids=start_ids)
##         for start_id in start_ids:
##             x = [x for x in scores if x[1] == start_id]
##             x = [x[0] for x in x]
##             if not x:  # wrong data type
##                 continue
##             score_cache[(index, start_id)] = min(x)

##     scored_combos = []
##     for x in itertools.permutations(data_indexes):
##         indexes = x[:len(start_ids)]  # assignment

##         compatible = True
##         total = 0
##         for i in range(len(indexes)):
##             data_index = indexes[i]
##             start_id = start_ids[i]
##             if (data_index, start_id) not in score_cache:
##                 compatible = False
##                 break
##             score = score_cache[(data_index, start_id)]
##             total += score
##         if not compatible:  # mismatched data types
##             continue
##         x = total, indexes
##         scored_combos.append(x)
##     return scored_combos


## def _print_global_incompatibilities(network, pipeline2scores):
##     # If there is an incompatibility seen across all pipelines, print
##     # it out.
##     # scores from bie3._score_start_nodes

##     incompatibility = {}  # (datatype_name, attr_name, network_value) -> count
##     for scores in pipeline2scores:
##         seen = {}  # only count once for each pipeline
##         for x in scores:
##             score, node_id, user_data, attr_values = x
##             dt_name = network.nodes[node_id].datatype.name
##             for name, netw_value, user_value in attr_values:
##                 if type(netw_value) is type([]):
##                     netw_value = tuple(netw_value)
##                 x = dt_name, name, netw_value
##                 if x in seen:
##                     continue
##                 seen[x] = 1
##                 incompatibility[x] = incompatibility.get(x, 0) + 1
##     printed = False
##     for (x, count) in incompatibility.iteritems():
##         if count < len(pipeline2scores):
##             continue
##         dt_name, attr_name, netw_value = x
##         prefix = ""
##         if not printed:
##             prefix = "Given these inputs, "
##             printed = True
##         x = netw_value
##         if type(x) is type(""):
##             x = '"%s"' % x
##         # If netw_value is a tuple, then don't do any special formatting.
##         print '%s%s.%s needs to be %s.' % (prefix, dt_name, attr_name, x)


def check_attributes_complete(
    network, user_options, paths, network_png, verbose):
    import sys
    from genomicode import parselib
    from Betsy import bie3

    # Make sure all required mattr are provided.  This can only be run
    # after the final network is generated.
    print "Making sure all required attributes (--mattr) are provided."
    sys.stdout.flush()

    assert paths
    
    all_missing = {}
    good_paths = []
    for p in paths:
        #path, start_ids, data_indexes = x
        #assert len(start_ids) == len(data_indexes)
        module_ids = [
            x for x in p.node_ids if
            isinstance(network.nodes[x], bie3.ModuleNode)]
        modules = [network.nodes[x] for x in module_ids]

        missing = {}
        opt2mods = get_required_option_names(modules)
        x = opt2mods.keys()
        x = [x for x in x if x not in user_options]
        for on in x:
            for mn in opt2mods[on]:
                missing[(mn, on)] = 1
        if not missing:
            good_paths.append(p)
            #x = path, start_ids, data_indexes
            #good_paths.append(x)
        all_missing.update(missing)

    if good_paths:
        if len(good_paths) < len(paths):
            num_removed = len(paths)-len(good_paths)
            assert all_missing
            if len(all_missing) == 1:
                mn, on = all_missing.keys()[0]
                x = ("Pruned %d pathways because --mattr option %r was "
                     "not provided." % (num_removed, on))
            else:
                names = sorted([x[1] for x in all_missing.keys()])
                x = parselib.pretty_list(names)
                x = ("Pruned %d pathways because the following --mattr "
                     "options were not provided: %s" % (num_removed, x))
            parselib.print_split(x, prefixn=2)
        return good_paths
    for (mn, on) in all_missing:
        print 'Missing --mattr: %s requires attribute "%s".' % (mn, on)

    x = [x.start_ids for x in paths]
    x = bie3._uniq(bie3._flatten(x))
    x = [x for x in x if x is not None]
    start_ids = x
    #start_ids = []
    #for x in paths:
    #    path, sids, inds = x
    #    start_ids.extend(sids)
    #start_ids = _uniq(start_ids)
    plot_network(
        network_png, network, user_options=user_options,
        highlight_green=start_ids, verbose=verbose)
    return []


def prune_pipelines(network, custom_attributes, paths):
    # Any any pipelines look weird, then remove them.
    import sys
    from Betsy import bie3

    print "Pruning redundant pipelines."
    sys.stdout.flush()
    num_paths_orig = len(paths)

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
    paths = _prune_by_custom_attributes(
        network, custom_attributes, paths, nodeid2parents)
    paths = _prune_superset_pipelines(network, paths)
    paths = _prune_alternate_attributes(
        network, custom_attributes, paths, nodeid2parents)
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
        return paths
    
    num_pruned = num_paths_orig - len(paths)
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
    # custom_attributes is a list of Attribute objects.
    from Betsy import bie3
    if not custom_attributes:
        return paths
    
    dname2attrs = {}  # datatype name -> attr name -> value
    for x in custom_attributes:
        if x.datatype.name not in dname2attrs:
            dname2attrs[x.datatype.name] = {}
        dname2attrs[x.datatype.name][x.name] = x.value

    # If a module's input datatype is different from the output
    # datatype, then the attributes of that input datatype should
    # match the user attributes.

    # Search through the network for data nodes that might be subject
    # to custom_attributes.
    all_node_ids = {}
    for x in paths:
        #node_ids, start_ids, data_indexes = x
        #all_node_ids.update(node_ids)
        all_node_ids.update(x.node_ids)
    module_ids = [x for x in all_node_ids
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
        # Find the data_node_ids that match a user attribute.
        x = nodeid2parents.get(module_id, [])
        x = [x for x in x if network.nodes[x].datatype.name in dname2attrs]
        data_node_ids.extend(x)
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


def _prune_alternate_attributes(
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
        x = _list_alternate_attributes(network, i, nodeid2parents)
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


def _list_alternate_attributes(network, data_id, nodeid2parents):
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
    
    #node_ids_1, start_ids_1, data_indexes_1 = path_1
    #node_ids_2, start_ids_2, data_indexes_2 = path_2

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


## def _find_alternate_attributes_old(network, nodeid2previds, path_1, path_2):
##     from Betsy import bie3
    
##     node_ids_1, start_ids_1, data_indexes_1 = path_1
##     node_ids_2, start_ids_2, data_indexes_2 = path_2

##     shared_ids = [x for x in node_ids_1 if x in node_ids_2]
##     unique_ids_1 = [x for x in node_ids_1 if x not in node_ids_2]
##     unique_ids_2 = [x for x in node_ids_2 if x not in node_ids_1]

##     # Look for a DataNode that is shared between the two pathways.
##     data_ids = [
##         x for x in shared_ids if isinstance(network.nodes[x], bie3.DataNode)]
##     # Must have at least two parents.
##     data_ids = [
##         x for x in data_ids if len(nodeid2previds.get(x, [])) >= 2]
##     if not data_ids:
##         return None
    
##     # Look for ModuleNodes that are unique for each pathway.
##     module_ids_1 = [x for x in unique_ids_1
##                     if isinstance(network.nodes[x], bie3.ModuleNode)]
##     module_ids_2 = [x for x in unique_ids_2
##                     if isinstance(network.nodes[x], bie3.ModuleNode)]
##     # This is slower.
##     #module_ids_1 = [x for x in node_ids_1 if x not in shared_ids and 
##     #                isinstance(network.nodes[x], bie3.ModuleNode)]
##     #module_ids_2 = [x for x in node_ids_2 if x not in shared_ids and 
##     #                isinstance(network.nodes[x], bie3.ModuleNode)]
##     if not module_ids_1 or not module_ids_2:
##         return None

##     # Check each data_id carefully.
##     for data_id in data_ids:
##         prev_ids = nodeid2previds[data_id]

##         # These parents must be unique to one of the pathways.
##         prev_ids = [
##             x for x in prev_ids if x in module_ids_1 or x in module_ids_2]
##         if len(prev_ids) < 2:
##             continue
        
##         # 2 parents must have Consequences SET_TO to the same
##         # attribute.
##         attr2values = {}  # name -> list of values to set to
##         for prev_id in prev_ids:
##             module_node = network.nodes[prev_id]
##             for cons in module_node.consequences:
##                 if cons.behavior != bie3.SET_TO:
##                     continue
##                 n, v = cons.name, cons.arg1
##                 assert v is not None
##                 if n not in attr2values:
##                     attr2values[n] = []
##                 if v not in attr2values[n]:
##                     attr2values[n].append(v)
##         # List of attributes with at least 2 values.
##         attrs = [n for (n, v) in attr2values.iteritems() if len(v) >= 2]
##         if not attrs:
##             continue

##         # Look for an attribute with parents from different pathways.
##         for name in attrs:
##             # Look for the modules that SET_TO this attribute.
##             good_prev_ids = []
##             for prev_id in prev_ids:
##                 module_node = network.nodes[prev_id]
##                 found = False
##                 for cons in module_node.consequences:
##                     if cons.behavior != bie3.SET_TO:
##                         continue
##                     if cons.name != name:
##                         continue
##                     found = True
##                     break
##                 if found:
##                     good_prev_ids.append(prev_id)
##             good_1 = [x for x in good_prev_ids if x in module_ids_1]
##             good_2 = [x for x in good_prev_ids if x in module_ids_2]
##             if not good_1 or not good_2:
##                 continue
            
##             # Found one!
##             module_id_1 = good_1[0]
##             module_id_2 = good_2[0]
##             # Get the values of the attributes set by each module.
##             module_1 = network.nodes[module_id_1]
##             module_2 = network.nodes[module_id_2]
##             x1 = [x for x in module_1.consequences if x.name == name]
##             x2 = [x for x in module_2.consequences if x.name == name]
##             assert len(x1) == 1
##             assert len(x2) == 1
##             attr_value_1 = x1[0].arg1
##             attr_value_2 = x2[0].arg1
##             x = name, module_id_1, attr_value_1, \
##                 module_id_2, attr_value_2, data_id
##             return x
##     return None


def _print_attributes_pruned(pruned_attributes):
    # Print out a message about what was deleted.
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
        all_values = sorted(_uniq(attr2values[name]))
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


def _prune_superset_pipelines(network, paths):
    # Remove pipelines that are just supersets of another pipeline.
    import itertools

    #path2length = [None] * len(paths)
    #for i, x in enumerate(paths):
    #    node_ids, start_ids, data_indexes = x
    #    path2length[i] = len(node_ids)
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
    
    #node_ids_1, start_ids_1, data_indexes_1 = path_1
    #node_ids_2, start_ids_2, data_indexes_2 = path_2

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
    # BAM (no) -> addgroup -> BAM (group=y) -> sort -> BAM (sort=y, group=y)
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
    # 
    # Parallel:
    # DataNode -> sort_coord -> mark_dup -> add_read_group -> DataNode
    #          -> add_read_group -> sort_coord -> mark_dup ->
    from Betsy import bie3
    
    ## If they're parallel, they should have the same number of nodes.
    ## The number of parallel nodes should be the same.
    # No.  One of them could have extra steps.

    path_1, path_2 = paths[path_id_1], paths[path_id_2]
    x = _compare_paths(network, path_1, path_2)
    shared_ids, unique_ids_1, unique_ids_2 = x
    
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
    bod1, bod2 = subpath_1.based_on_data, subpath_2.based_on_data
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


# Should optimize this function.
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
    

## def print_input_nodes(
##     network, required_datatypes, excluded_datatypes, custom_attributes,
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


## def is_complete_input(network, user_attributes, start_node_ids):
##     # Return a boolean indicating whether this is a valid set of start
##     # nodes.  (SLOW).
##     from Betsy import bie3

##     start_node_ids = sorted(start_node_ids)
##     inputs = bie3.get_inputs(network, user_attributes)
##     for input_ in inputs:
##         # Bug: should check if input_ is subset of start_node_ids.
##         if start_node_ids == sorted(input_):
##             return True
##     return False
##     #dt2inputs = bie3.group_inputs_by_datatype(network, inputs)
##     #for i, dt in enumerate(sorted(dt2inputs)):
##     #    x = [x.name for x in dt]
##     #    for j, dtinput in enumerate(dt2inputs[dt]):
##     #        required = sorted(list(set(dtinput)))
##     #        if start_node_ids == required:
##     #            return True
##     #return False


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
    # Should merge with bie3.
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
    ## TODO: Make sure args.exclude_input is valid.
    ## args.exclude_input

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

    # test outtype and build the list of custom_attributes.
    custom_attributes = []  # List of bie3.Attribute objects.
    if outtype:
        #print "Parsing output options."

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
    # Step 4: Make sure the inputs can be found in the network.
    # Step 5: Make sure the inputs can generate a pipeline.
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
    if not paths:
        return
    ## # Step 4: Make sure the data types provided will generate a viable
    ## # path through the network.
    ## paths = check_datatypes_complete(
    ##     network, user_attributes, user_options,
    ##     in_data_nodes, args.network_png, verbose)
    ## if not paths:
    ##     return
    ## # Step 5: Make sure the inputs can generate a path through the
    ## # network.
    ## paths = check_inputs_complete(
    ##     network, user_options, in_data_nodes, paths, args.network_png, verbose)
    ## if not paths:
    ##     return
    # Step 6: Make sure required attributes are given.
    paths = check_attributes_complete(
        network, user_options, paths, args.network_png, verbose)
    if not paths:
        return
    # DEBUG: Print out each of the pipelines.
    #plot_pipelines(
    #    "pipeline", network, paths, user_options, max_pipelines=16,
    #    verbose=True)
    # Step 7: Prune undesired pipelines.
    paths = prune_pipelines(network, custom_attributes, paths)
    if not paths:
        return
    # DEBUG: Print out each of the pipelines.
    #plot_pipelines(
    #    "pipeline", network, paths, user_options, max_pipelines=16,
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
            # Show transitions
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
