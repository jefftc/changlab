"""

Functions:
run_pipeline
run_module

"""
# _get_available_input_combinations
# _write_parameter_file

VERSION = 5

BETSY_PARAMETER_FILE = "BETSY_parameters.txt"

DEBUG_POOL = {}
def run_pipeline(
    network, in_datas, custom_attributes, user_options, paths, user=None,
    job_name='', clean_up=True, num_cores=8):
    # Run the pipeline that is indicated by the network.  Returns a
    # tuple of (dictionary of node_id -> IdentifiedDataNode, output
    # filename).  Returns None if not successful.
    # Can raise an exception if there was an error in one of the
    # modules, or if there is no path through the network.
    #
    # in_datas         List of IdentifiedDataNodes.
    # user_attributes  From --dattr.  AttributeDef
    # user_options     From --mattr.  OptionDef
    # paths            List of (node_ids, start_ids).
    global DEBUG_POOL
    
    import os
    import getpass
    import logging
    import time

    from Betsy import bie3
    from Betsy import config
    from Betsy import module_utils

    user = user or getpass.getuser()
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    # Is this thread-safe?
    LOG_FILENAME = os.path.join(output_path, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)

    # Make a list of the valid node_ids in the pipeline.
    x = [x.node_ids for x in paths]
    x = bie3._uniq(bie3._flatten(x))
    path_ids = x    # list of node_ids in the pipeline.

    # Add the start_ids to the stack.  Each member of the stack can be
    # a tuple of:
    # 1.  (IdentifiedDataNode, node_id, None, transitions)
    # 2.  (ModuleNode, node_id, None, transitions)
    # 3.  (ModuleNode, node_id, antecedent_ids, transitions)
    #     Keep track of which set of antecedents to run.
    # transitions is a dictionary of (node_id, next_node_id) -> 1
    stack = []
    for p in paths:
        for index, id_ in enumerate(p.start_ids):
            assert id_ is not None
            node = in_datas[index]
            x = node, id_, None, {}
            if x not in stack:
                stack.append(x)

    # Keep track of nodes that have already been generated.
    pool = {}         # dict of node_id -> IdentifiedDataNode or ModuleNode
    pipeline = []     # list of the names of the modules run.
    next_node = None
    
    # Keep track of the number of times a module can't be run because
    # it doesn't have enough input data.  If this exceeds the total
    # number of modules, then quit.
    num_failures = 0
    success = False
    while stack:
        DEBUG_POOL = pool
        assert num_failures < len(stack), \
               "Inference error: No more nodes to run."
        node, node_id, more_info, transitions = stack.pop()
        if node_id not in path_ids:  # ignore if not in pipeline
            continue
        if node_id == 0:
            pool[node_id] = node
            success = True
            break
        elif isinstance(node, module_utils.IdentifiedDataNode):
            pool[node_id] = node
            assert node_id != 0
            # Add the next modules into the stack, if not already there.
            on_stack = [x[1] for x in stack]
            for next_id in network.transitions[node_id]:
                next_node = network.nodes[next_id]
                assert isinstance(next_node, bie3.ModuleNode)
                if next_id in on_stack:
                    continue
                # Module updates the transitions based on which set of
                # antecedent IDs are used.
                stack.append((next_node, next_id, None, transitions))
            # If I've added new modules, then reset the failures counter.
            num_failures = 0
        elif isinstance(node, bie3.ModuleNode) and more_info is None:
            # If the input data for this module doesn't exist, then
            # just try it again later.
            all_antecedent_ids = _get_available_input_combinations(
                network, node_id, custom_attributes, pool)
            if not all_antecedent_ids:
                # No sets of inputs are ready to run.  Put back to the
                # bottom of the stack and try again later.
                stack.insert(0, (node, node_id, more_info, transitions))
                num_failures += 1
                continue
            for antecedent_ids in all_antecedent_ids:
                assert len(node.in_datatypes) == len(antecedent_ids)
                x = node, node_id, antecedent_ids, transitions
                stack.append(x)
            # If I've added new analyses to run, then reset the
            # failures counter.
            num_failures = 0
        elif isinstance(node, bie3.ModuleNode):
            # Run this module.
            pool[node_id] = node
            pipeline.append(node.name)
            antecedent_ids = more_info
            assert len(node.in_datatypes) == len(antecedent_ids)
            x = run_module(
                network, pipeline, node_id, antecedent_ids,
                user_options, pool, user, job_name, clean_up=clean_up,
                num_cores=num_cores)
            if x is None:
                # Can happen if this module has already been run.  It
                # might've gotten added to the stack because there are
                # many input nodes that can go into this.  This should
                # not count as a failure.
                #num_failures += 1
                continue
            next_node, next_id = x
            # Successfully completed this module.  Reset num_failures.
            num_failures = 0
            #if next_id == 0:
            #    success = True
            #    break
            assert next_id is not None
            trans = transitions.copy()
            for x in antecedent_ids:
                trans[(x, node_id)] = 1
            trans[(node_id, next_id)] = 1
            stack.append((next_node, next_id, None, trans))
        else:
            raise AssertionError

    if success:
        msg = "Completed"
        #print "[%s]  %s" % (time.strftime('%I:%M %p'), msg)
        print "[%s]  %s" % (time.strftime('%a %I:%M %p'), msg)
        #print next_node.identifier
        #sys.stdout.flush()
        assert next_node
        return pool, transitions, next_node.identifier
    print "This pipeline has completed unsuccessfully."
    return None

    #if flag and next_node and module_utils.exists_nz(next_node.identifier):
    #if next_node and module_utils.exists_nz(next_node.identifier):
    #    msg = "Completed successfully and generated a file:"
    #    print "[%s]  %s" % (time.strftime('%I:%M %p'), msg)
    #    print next_node.identifier
    #    sys.stdout.flush()
    #    return next_node.identifier
    #else:
    #    print 'This pipeline has completed unsuccessfully'
    #    raise ValueError('there is no output for this pipeline')
    #return None


def run_module(
    network, pipeline, module_id, input_ids, 
    all_user_options, pool, user, job_name='', clean_up=True, num_cores=8):
    # Return tuple of (IdentifiedDataNode, node_id) for the node that
    # was created.  Returns None if this module fails (no compatible
    # output nodes, or all output nodes already generated).

    import os
    import sys
    import time
    import tempfile
    import shutil
    import logging
    import random

    from genomicode import filelib
    from genomicode import parselib
    from Betsy import config
    from Betsy import module_utils
    from Betsy import bie3

    assert user
    output_path = config.OUTPUTPATH
    assert os.path.exists(output_path), \
           'the output_path %s does not exist' % output_path

    # If there are no possible nodes to generate, return an empty list.
    assert network.transitions[module_id]
    x = network.transitions[module_id]
    x = [x for x in x if x not in pool]
    if not x:
        return None

    # Import module.
    module_node = network.nodes[module_id]
    module_name = module_node.name
    assert len(module_node.in_datatypes) == len(input_ids)
    
    # If module is missing, will raise ImportError with decent error
    # message.
    x = __import__(
        'modules.'+module_name, globals(), locals(), [module_name], -1)
    module = x.Module()

    # Get the antecedents from the pool.
    in_identified_data_nodes = [pool[x] for x in input_ids]
    # Same as in_identified_data_nodes, but lists of length 1 are
    # turned into single objects.
    antecedents = in_identified_data_nodes
    if len(antecedents) == 1:
        antecedents = antecedents[0]

    # Get the user_options.  all_user_options contains all options
    # provided by the user.  Pull out the ones relevant for this
    # module.  Use the defaults when necessary.
    user_options = {}
    for option in module_node.option_defs:
        value = all_user_options.get(option.name)
        if value is None:
            value = option.default
        assert value is not None, "Missing input: %s" % option.name
        user_options[option.name] = value

    # ModuleNodes can point to multiple DataNodes.  They should always
    # be the same type, but their attributes may be different (e.g. is
    # logged or not).  Figure out which one is generated by this
    # module.

    # Create a list of all the output nodes that can be generated by
    # this set of inputs.
    x = [x.data for x in in_identified_data_nodes]
    out_data_nodes = bie3._fc_to_outputs(module_node, x)
    # Make sure the DataTypes are the same.
    assert out_data_nodes
    x = sorted([x.datatype for x in out_data_nodes])
    assert x[0] == x[-1], "ModuleNode points to different DataTypes."

    # Optimization: set_out_attributes can be computationally
    # expensive, so assume that the module will set the output
    # attributes for all the out_nodes identically.  For example, the
    # out_nodes may differ by logged value.  However, the module will
    # set all the out_nodes to the same logged value, so they'll be
    # identical.
    # TODO: Cache this result.  Some can take a long time,
    # e.g. finding out how paired reads are oriented.
    i = 0
    while i < len(out_data_nodes):
        out_data_node = out_data_nodes[i]
        attr = module.set_out_attributes(antecedents, out_data_node.attributes)
        out_data_node.attributes = attr
    
        # This might generate duplicate out_data_nodes.  Get rid of them.
        j = i+1
        while j < len(out_data_nodes):
            if out_data_node == out_data_nodes[j]:
                del out_data_nodes[j]
            else:
                j += 1
        i += 1

    # Figure out which node in the network is compatible with out_node.
    compatible = []  # list of (out_data_node, next_ids)
    for next_id in network.transitions[module_id]:
        # If this has already been run, then ignore.
        if next_id in pool:
            continue
        next_node = network.nodes[next_id]
        for out_data_node in out_data_nodes:
            # out_attributes should be subset of next_data.attributes.
            if not bie3._is_data_compatible(out_data_node, next_node):
                continue
            compatible.append((out_data_node, next_id))
    # If there are no compatible out nodes, then return None.
    if not compatible:
        return None
    # If there are multiple compatible nodes, arbitrarily use the first one.
    assert len(compatible) >= 1
    out_data_node, next_id = compatible[0]
    
    # Set up the directories and outfile.
    # Unfortunately, can't use timestamp in pathname, or else this
    # will never re-use prior analyses.  Have to be more clever about
    # this.
    h = module.hash_input(
        module_name, antecedents, out_data_node.attributes, user_options)
    ## Get milliseconds.
    #x = time.time()
    #ms = int((x-int(x))*100)
    #ts = time.strftime("%y%m%d.%H%M%S", time.localtime())
    #x = "%s.%02d__%s__B%03d__%s" % (ts, ms, module_name, VERSION, h)
    x = "%s__B%03d__%s" % (module_name, VERSION, h)
    result_dir = os.path.join(output_path, x)
    outfile = module.name_outfile(antecedents, user_options)

    # Create the IdentifiedDataNode that will be the output once the
    # module has been run.
    full_outfile = os.path.join(result_dir, outfile)
    #x = bie3.DataNode(network.nodes[next_id].datatype, **out_data.attributes)
    out_identified_data_node = module_utils.IdentifiedDataNode(
        out_data_node, full_outfile)

    
    #time_str = "[%s]" % time.strftime('%I:%M %p')
    time_str = "[%s]" % time.strftime('%a %I:%M %p')
    
    # Check if this has already been run.
    if _is_module_output_complete(result_dir):
        filename = os.path.join(result_dir, BETSY_PARAMETER_FILE)
        assert os.path.exists(filename)
        params = _read_parameter_file(filename)
        run_time = params["elapsed_pretty"]
        if run_time == "instant":
            x = "ran instantly"
        else:
            x = "took %s" % run_time
        x = "%s  %s (CACHED, previously %s)" % (time_str, module_name, x)
        #parselib.print_split(x, prefixn=2)
        print x
        sys.stdout.flush()
        return out_identified_data_node, next_id

    # Running this module now.
    x = "%s  %s" % (time_str, module_name)
    #parselib.print_split(x, prefixn=2)
    print x
    sys.stdout.flush()

    # Run the module in a temporary directory.
    completed_successfully = False
    cwd = os.getcwd()
    try:
        temp_dir = tempfile.mkdtemp(dir=output_path)
        assert temp_dir != result_dir
        os.chdir(temp_dir)
        temp_outfile = os.path.join(temp_dir, outfile)
        start_time = time.localtime()
        try:
            module.run(
                network, antecedents, out_data_node.attributes, user_options,
                num_cores, temp_outfile)
        except (SystemError, KeyboardInterrupt, MemoryError), x:
            raise
        except Exception, x:
            logging.exception("Exception in module: %s" % module_name)
            raise
        end_time = time.localtime()

        # Write stdout.txt to indicate success.
        assert filelib.fp_exists_nz(temp_outfile), "no file generated"
        success_file = os.path.join(temp_dir, 'stdout.txt')
        open(success_file, 'w').write('success\n')

        # Write parameters.
        x = os.path.join(temp_dir, BETSY_PARAMETER_FILE)
        _write_parameter_file(
            x, pipeline, antecedents, out_data_node.attributes, user_options,
            outfile, start_time, end_time, user, job_name)

        # Move files to the real directory.  If someone else is
        # currently copying files to the exact same directory, then
        # wait for them to finish.
        success_file = os.path.join(result_dir, "stdout.txt")
        copy_files = True
        if os.path.exists(result_dir):
            # Define MODULE_TIMEOUT
            # BUG: This will fail if it takes longer than
            # MODULE_TIMEOUT to copy the files.
            timeout = int(config.MODULE_TIMEOUT)
            # Add a random factor, so two processes that were started
            # at the same time won't collide.
            end_wait = time.time() + timeout + random.randint(0, 59)
            while time.time() < end_wait:
                if os.path.isfile(success_file):
                    # Other process finished. Don't copy.
                    copy_files = False
                    break
                time.sleep(1)
        if copy_files:
            if os.path.exists(result_dir):
                shutil.rmtree(result_dir)
            #shutil.copytree(temp_dir, result_dir)
            shutil.move(temp_dir, result_dir)
        completed_successfully = True
    finally:
        os.chdir(cwd)
        if os.path.exists(temp_dir):
            if clean_up or completed_successfully:
                shutil.rmtree(temp_dir)

    return out_identified_data_node, next_id


def _is_module_output_complete(path):
    # Return a True or False indicating whether path contains the
    # complete results of a previously run module.
    import os

    if not os.path.exists(path):
        return False

    # Make sure "stdout.txt" exists.
    x = os.path.join(path, "stdout.txt")
    if not os.path.exists(x):
        return False
    # Make sure BETSY_PARAMETER_FILE exists.
    x = os.path.join(path, BETSY_PARAMETER_FILE)
    if not os.path.exists(x):
        return False
    # Make sure no broken symlinks.
    for x in os.walk(path, followlinks=True):
        dirpath, dirnames, filenames = x
        for x in dirnames + filenames:
            x = os.path.join(dirpath, x)
            # Returns False for broken links.
            if not os.path.exists(x):
                return False
    return True


def _get_available_input_combinations(
    network, module_id, custom_attributes, pool):
    # Return a list of tuples indicating the node_ids that are:
    # 1.  Valid sets of input data.
    # 2.  Available in the pool.
    # node_ids are in the same order as the module.datatypes.
    import bie3

    # Make a list of every possible combination of inputs that goes
    # into this module.
    x = bie3._bc_to_input_ids(network, module_id, custom_attributes)
    all_input_ids = x

    # See if we have the data to run this module.
    module = network.nodes[module_id]
    available = []
    for input_ids in all_input_ids:
        assert len(module.in_datatypes) == len(input_ids)
        x = [x for x in input_ids if x in pool]
        if len(x) == len(input_ids):
            available.append(x)
    return available


def _pretty_time_delta(delta):
    days, x = divmod(delta, 60*60*24)
    hours, x = divmod(x, 60*60)
    minutes, seconds = divmod(x, 60)

    if not days and not hours and not minutes and seconds < 1:
        return "instant"
    if not days and not hours and not minutes:
        return "%d s" % seconds
    if not days and not hours:
        x = minutes + seconds/60.
        return "%.1f mins" % x
    if not days:
        x = hours + minutes/60. + seconds/3600.
        return "%.1f hrs" % x

    day_or_days = "day"
    if days > 1:
        day_or_days = "days"
    x = hours + minutes/60. + seconds/3600.
    return "%d %s, %.1f hours" % (days, day_or_days, x)


def _write_parameter_file(
    parameter_file, pipeline, antecedents, out_attributes, user_options,
    outfile, start_time, end_time, user, job_name):
    import json
    import time
    import operator

    # Time format.
    FMT = "%a %b %d %H:%M:%S %Y"

    elapsed = time.mktime(end_time) - time.mktime(start_time)

    params = {}
    params["pipeline"] = pipeline
    if not operator.isSequenceType(antecedents):
        antecedents = [antecedents]
    ante = [(x.data.datatype.name, x.identifier) for x in antecedents]
    params["antecedents"] = ante
    params["out_attributes"] = out_attributes
    params["user_options"] = user_options
    params["outfile"] = outfile
    params["start_time"] = time.strftime(FMT, start_time)
    params["end_time"] = time.strftime(FMT, end_time)
    params["elapsed"] = elapsed
    params["elapsed_pretty"] = _pretty_time_delta(elapsed)
    params["user"] = user
    params["job_name"] = job_name
        
    x = json.dumps(params, sort_keys=True, indent=4)
    x = x + "\n"
    open(parameter_file, 'w').write(x)

def _read_parameter_file(filename):
    import json
    return json.loads(open(filename).read())
