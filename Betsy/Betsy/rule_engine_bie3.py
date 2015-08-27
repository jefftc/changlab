"""

Functions:
run_pipeline
run_module

"""
# _get_available_input_combinations
# _write_parameter_file

# TODO: generate figure that show based on the original network all
# the nodes that were run.

VERSION = 4

def run_pipeline(
    network, in_datas, user_attributes, user_options, user=None,
    job_name='', num_cores=8):
    # Run the pipeline that is indicated by the network.  Returns a
    # filename if successful or None if not.
    #
    # in_datas         List of IdentifiedDataNodes.
    # user_attributes  From --dattr.  AttributeDef
    # user_options     From --mattr.  OptionDef
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

    # Each member of the stack can be a tuple of:
    # 1.  (IdentifiedDataNode, node_id, None)
    # 2.  (ModuleNode, node_id, None)
    # 3.  (ModuleNode, node_id, antecedent_ids)
    stack = []
    for node in in_datas:
        data_node = node.data
        start_node_ids = bie3._find_start_nodes(network, data_node)
        assert start_node_ids, \
               '%s cannot matched to network' % data_node.datatype.name
        for start_node_id in start_node_ids:
            stack.append((node, start_node_id, None))
            #pool[start_node_id] = idata
            
    # Keep track of nodes that have already been generated.
    pool = {}      # dict of node_id -> IdentifiedDataNode
    pipeline = []  # list of the names of the modules run.
    next_node = None
    # Keep track of the number of times a module can't be run because
    # it doesn't have enough input data.  If this exceeds the total
    # number of modules, then quit.
    num_failures = 0
    success = False
    while stack:
        assert num_failures < len(stack), \
               "Inference error: No more nodes to run."
        node, node_id, more_info = stack.pop()
        if node_id == 0:
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
                stack.append((next_node, next_id, None))
            # If I've added new modules, then reset the failures counter.
            num_failures = 0
        elif isinstance(node, bie3.ModuleNode) and more_info is None:
            module_id = node_id

            # If the input data for this module doesn't exist, then
            # just try it again later.
            all_antecedent_ids = _get_available_input_combinations(
                network, module_id, user_attributes, pool)
            if not all_antecedent_ids:
                # Put back to the bottom of the stack.
                stack.insert(0, (node, node_id, more_info))
                num_failures += 1
                continue
            for antecedent_ids in all_antecedent_ids:
                x = node, node_id, antecedent_ids
                stack.append(x)
            # If I've added new analyses to run, then reset the
            # failures counter.
            num_failures = 0
        elif isinstance(node, bie3.ModuleNode):
            # Run this module.
            pipeline.append(node.name)
            x = run_module(
                network, pipeline, module_id, antecedent_ids,
                user_attributes, user_options, 
                pool, user, job_name, num_cores=num_cores)
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
            stack.append((next_node, next_id, None))
            
            #flag = False
            #if not out_nodes:
            #    break
            #for x in out_nodes:
            #    next_node, next_id = x
            #    if next_id == 0:
            #        flag = True
            #        break
            #    elif next_id is None:
            #        raise ValueError('cannot match the output node')
            #    stack.append((next_node, next_id))
            #if flag:
            #    break
        else:
            raise AssertionError

    if success:
        msg = "Completed"
        print "[%s]  %s" % (time.strftime('%I:%M %p'), msg)
        #print next_node.identifier
        #sys.stdout.flush()
        assert next_node
        return next_node.identifier
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
    network, pipeline, module_id, input_ids, user_attributes,
    all_user_options, pool, user, job_name='', clean_up=True, num_cores=8):
    # Return tuple of (IdentifiedDataNode, node_id) for the node that
    # was created.  Returns None if this module fails (no compatible
    # output nodes, or all output nodes already generated).

    # TODO: Remove user_attributes.  Used in the generation of the
    # network.  Not necessary anymore.
    import os
    import sys
    import time
    import tempfile
    import shutil
    import logging
    
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
    x = __import__(
        'modules.'+module_name, globals(), locals(), [module_name], -1)
    module = x.Module()

    # Running this module now.
    print "[%s]  %s" % (time.strftime('%I:%M %p'), module_name)
    sys.stdout.flush()

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
    out_data_nodes = bie3._forwardchain_to_outputs(module_node, x)
    # Make sure the DataTypes are the same.
    assert out_data_nodes
    x = sorted([x.datatype for x in out_data_nodes])
    assert x[0] == x[-1], "ModuleNode points to different DataTypes."
    out_data_node = out_data_nodes[0]

    # Optimization: set_out_attributes can be computationally
    # expensive, so assume that the module will set the output
    # attributes for all the out_nodes identically.  For example, the
    # out_nodes may differ by logged value.  However, the module will
    # set all the out_nodes to the same logged value, so they'll be
    # identical.
    attr = module.set_out_attributes(antecedents, out_data_node.attributes)
    out_data_node.attributes = attr

    # Figure out which node in the network is compatible with out_node.
    compatible = []  # list of next_ids
    for next_id in network.transitions[module_id]:
        # If this has already been run, then ignore.
        if next_id in pool:
            continue
        next_node = network.nodes[next_id]
        # out_attributes should be subset of next_data.attributes.
        if not bie3._is_data_compatible_with(out_data_node, next_node):
            continue
        compatible.append(next_id)
    # If there are no compatible out nodes, then return None.
    if not compatible:
        return None
    assert len(compatible) == 1
    next_id = compatible[0]
    
    # Set up the directories and outfile.
    # Unfortunately, can't use timestamp in pathname, or else this
    # will never re-use prior analyses.  Have to be more clever about
    # this.
    h = module.hash_input(
        pipeline, antecedents, out_data_node.attributes, user_options)
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

    # Check if this has already been run.
    if os.path.exists(os.path.join(result_dir, 'stdout.txt')):
        return out_identified_data_node, next_id

    # Run the module in a temporary directory.
    cwd = os.getcwd()
    try:
        temp_dir = tempfile.mkdtemp()
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
        assert module_utils.exists_nz(temp_outfile), "no file generated"
        success_file = os.path.join(temp_dir, 'stdout.txt')
        open(success_file, 'w').write('success')

        # Write parameters.
        x = os.path.join(temp_dir, 'BETSY_parameters.txt')
        _write_parameter_file(
            x, pipeline, antecedents, out_data_node.attributes, user_options,
            outfile, start_time, end_time, user, job_name)

        # Move files to the real directory.  If someone else is
        # currently copying files to the exact same directory, then
        # wait for them to finish.
        copy_files = True
        if os.path.exists(result_dir):
            timeout = int(config.MODULE_TIMEOUT)
            end_wait = time.localtime() + timeout
            while time.localtime() < end_wait:
                if os.path.isfile(success_file):
                    # Other process finished. Don't copy.
                    copy_files = False
                    break
                time.sleep(1)
        if copy_files:
            if os.path.exists(result_dir):
                shutil.rmtree(result_dir)
            shutil.copytree(temp_dir, result_dir)
    finally:
        os.chdir(cwd)
        if clean_up and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

    return out_identified_data_node, next_id


def _get_available_input_combinations(
    network, module_id, user_attributes, pool):
    # Return a list of tuples indicating the node_ids that are:
    # 1.  Valid sets of input data.
    # 2.  Available in the pool.
    # node_ids are in the same order as the module.datatypes.
    import bie3

    # Make a list of every possible combination of inputs that goes
    # into this module.
    prev_ids = network.points_to(module_id)
    all_input_ids = bie3._get_valid_input_combinations(
        network, module_id, prev_ids, user_attributes)

    # See if we have the data to run this module.
    available = []
    for input_ids in all_input_ids:
        x = [x for x in input_ids if x in pool]
        if len(x) == len(input_ids):
            available.append(x)
    return available


## def get_out_node(
##     working_dir, network, module_id, module_file, parameters, pool,
##     user_options, user_attributes):
##     # Returns a list of tuples (out_node, out_id, out_identifier).
##     import os
##     from Betsy import bie3
##     from Betsy import rulebase
##     from Betsy import module_utils

##     # Need to change directory because find_antecedents checks for the
##     # existence of the file.  Assumes the files are in the current
##     # directory.  Should fix this by making the identifiers for
##     # DataObject the full path name.
##     current_dir = os.getcwd()
##     try:
##         os.chdir(working_dir)
##         data_node = module_file.find_antecedents(
##             network, module_id, parameters, user_attributes, pool)
##         outfile = module_file.name_outfile(data_node, user_options)
##     finally:
##         os.chdir(current_dir)
##     next_possible_ids = network.transitions[module_id]
##     fn = getattr(
##         rulebase, network.nodes[next_possible_ids[0]].datatype.name)
##     out_node = bie3.DataNode(fn, **parameters)
##     out_object = module_utils.IdentifiedDataNode(
##         out_node, os.path.join(working_dir, outfile))
##     out_id = None
##     result = []
##     for next_id in next_possible_ids:
##         if compare_two_dict(parameters, network.nodes[next_id].attributes):
##             out_id = next_id
##             result.append((out_object, out_id))
##             #result.append((out_node, out_id, outfile))
##     return result


## def compare_two_dict(dict_A, dict_B):
##     if len(dict_A) != len(dict_B):
##         return False
##     if set(dict_A) - set(dict_B):
##         return False
##     for key in dict_A:
##         setA = dict_A[key]
##         setB = dict_B[key]
##         if isinstance(setA, str):
##             setA = [setA]
##         if isinstance(setB, str):
##             setB = [setB]
##         if (not set(setA).issubset(set(setB)) and
##             not set(setB).issubset(set(setA))):
##             return False
##     return True


## def copy_result_folder(working_dir, temp_dir, temp_outfile):
##     """copy result files in temp folder to result folder,
##       if fails, delete result folder"""
##     import os
##     import shutil
    
##     if not os.path.exists(working_dir):
##         os.mkdir(working_dir)
##     try:
##         shutil.copyfile(os.path.join(temp_dir, 'Betsy_parameters.txt'),
##                         os.path.join(working_dir, 'Betsy_parameters.txt'))
##         if os.path.isdir(os.path.join(temp_dir, temp_outfile)):
##             shutil.copytree(os.path.join(temp_dir, temp_outfile),
##                             os.path.join(working_dir, temp_outfile))
##         else:
##             #change to copy other files
##             all_files = os.listdir(temp_dir)
##             extra_files = list(set(all_files) - set(
##                 [temp_outfile, 'Betsy_parameters.txt', 'stdout.txt']))
##             shutil.copyfile(os.path.join(temp_dir, temp_outfile),
##                             os.path.join(working_dir, temp_outfile))
##             for extra_file in extra_files:
##                 if os.path.isfile(os.path.join(temp_dir, extra_file)):
##                     shutil.copyfile(os.path.join(temp_dir, extra_file),
##                                     os.path.join(working_dir, extra_file))
##                 elif os.path.isdir(os.path.join(temp_dir, extra_file)):
##                     shutil.copytree(os.path.join(temp_dir, extra_file),
##                                     os.path.join(working_dir, extra_file))
##         shutil.copyfile(os.path.join(temp_dir, 'stdout.txt'),
##                         os.path.join(working_dir, 'stdout.txt'))
##     finally:
##         if not os.path.isfile(os.path.join(working_dir, 'stdout.txt')):
##             shutil.rmtree(working_dir)
##             raise ValueError('copy fails')



def _write_parameter_file(
    parameter_file, pipeline, antecedents, out_attributes, user_options,
    outfile, start_time, end_time, user, job_name):
    import json
    import time
    import operator

    FMT = "%a %b %d %H:%M:%S %Y"

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
    params["user"] = user
    params["job_name"] = job_name
        
    x = json.dumps(params, sort_keys=True, indent=4)
    open(parameter_file, 'w').write(x)
