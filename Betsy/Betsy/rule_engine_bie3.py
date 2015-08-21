"""

Functions:
run_pipeline
run_module

has_input_data
write_parameter_file


user_attributes  from --dattr  AttributeDef(?)
user_options     from --mattr  OptionDef(?)

TODO: when running pipeline, should show all possible options
TODO: run_pipeline should not print to screen.

"""

VERSION = 4

def run_pipeline(
    network, in_datas, user_attributes, user_options, user=None,
    job_name='', num_cores=8):
    # Run the pipeline that is indicated by the network.  Returns a
    # filename if successful or None if not.
    #
    # in_datas         List of IdentifiedDataNodes.
    # user_attributes  from --dattr  AttributeDef(?)
    # user_options     from --mattr  OptionDef(?)
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

    # What is the pool for?
    stack = []  # list of (IdentifiedDataNode or ModuleNode, node_id) to be run
    for node in in_datas:
        data_node = node.data
        start_node_ids = bie3._find_start_nodes(network, data_node)
        assert start_node_ids, \
               '%s cannot matched to network' % data_node.datatype.name
        for start_node_id in start_node_ids:
            stack.append((node, start_node_id))
            #pool[start_node_id] = idata
            
    # Keep track of nodes that have already been generated.
    pool = {}      # dict of node_id -> IdentifiedDataNode
    pipeline = []  # list of ???
    next_node = None
    num_failures = 0
    #flag = False
    success = False
    while stack:
        assert num_failures < len(stack)
        node, node_id = stack.pop()
        if isinstance(node, module_utils.IdentifiedDataNode):
            pool[node_id] = node
            assert node_id != 0
            #if node_id == 0:
            #    # When does this occur?  Failed?
            #    return None
            # Add the next modules into the stack.
            for next_id in network.transitions[node_id]:
                next_node = network.nodes[next_id]
                assert isinstance(next_node, bie3.ModuleNode)
                stack.append((next_node, next_id))
        elif isinstance(node, bie3.ModuleNode):
            module_id = node_id

            # If the input data for this module doesn't exist, then
            # just try it again later.
            if not has_input_data(network, module_id, user_attributes, pool):
                # Put back to the back of the stack.
                stack.insert(0, (node, node_id))
                num_failures += 1
                continue

            # Run this module.
            pipeline.append(node.name)
            x = run_module(
                network, pipeline, module_id, user_attributes, user_options, 
                pool, user, job_name, num_cores=num_cores)
            if x is None:
                continue
            next_node, next_id = x
            if next_id == 0:
                success = True
                break
            assert next_id is not None
            stack.append((next_node, next_id))
            
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

    if success:
        msg = "Completed"
        print "[%s]  %s" % (time.strftime('%I:%M %p'), msg)
        #print next_node.identifier
        #sys.stdout.flush()
        return next_node.identifier
    print 'This pipeline has completed unsuccessfully'
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
    network, pipeline, module_id, user_attributes, all_user_options, pool,
    user, job_name='', clean_up=True, num_cores=8):
    # return tuple of (node, node_id) or None, if the
    # module failed.
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

    # Get the user_options.  all_user_options contains all the ones
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
    # logged or not).  Figure out which one fits this module.
    next_ids = network.transitions[module_id]
    assert next_ids
    # Make sure the DataTypes are the same.
    x = sorted([network.nodes[x].datatype for x in next_ids])
    assert x[0] == x[-1], "ModuleNode points to different DataTypes."
    # Find a compatible one.
    compatible = []  # list of (next_id, out_attributes)
    for next_id in network.transitions[module_id]:
        # Get the out_attributes and antecedents.
        attr1 = network.nodes[next_id].attributes
        antecedents = module.find_antecedents(
            network, module_id, attr1, user_attributes, pool)
        # Let the module set it based on the antecedents.  E.g. for
        # modules like is_logged.
        attr2 = module.set_out_attributes(antecedents, attr1.copy())
        if attr1 != attr2:
            continue
        x = next_id, attr2
        compatible.append(x)
    if not compatible:
        return None
    assert len(compatible) == 1
    next_id, out_attributes = compatible[0]

    # Set up the directories and outfile.
    # Unfortunately, can't use timestamp in pathname, or else this
    # will never re-use prior analyses.  Have to be more clever about
    # this.
    ## Get milliseconds.
    h = module.hash_input(pipeline, antecedents, out_attributes, user_options)
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
    x = bie3.DataNode(network.nodes[next_id].datatype, **out_attributes)
    out_node = module_utils.IdentifiedDataNode(x, full_outfile)

    # Check if this has already been run.
    if os.path.exists(os.path.join(result_dir, 'stdout.txt')):
        return out_node, next_id

    # Run the module in a temporary directory.
    cwd = os.getcwd()
    try:
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        temp_outfile = os.path.join(temp_dir, outfile)
        start_time = time.localtime()
        try:
            module.run(
                network, antecedents, out_attributes, user_options, num_cores,
                temp_outfile)
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
        write_parameter_file(
            x, pipeline, antecedents, out_attributes, user_options,
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

    return out_node, next_id


def has_input_data(network, module_id, user_attributes, pool):
    # Return a boolean indicating whether this module has all the
    # input data needed.
    import bie3

    # Make a list of every possible combination of inputs that goes
    # into this module.
    prev_ids = network.points_to(module_id)
    all_input_ids = bie3._get_valid_input_combinations(
        network, module_id, prev_ids, user_attributes)

    # See if we have the data to run this module.
    for input_ids in all_input_ids:
        x = [x for x in input_ids if x in pool]
        if len(x) == len(input_ids):
            return True
    return False


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
##     # TODO: make identifiers for DataObject the full path name.
##     # TODO: don't change directory
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



def write_parameter_file(
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
