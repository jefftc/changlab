"""

Functions:
run_pipeline
run_module

have_input_data

make_module_wd_name
get_out_node
create_out_attributes

compare_two_dict
copy_result_folder

TODO: does module.find_antecedents need out_attributes?
TODO: Rename DataObject.  Figure out better way to handle this.
TODO: when running pipeline, should show all possible options

TODO: run_pipeline should not print to screen.



user_input  something from the user that is passed to the module.  OptionDef
Should be renamed user_options

user_attributes  from --dattr  AttributeDef(?)
user_options     from --mattr  OptionDef(?)


Functions in each module:

run(network, antecedents, out_attributes, user_options, num_cores)
  network         Network object
  antecedents     Either a DataObject, or a sequence of DataObjects.
  out_attributes  The desired attributes of the output object.
  user_options    Options for the module provided by the user.
  num_cores
  Returns a DataObject.

find_antecedents(network, module_id, out_attributes, user_attributes, pool)
  Gets from the pool.
  Returns a Data Object, or a sequence of DataObjects.
TODO: check module_utils.get_identifier, when you need multiple antecedents

set_out_attributes(antecedents, out_attributes)
Takes the default attributes (as a dictionary) of the output node, and
sets it based on the input data.  This is for modules that can
determine the output based on the input data, e.g. a module that
checks if a user's input file is logged.

name_outfile(antecedents, user_options)
Come up with a human readable name for the output file.
TODO: Should return just file, rather than the full path.
  
make_unique_hash(pipeline, antecedents, out_attributes, user_options)
TODO: Rename to hash_module???
TODO: rename data_node to antecedents?

TODO: why can't get_antecedents be done automatically based on the rules?


"""

def make_module_wd_name(
    network, module_file, module_id, user_options, pipeline, pool,
    out_attributes, user_attributes):
    # what's the difference between current_attributes and user_attributes?
    
    antecedents = module_file.find_antecedents(
        network, module_id, out_attributes, user_attributes, pool)
    uhash = module_file.make_unique_hash(
        pipeline, antecedents, out_attributes, user_options)
    
    module_name = network.nodes[module_id].name
    working_dir = module_name + '_BETSYHASH1_' + uhash
    return working_dir


## def make_module_wd(working_dir):
##     import os
    
##     os.mkdir(working_dir)


def copy_result_folder(working_dir, temp_dir, temp_outfile):
    """copy result files in temp folder to result folder,
      if fails, delete result folder"""
    import os
    import shutil
    
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    try:
        shutil.copyfile(os.path.join(temp_dir, 'Betsy_parameters.txt'),
                        os.path.join(working_dir, 'Betsy_parameters.txt'))
        if os.path.isdir(os.path.join(temp_dir, temp_outfile)):
            shutil.copytree(os.path.join(temp_dir, temp_outfile),
                            os.path.join(working_dir, temp_outfile))
        else:
            #change to copy other files
            all_files = os.listdir(temp_dir)
            extra_files = list(set(all_files) - set(
                [temp_outfile, 'Betsy_parameters.txt', 'stdout.txt']))
            shutil.copyfile(os.path.join(temp_dir, temp_outfile),
                            os.path.join(working_dir, temp_outfile))
            for extra_file in extra_files:
                if os.path.isfile(os.path.join(temp_dir, extra_file)):
                    shutil.copyfile(os.path.join(temp_dir, extra_file),
                                    os.path.join(working_dir, extra_file))
                elif os.path.isdir(os.path.join(temp_dir, extra_file)):
                    shutil.copytree(os.path.join(temp_dir, extra_file),
                                    os.path.join(working_dir, extra_file))
        shutil.copyfile(os.path.join(temp_dir, 'stdout.txt'),
                        os.path.join(working_dir, 'stdout.txt'))
    finally:
        if not os.path.isfile(os.path.join(working_dir, 'stdout.txt')):
            shutil.rmtree(working_dir)
            raise ValueError('copy fails')


def compare_two_dict(dict_A, dict_B):
    if len(dict_A) != len(dict_B):
        return False
    if set(dict_A) - set(dict_B):
        return False
    for key in dict_A:
        setA = dict_A[key]
        setB = dict_B[key]
        if isinstance(setA, str):
            setA = [setA]
        if isinstance(setB, str):
            setB = [setB]
        if (not set(setA).issubset(set(setB)) and
            not set(setB).issubset(set(setA))):
            return False
    return True


# What does this do?
def get_out_node(
    working_dir, network, module_id, module_file, parameters, pool,
    user_input, user_attributes):
    # Returns a list of tuples (out_node, out_id, out_identifier).
    import os

    import bie3
    import rulebase
    from module_utils import DataObject    

    # Need to change directory because find_antecedents checks for the
    # existence of the file.  Assumes the files are in the current
    # directory.  Should fix this by making the identifiers for
    # DataObject the full path name.
    # TODO: make identifiers for DataObject the full path name.
    # TODO: don't change directory
    current_dir = os.getcwd()
    try:
        os.chdir(working_dir)
        data_node = module_file.find_antecedents(
            network, module_id, parameters, user_attributes, pool)
        outfile = module_file.name_outfile(data_node, user_input)
    finally:
        os.chdir(current_dir)
    next_possible_ids = network.transitions[module_id]
    fn = getattr(
        rulebase, network.nodes[next_possible_ids[0]].datatype.name)
    out_node = bie3.Data(fn, **parameters)
    out_object = DataObject(out_node, outfile)
    out_id = None
    result = []
    for next_id in next_possible_ids:
        if compare_two_dict(parameters, network.nodes[next_id].attributes):
            out_id = next_id
            result.append((out_object, out_id))
            #result.append((out_node, out_id, outfile))
    return result


def create_out_attributes(
    network, module_id, module_file, pool, user_attributes):
    # return out_attributes, if the output already generated, return None
    # Does this for the next Data node that hasn't been calculated.
    next_ids = network.transitions[module_id]
    for next_id in next_ids:
        if next_id in pool:
            continue
        out_attributes = network.nodes[next_id].attributes
        data_node = module_file.find_antecedents(
            network, module_id, out_attributes, user_attributes, pool)
        parameters = module_file.set_out_attributes(data_node, out_attributes)
        return parameters
    return None


## def choose_next_module(network, node_id, pool, user_attributes):
##     # choose the next module and return a list of (module, id)
##     import bie3

##     assert node_id in network.transitions
##     #if node_id not in network.transitions:
##     #    return False
##     result = []
##     for next_id in network.transitions[node_id]:
##         next_node = network.nodes[next_id]
##         assert isinstance(next_node, bie3.Module), \
##                'next node supposed to be a module'
##         # Why does this have to be here?
##         #if have_input_data(network, next_node_id, pool, user_attributes):
##         #    result.append((next_node, next_node_id))
##         result.append((next_node, next_id))
##     #result.sort(key=lambda x: x[1], reverse=False)
##     return result


def have_input_data(network, module_id, user_attributes, pool):
    # Return a boolean if the data exists for the module to be run.
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


def run_module(
    network, pipeline_sequence, module_id, user_attributes, user_options, pool,
    user, job_name='', clean_up=True, num_cores=8):
    # return out_nodes, if module fails, the module itself will raise a break,
    # if cannot find a outnode inside the network, will return a empty list
    import os
    import sys
    import getpass
    import time
    import tempfile
    import shutil
    import logging
    
    import config
    import module_utils

    assert user
    output_path = config.OUTPUTPATH
    assert os.path.exists(output_path), \
           'the output_path %s does not exist' % output_path
    assert network.transitions[module_id]

    module_node = network.nodes[module_id]
    module_name = module_node.name
    module = __import__(
        'modules.'+module_name, globals(), locals(), [module_name], -1)

    # Pull out the user inputs that are used for this module.
    sub_user_input = {}
    for option in module_node.option_defs:
        if option.default is None:
            assert option.name in user_options, \
                   'Missing input: %s' % option.name
        if option.name in user_options:
            sub_user_input[option.name] = user_options[option.name]
            
    out_attributes = create_out_attributes(
        network, module_id, module, pool, user_attributes)
    if out_attributes is None:
        # When does this happen?  When this module doesn't point to
        # anything.
        return []
    
    print "[%s]  %s" % (time.strftime('%I:%M %p'), module_name)
    #print '[' + time.strftime('%l:%M%p') + '].' + module_name
    sys.stdout.flush()

    x = make_module_wd_name(
        network, module, module_id, sub_user_input, pipeline_sequence, pool,
        out_attributes, user_attributes)
    working_dir = os.path.join(output_path, x)
    # make name of outfile
    data_node = module.find_antecedents(
        network, module_id, out_attributes, user_attributes, pool)
    outfile = os.path.split(module.name_outfile(data_node, sub_user_input))[-1]
    temp_dir = ''
    #temp_outfile = ''

    # If this has already been generated, then
    # XXX
    
    if not os.path.exists(os.path.join(working_dir, 'stdout.txt')):
        #if no result has been generated, create temp folder and run analysis
        cwd = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        try:
            os.chdir(temp_dir)
            starttime = time.strftime(module_utils.FMT, time.localtime())
            #out_node = module.run(
            #    data_node, out_attributes, sub_user_input, network, num_cores)
            out_node = module.run(
                network, data_node, out_attributes, sub_user_input, num_cores)
            module_utils.write_Betsy_parameters_file(
                out_node.data.attributes, data_node, outfile, sub_user_input,
                pipeline_sequence, starttime, user, job_name)
            f = file(os.path.join(temp_dir, 'stdout.txt'), 'w')
            f.write('This module runs successully.')
            f.close()
        except:
            # BUG: Should try to figure out what kind of exception this is.
            # Also, should let MemoryErrors and System errors go.
            logging.exception('Got exception on %s .run()' % module_name)
            raise
        finally:
            os.chdir(cwd)
        try:
            # found if the same analysis has been run and wait for it finished
            try:
                os.mkdir(working_dir)
            except OSError:
                try:
                    module_timeout = int(config.MODULE_TIMEOUT)
                except AttributeError:
                    module_timeout = 60
                i = 0
                while i <= module_timeout:
                    if os.path.isfile(os.path.join(working_dir, 'stdout.txt')):
                        break
                    i = i + 1
                    time.sleep(1)
            # if previous analysis not working, copy the current
            # results to working_dir
            if (module_utils.exists_nz(os.path.join(temp_dir, outfile)) and
                not os.path.exists(os.path.join(working_dir, 'stdout.txt'))):
                copy_result_folder(working_dir, temp_dir, outfile)
        finally:
            if clean_up and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    out_nodes = get_out_node(
        working_dir, network, module_id, module, out_attributes, pool,
        sub_user_input, user_attributes)
    assert out_nodes, 'module %s fails' % module_node.name
    return out_nodes


def run_pipeline(
    network, in_objects, user_attributes, user_inputs, user=None,
    job_name='', num_cores=8):
    # Returns a filename or None.
    import os
    import sys
    import getpass
    import logging
    import time

    import bie3
    import config
    import module_utils

    user = user or getpass.getuser()
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    # Is this thread-safe?
    LOG_FILENAME = os.path.join(output_path, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)

    # What is the pool for?
    stack = []  # list of (DataObject or Module, node_id) to be run.
    for data_object in in_objects:
        data_node = data_object.data
        start_node_ids = bie3._find_start_nodes(network, data_node)
        assert start_node_ids, \
               '%s cannot matched to network' % data_node.datatype.name
        for start_node_id in start_node_ids:
            stack.append((data_object, start_node_id))
            #pool[start_node_id] = idata
            
    # Keep track of nodes that have already been generated.
    pool = {}      # dict of node_id -> (Data, identifier) # XXX implement!
    pipeline = []  # list of ???
    next_node = None
    num_failures = 0
    flag = False
    while stack:
        assert num_failures < len(stack)
        data_object, node_id = stack.pop()
        if isinstance(data_object, module_utils.DataObject):
            pool[node_id] = data_object
            if node_id == 0:
                # When does this occur?  Failed?
                return None
            # Add the next modules into the stack.
            for next_id in network.transitions[node_id]:
                next_node = network.nodes[next_id]
                assert isinstance(next_node, bie3.Module)
                stack.append((next_node, next_id))
            #x = choose_next_module(network, node_id, pool, user_attributes)
            #stack.extend(x)
        elif isinstance(data_object, bie3.Module):
            # Should rename data_object data_or_module.
            module_id = node_id

            # If the input data for this module doesn't exist, then
            # just try it again later.
            if not have_input_data(network, module_id, user_attributes, pool):
                # Put back to the back of the stack.
                stack.insert(0, (data_object, node_id))
                num_failures += 1
                continue

            # Run this module.
            pipeline.append(data_object.name)
            #print module_id
            out_nodes = run_module(
                network, pipeline, module_id, user_attributes, user_inputs, 
                pool, user, job_name, num_cores=num_cores)
            
            flag = False
            if not out_nodes:
                break
            for x in out_nodes:
                next_node, next_id = x
                if next_id == 0:
                    flag = True
                    break
                elif next_id is None:
                    raise ValueError('cannot match the output node')
                stack.append((next_node, next_id))
            if flag:
                break
            
    if flag and next_node and module_utils.exists_nz(next_node.identifier):
        msg = "Completed successfully and generated a file:"
        print "[%s]  %s" % (time.strftime('%I:%M %p'), msg)
        print next_node.identifier
        sys.stdout.flush()
        return next_node.identifier
    else:
        print 'This pipeline has completed unsuccessfully'
        raise ValueError('there is no output for this pipeline')
    return None
