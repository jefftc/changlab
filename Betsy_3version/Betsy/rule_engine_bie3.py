#rule_engine_bie3.py
import rulebase
import config
import bie3
import logging
import sys
import tempfile
import shutil
import module_utils
from  module_utils import DataObject as DataObject
import os
import time
import getpass
from time import strftime, localtime


"""
Functions:
make_module_wd_name
make_module_wd
copy_result_folder
compare_two_dict
get_out_node
create_out_attributes
choose_next_module
test_require_data
run_module
run_pipeline
"""


def make_module_wd_name(network, module_file, module_id,
                        user_input, pipeline_sequence,
                        pool, current_attributes):
    module_name = network.nodes[module_id].name
    data_object = module_file.find_antecedents(network, module_id,
                                               pool, current_attributes)
    hash_string = module_file.make_unique_hash(
        data_object, pipeline_sequence, current_attributes, user_input)
    working_dir = module_name + '_BETSYHASH1_' + hash_string
    return working_dir


def make_module_wd(working_dir):
    os.mkdir(working_dir)


def copy_result_folder(working_dir, temp_dir, temp_outfile):
    """copy result files in temp folder to result folder,
      if fails, delete result folder"""
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    try:
        shutil.copyfile(os.path.join(temp_dir, 'Betsy_parameters.txt'),
                        os.path.join(working_dir, 'Betsy_parameters.txt'))
        if os.path.isdir(temp_outfile):
            shutil.copytree(os.path.join(temp_dir, temp_outfile),
                            os.path.join(working_dir, temp_outfile))
        else:
            #change to copy other files
            all_files = os.listdir(temp_dir)
            extra_files = list(set(all_files) - set([temp_outfile,
                                                     'Betsy_parameters.txt',
                                                     'stdout.txt']))
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


def get_out_node(network, module_id, module_file, parameters,
                pool, user_input):
    data_node = module_file.find_antecedents(network,
                                             module_id, pool, parameters)
    outfile = module_file.name_outfile(data_node, user_input)
    next_possible_ids = network.transitions[module_id]
    out_object = None
    fn = getattr(rulebase, network.nodes[next_possible_ids[0]].datatype.name)
    out_node = bie3.Data(fn, **parameters)
    out_object = DataObject(out_node, outfile)
    out_id = None
    result = []
    for next_id in next_possible_ids:
        if compare_two_dict(parameters, network.nodes[next_id].attributes):
            out_id = next_id
            result.append((out_object, out_id))
    return result


def create_out_attributes(network, module_id, module_file, pool):
    # return out_attributes, if the output already generated, return None
    next_ids = network.transitions[module_id]
    for next_id in next_ids:
        if next_id in pool:
            continue
        out_attributes = network.nodes[next_id].attributes
        data_node = module_file.find_antecedents(
            network, module_id, pool, out_attributes)
        parameters = module_file.get_out_attributes(
            out_attributes, data_node)
        return parameters
    return None


def choose_next_module(network, node_id, pool):
    # choose the next module and return a list of (module,id)
    if node_id not in network.transitions:
        return False
    next_node_ids = network.transitions[node_id]
    result = []
    for next_node_id in next_node_ids:
        next_node = network.nodes[next_node_id]
        assert isinstance(next_node, bie3.Module),'next node supposed to be a module'
        if test_require_data(network, next_node_id, pool):
            result.append((next_node, next_node_id))
    result.sort(key=lambda x: x[1], reverse=False)
    return result


def test_require_data(network, module_id, pool):
    # test if the required data for the module is met.
    require_id = []
    for key in network.transitions:
        if module_id in network.transitions[key]:
            require_id.append(key)
    combine_ids = bie3._get_valid_input_combinations(
        network, module_id, require_id)
    for combine_id in combine_ids:
        flag = True
        for i in combine_id:
            if i not in pool:
                flag = False
        if flag:
            return True
    return False


def run_module(network, module_id, pool, user_inputs, pipeline_sequence,
               user=getpass.getuser(), job_name='', clean_up=True):
    output_path = config.OUTPUTPATH
    assert os.path.exists(output_path), (
        'the output_path %s does not exist' % output_path)
    # get module
    module_node = network.nodes[module_id]
    sub_user_input = {}
    if module_node.user_inputs:
        for user_in in module_node.user_inputs:
            if user_in.name in user_inputs:
                sub_user_input[user_in.name] = user_inputs[user_in.name]
            if user_in.default is None:
                assert user_in.name in user_inputs, (
                    '%s should be specified' % user_in.name)
    module_name = module_node.name
    module = __import__('modules.' + module_name, globals(),
                        locals(), [module_name], -1)
    out_attributes = create_out_attributes(
        network, module_id, module, pool)
    if not out_attributes:
        return []
    print '[' + time.strftime('%l:%M%p') + '].' + module_name
    working_dir = os.path.join(output_path, make_module_wd_name(
        network, module, module_id, sub_user_input, pipeline_sequence,
        pool, out_attributes))
    # make name of outfile
    data_node = module.find_antecedents(network, module_id,
                                        pool, out_attributes)
    outfile = module.name_outfile(data_node, sub_user_input)
    temp_dir = ''
    temp_outfile = ''
    if not os.path.exists(os.path.join(working_dir, 'stdout.txt')):
        #if no result has been generated, create temp folder and run analysis
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        try:
            starttime = strftime(module_utils.FMT, localtime())
            out_node = module.run(data_node, out_attributes,
                                  sub_user_input, network)
            module_utils.write_Betsy_parameters_file(out_node.data.attributes,
                                                     data_node, sub_user_input,
                                                     pipeline_sequence,
                                                     starttime, user, job_name)
            if out_node:
                f = file(os.path.join(temp_dir, 'stdout.txt'), 'w')
                f.write('This module runs successully.')
                f.close()
            else:
                return False
        except:
            logging.exception('Got exception on %s .run()' % module_name)
            raise
        temp_outfile = os.path.split(outfile)[-1]
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
            if (module_utils.exists_nz(temp_outfile) and not
                os.path.exists(os.path.join(working_dir, 'stdout.txt'))):
                copy_result_folder(working_dir, temp_dir, temp_outfile)
        finally:
            if clean_up and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    os.chdir(working_dir)
    out_nodes = get_out_node(
        network, module_id, module, out_attributes, pool, sub_user_input)
    assert out_nodes, 'module %s fails' % module_node.name
    return out_nodes


def run_pipeline(network, in_objects, user_inputs,
                 user=getpass.getuser(), job_name=''):
    working_dir = os.getcwd()
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    LOG_FILENAME = os.path.join(output_path, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    try:
        pipeline_sequence = []
        stack_list = []
        pool = {}
        for data_object in in_objects:
            data_node = data_object.data
            start_node_ids = bie3._find_start_nodes(network, data_node)
            for start_node_id in start_node_ids:
                stack_list.append((data_object, start_node_id))
                pool[start_node_id] = data_object
        num_failures = 0
        while stack_list:
            assert num_failures < len(stack_list)
            data_object, node_id = stack_list.pop()
            data_node = data_object
            if isinstance(data_node, DataObject):
                data_node = data_node.data
            if isinstance(data_node, bie3.Data):
                pool[node_id] = data_object
                result = choose_next_module(network, node_id, pool)
                if result:
                    result.sort(key=lambda x: x[1])
                    for x in result:
                        module, module_id = x
                        stack_list.append((module, module_id))
            elif isinstance(data_node, bie3.Module):
                module_id = node_id
                test_required = test_require_data(network, module_id, pool)
                if test_required:
                    pipeline_sequence.append(module.name)
                    out_nodes = run_module(network, module_id, pool,
                                           user_inputs, pipeline_sequence,
                                           user, job_name)
                    for x in out_nodes:
                        next_node, next_id = x
                        if next_id == 0:
                            if module_utils.exists_nz(next_node.identifier):
                                print ('[' + time.strftime('%l:%M%p') +
                                       '] Completed successfully and ' +
                                       'generated a file:')
                                print  next_node.identifier + '\r'
                                print '\r'
                                sys.stdout.flush()
                                return next_node.identifier
                            else:
                                print 'This pipeline has\
                                       completed unsuccessfully'
                                raise ValueError(
                                    'there is no output for this pipeline')
                            return None
                        elif next_id is None:
                            raise ValueError('cannot match the output node')
                        stack_list.append((next_node, next_id))
                else:
                    if stack_list:
                        stack_list.insert(0, (data_node, node_id))
                        num_failures += 1
    except Exception, x:
            raise
    finally:
        os.chdir(working_dir)
    return True
