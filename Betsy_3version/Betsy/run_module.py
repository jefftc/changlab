#run_module.py
import pickle
import os
import getpass
import argparse
import config
import time
from  module_utils import DataObject as DataObject
import rulebase
import bie3


def main():
    parser = argparse.ArgumentParser(description='run the module')
    parser.add_argument('--network', dest='network', type=str,
                        help='network file',
                        default=None)
    parser.add_argument('--module_id',dest='module_id', type=str,
                        help='module_id')
    parser.add_argument('--pool',dest='pool', type=str,
                        default=None, help='pool')
    parser.add_argument('--user_inputs', dest='user_inputs',type=str,
                        default=None,
                        help='user_input')
    parser.add_argument('--pipeline_sequence', dest='pipeline_sequence', default=None,
                        type=str,help='pipeline sequence')
    parser.add_argument('--user', dest='user', default=False,
                        type=str,help='user')
    parser.add_argument('--job_name', dest='job_name', default=None,
                        type=str,help='job name')
    parser.add_argument('--clean_up', dest='clean_up', default=False,
                        action='store_true',help='clean up the job files')
    parser.add_argument('--out', dest='out', default=None,
                        type=str,help='outfile')
    args = parser.parse_args()
    x = convert_input_to_objects(args.network,args.pool,args.user_inputs,
                                 args.pipeline_sequence)
    network,pool,user_inputs, pipeline_sequence, = x
    out_node = run_module(network, int(args.module_id), pool, user_inputs,
                          pipeline_sequence, args.user, args.job_name, args.clean_up)
    f = file(args.out,'wb')
    pickle.dump(out_node,f)
    f.close()


def convert_input_to_objects(network_file, pool_file,
                             user_inputs_file, pipeline_sequence_file):
    assert os.path.exists(network_file)
    assert os.path.exists(pool_file)
    assert os.path.exists(user_inputs_file)
    assert os.path.exists(pipeline_sequence_file)
    f = file(network_file,'r')
    network = pickle.load(f)
    f.close()
    f = file(pool_file,'r')
    pool = pickle.load(f)
    f.close()
    f = file(user_inputs_file,'r')
    user_inputs = pickle.load(f)
    f.close()
    f = file(pipeline_sequence_file,'r')
    pipeline_sequence = pickle.load(f)
    f.close()
    return network, pool, user_inputs, pipeline_sequence


def run_module(network, module_id, pool, user_inputs, pipeline_sequence,
               user=getpass.getuser(), job_name='', clean_up=True):
    current_dir = os.getcwd()
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
    if out_attributes is None:
        return []
    print '[' + time.strftime('%l:%M%p') + '].' + module_name
    working_dir = os.path.join(output_path, make_module_wd_name(
        network, module, module_id, sub_user_input, pipeline_sequence,
        pool, out_attributes))
    # make name of outfile
    data_node = module.find_antecedents(network, module_id,
                                        pool, out_attributes)
    outfile = os.path.split(module.name_outfile(data_node, sub_user_input))[-1]
    temp_dir = ''
    temp_outfile = ''
    if not os.path.exists(os.path.join(working_dir, 'stdout.txt')):
        #if no result has been generated, create temp folder and run analysis
        temp_dir = tempfile.mkdtemp()
        try:
            os.chdir(temp_dir) 
            starttime = strftime(module_utils.FMT, localtime())
            out_node = module.run(data_node, out_attributes,
                                  sub_user_input, network)
            if not out_node:
                return False
            module_utils.write_Betsy_parameters_file(out_node.data.attributes,
                                                     data_node, sub_user_input,
                                                     pipeline_sequence,
                                                     starttime, user, job_name)
            f = file(os.path.join(temp_dir, 'stdout.txt'), 'w')
            f.write('This module runs successully.')
            f.close()
        except:
            logging.exception('Got exception on %s .run()' % module_name)
            raise
        finally:
            os.chdir(current_dir)
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
            if (module_utils.exists_nz(os.path.join(temp_dir,outfile)) and not
                os.path.exists(os.path.join(working_dir, 'stdout.txt'))):
                copy_result_folder(working_dir, temp_dir, outfile)
        finally:
            if clean_up and os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
    out_nodes = get_out_node(working_dir,
        network, module_id, module, out_attributes, pool, sub_user_input)
    assert out_nodes, 'module %s fails' % module_node.name
    return out_nodes


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


def get_out_node(working_dir, network, module_id, module_file, parameters,
                pool, user_input):
    current_dir = os.getcwd()
    try:
        os.chdir(working_dir)
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
    finally:
        os.chdir(current_dir)

        
def copy_result_folder(working_dir, temp_dir, temp_outfile):
    """copy result files in temp folder to result folder,
      if fails, delete result folder"""
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    try:
        shutil.copyfile(os.path.join(temp_dir, 'Betsy_parameters.txt'),
                        os.path.join(working_dir, 'Betsy_parameters.txt'))
        if os.path.isdir(os.path.join(temp_dir,temp_outfile)):
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


if __name__ == '__main__':
    main()
