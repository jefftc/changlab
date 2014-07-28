#rule_engine_bie3_paralleling.py
import config
import bie3
import logging
import sys
import module_utils
from  module_utils import DataObject as DataObject
import os
import time
import getpass
from time import strftime, localtime
import module_runner_parallel as module_runner
from Betsy import module_utils 
import tempfile
import pickle


"""
Functions:
save_objects
choose_next_module
test_require_data
run_pipeline
"""
class TempFileManager:
    def __init__(self):
        self.files = []
    def make(self):
        fd,filename = tempfile.mkstemp()
        self.files.append(filename)
        return fd, filename
    def flush(self):
        for temp in self.files:
            if os.path.exists(temp):
                os.remove(temp)
        self.files = []
    def __del__(self):
        self.flush()
        
def save_objects(network,pool,user_inputs,pipeline_sequence, tempfile_maganer):
    fd1, network_dat = tempfile_maganer.make()
    fd2, pool_dat = tempfile_maganer.make()
    fd3, user_inputs_dat = tempfile_maganer.make()
    fd4, pipeline_sequence_dat = tempfile_maganer.make()
    f1 = file(network_dat,'wb')
    pickle.dump(network, f1)
    f1.close()
    os.close(fd1)
    f2 = file(pool_dat, 'wb')
    pickle.dump(pool,f2)
    f2.close()
    os.close(fd2)
    f3= file(user_inputs_dat,'wb')
    pickle.dump(user_inputs,f3)
    f3.close()
    os.close(fd3)
    f4 = file(pipeline_sequence_dat,'wb')
    pickle.dump(pipeline_sequence,f4)
    f4.close()
    os.close(fd4)
    x = (network_dat, pool_dat, user_inputs_dat, pipeline_sequence_dat)
    return x


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



def run_pipeline(network, in_objects, user_inputs,
                 user=getpass.getuser(), job_name=''):
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    LOG_FILENAME = os.path.join(output_path, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
    pipeline_sequence = []
    stack_list = []
    stack_jobs = []
    pool = {}
    next_node = None
    tempfile_manager = TempFileManager()
    for data_object in in_objects:
        data_node = data_object.data
        start_node_ids = bie3._find_start_nodes(network, data_node)
        assert start_node_ids,'%s cannot matched to network' % data_node.datatype.name
        assert len(start_node_ids)==1,'data match to more than one node'
        start_node_id = start_node_ids[0]
        stack_list.append((data_object, start_node_id))
        pool[start_node_id] = data_object
    num_failures = 0
    try:
        while stack_list:
            assert num_failures < len(stack_list)
            x = stack_list.pop()
            if isinstance(x, tuple):
                data_object, node_id = x
                if isinstance(data_object, DataObject):
                    pool[node_id] = data_object
                    if node_id == 0:
                        return None
                    next_modules = choose_next_module(network, node_id, pool)
                    stack_list.extend(next_modules)
                elif isinstance(data_object, bie3.Module):
                    module_id = node_id
                    test_required = test_require_data(network, module_id, pool)
                    if not test_required:
                        stack_list.insert(0, (data_object, node_id))
                        num_failures += 1
                        continue
                    pipeline_sequence.append(data_object.name)
                    x = save_objects(network, pool, user_inputs, pipeline_sequence,tempfile_manager)
                    network_dat, pool_dat, user_inputs_dat, pipeline_sequence_dat = x
                    job = module_runner.run_module(network_dat, module_id, pool_dat,
                                           user_inputs_dat, pipeline_sequence_dat,
                                           user, job_name)
                    stack_list.append(job)
                    
                else:
                    raise Exception
            elif isinstance(x, module_runner.ModuleJob):
                job = x
                if job.outdata not in tempfile_manager.files:
                    tempfile_manager.files.append(job.outdata)
                if module_runner.get_run_time(job)>100:
                    for in_dat in job.input_data:
                            assert os.path.exists(in_dat)
                            os.remove(in_dat)
                    raise ValueError('job %s is running too long'%job.module_name)
                
                if module_runner.is_module_running(job) in ['pending','running']:
                    stack_list.insert(0,job)
                    continue
                elif module_runner.is_module_running(job)=='completed':
                    for in_dat in job.input_data:
                        assert os.path.exists(in_dat)
                        os.remove(in_dat)
                    out_nodes = module_runner.get_module_results(job,clean=True)
                    for x in out_nodes:
                        next_node, next_id = x
                        if next_id == 0:
                            break
                        elif next_id is None:
                            raise ValueError('cannot match the output node')
                        stack_list.append((next_node, next_id))
                else:
                    raise ValueError('running status is not valid')
            else:
                 raise Exception
        if next_node and module_utils.exists_nz(next_node.identifier):
            print ('[' + time.strftime('%l:%M%p') +
                   '] Completed successfully and ' +
                   'generated a file:')
            print  next_node.identifier + '\r'
            print '\r'
            sys.stdout.flush()
            return next_node.identifier
        else:
            print 'This pipeline has completed unsuccessfully'
            raise ValueError(
                'there is no output for this pipeline')
        return None
    except:
        raise 
    finally:
        tempfile_manager.flush()
