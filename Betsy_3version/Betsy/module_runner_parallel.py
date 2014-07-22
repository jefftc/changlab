#module_runner_parallel.py

import os
import subprocess
import random
import string
import sys
import time
import pickle
import tempfile
import config
import module_utils as module_utils
import bie3 as bie3
sys.modules['module_utils'] = module_utils
sys.modules['bie3'] = bie3
run_module_path = config.RUN_MODULE_PATH
queue_path = config.QUEUE_PATH

class ModuleJob:
    def __init__(self, module_jobname, outdata,
                 module_name, start_time, datfiles):
        self.module_jobname = module_jobname
        self.outdata = outdata
        self.module_name = module_name
        self.start_time = start_time
        self.input_data = datfiles



def run_module(network_dat, module_id, pool_dat, user_inputs_dat,
               pipeline_sequence_dat, user, job_name, clean_up=True):
    assert os.path.exists(network_dat)
    assert os.path.exists(pool_dat)
    assert os.path.exists(user_inputs_dat)
    assert os.path.exists(pipeline_sequence_dat)
    module_jobname = ''.join(random.choice(string.ascii_uppercase+string.digits)
                      for x in range(6))
    outfile = module_jobname+'.dat'
    job_command = [ 'python',run_module_path,
                       '--network',network_dat,
                       '--module_id',str(module_id),'--pool',pool_dat,
                       '--user_inputs',user_inputs_dat,'--pipeline_sequence',
                       pipeline_sequence_dat,'--user',user,'--out',outfile]
    if job_name:
        job_command.extend(['--job_name',job_name,])
    if clean_up:
        job_command.append('--clean_up')
    job_command = ' '.join(job_command)
    command = ['python',queue_path,'-j',module_jobname,job_command, '--user', user]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    message = process.communicate()
    error_message = message[1]
    if error_message:
            raise ValueError(error_message)
    handle = open(network_dat,'r')
    network = pickle.load(handle)
    handle.close()
    job = ModuleJob(module_jobname,outfile,network.nodes[module_id].name,
                    time.time(),[network_dat,pool_dat,user_inputs_dat,
                                 pipeline_sequence_dat])
    return job
    
        

def is_module_running(Job):
    command = ['python', queue_path, '--list']
    process = subprocess.Popen(command,shell=False,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    x,error_message = process.communicate()
    if error_message:
        raise ValueError(error_message)
    status = None
    for jobs in x.split('\n'):
        words = jobs.split('\t')
        jobid = words[0]
        jobname = words[1]
        if jobname == Job.module_jobname:
            status = words[2]
            break
    return status



def get_module_results(job,clean=False):
    outfile = job.outdata
    assert os.path.exists(outfile)
    f = open(outfile,'rb')
    out_data = pickle.load(f)
    f.close()
    if clean:
        os.remove(outfile)
    print '[' + time.strftime('%l:%M%p') + '].' + job.module_name
    return out_data


def get_run_time(job):
    current = time.time()
    elapse = current - job.start_time               
    return elapse

