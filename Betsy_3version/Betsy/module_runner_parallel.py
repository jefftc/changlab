#module_runner_parallel.py
import pickle
import os
import subprocess
import random
import string
import module_utils as module_utils
import bie3 as bie3
import sys
import time
sys.modules['module_utils'] = module_utils
sys.modules['bie3'] = bie3

run_module_path = '/home/xchen/chencode/Betsy_3version/Betsy/run_module.py'
queue_path = '/home/xchen/chencode/queue/queue.py'

class ModuleJob:
    def __init__(self, module_jobname, outdata, module_name,start_time):
        self.module_jobname = module_jobname
        self.outdata = outdata
        self.module_name = module_name
        self.start_time = start_time

def run_module(network, module_id, pool, user_inputs,
               pipeline_sequence, user, job_name, clean_up=True):
    f = file('network.dat','wb')
    pickle.dump(network, f)
    f.close()
    f = file('pool.dat', 'wb')
    pickle.dump(pool,f)
    f.close()
    f = file('user_inputs.dat','wb')
    pickle.dump(user_inputs,f)
    f.close()
    f = file('pipeline_sequence.dat','wb')
    pickle.dump(pipeline_sequence,f)
    f.close()
    module_jobname = ''.join(random.choice(string.ascii_uppercase+string.digits)
                          for x in range(6))
    outfile = module_jobname+'.dat'
    job_command = [ 'python',run_module_path,
                       '--network','network.dat',
                       '--module_id',str(module_id),'--pool','pool.dat',
                       '--user_inputs','user_inputs.dat','--pipeline_sequence',
                       'pipeline_sequence.dat','--user',user,'--out',outfile]
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
    job = ModuleJob(module_jobname,outfile,network.nodes[module_id].name,time.time())
    return job


def is_module_running(Job):
    command = ['python', queue_path, '--list']
    handle = file('list.txt','w')
    process = subprocess.Popen(command,shell=False,
                                stdout=handle,
                                stderr=subprocess.PIPE)
    handle.close()
    error_message = process.communicate()[1]
    if error_message:
            raise ValueError(error_message)
    status = None
    with open('list.txt') as FileObj:
        for jobs in FileObj:
            words = jobs.split('\t')
            jobid = words[0]
            jobname = words[1]
            if jobname == Job.module_jobname:
                status = words[2]
                break
    return status
##    if status == 'running':
##        return True
##    return False  #[running,completed,killed,pending]


def get_module_results(job):
    outfile = job.module_jobname+'.dat'
    assert os.path.exists(outfile)
    f=open(outfile,'rb')
    out_data = pickle.load(f)
    f.close()
    os.remove(outfile)
    print '[' + time.strftime('%l:%M%p') + '].' + job.module_name
    return out_data


def get_run_time(job):
    current = time.time()
    elapse = current - job.start_time               
    return elapse
