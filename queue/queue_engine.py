#queue_engine.py
import subprocess
import tempfile
import os
import datetime
"""

job_submit(job_name,command_line)

job_kill(job_number)

queue_status()

job_name2job_number(job_name_of_job_number)
"""

Log_file = '/home/xchen/chencode/queue/log.txt'
processing_info='/home/xchen/chencode/queue/processing_info/'

def hash_command(time, command_line):
    from hashlib import md5
    hashstring=time+command_line
    hash = md5()
    hash.update(hashstring)
    hash_result = hash.hexdigest()
    return hash_result


def job_submit(job_name, command_line,user):
    f = file(Log_file,'r')
    text = f.readlines()
    f.close()
    jobnames = []
    for line in text[1:]:
        (job_no, jobname_log, status_log,
        started,ended,submitted,command,log_user) = line.split('\t')
        jobnames.append(jobname_log)
    assert job_name not in jobnames,(
        'A job name %s already exisits, please give a different job name'
        %job_name)
    submitted = datetime.datetime.now()
    submitted = submitted.strftime('%m/%d%l:%M%p')
    info_file = os.path.join(processing_info,
                             hash_command(submitted,command_line) + '.txt')
    try:
        filename = tempfile.mktemp()
        f = file(filename,'w')
        f.write(command_line)
        f.write(' &> '+ info_file)
        f.close()
        command = ['batch','-f', filename]
        process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        job_information = process.communicate()[1]
        job_number = job_information.split()[1]
        add2log((job_number, job_name,
                 'pending',' ',' ', submitted, str(command_line),user))
        return job_number
    finally:
        os.remove(filename)

        
def add2log(x):
    job_number,job_name,status,started,ended,submitted, command, user = x
    f = file(Log_file,'a')
    f.write('\t'.join([str(job_number),job_name,
                       status,started,ended,submitted,command,user])+'\n')
    f.close()


def job_kill(job_number=None,kill_all=False):
    #check the current jobs
    update_dict = check_status()
    #kill the job
    if kill_all:
        command = ['atrm']+ [str(job) for job in update_dict]
    elif job_number:
        command = ['atrm', str(job_number)]
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    message = process.communicate()
    #update the log file
    if job_number and str(job_number) in update_dict:
        update_log(job_number,'killed')
    if kill_all:
        for job in update_dict:
            update_log(str(job),'killed')

    
def check_status():
    update_dict = {}
    command = ['atq']
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    job_information = process.communicate()[0]
    jobs = job_information.split('\n')
    jobs = [job for job in jobs if job]
    for line in jobs:
        job_num,job_info = line.split('\t')
        job_status = 'pending'
        if '=' in job_info:
            job_status = 'running'
        update_dict[job_num] = job_status
    return update_dict


def update_log(job_number=None, status=None):
    f = file(Log_file,'r')
    text = f.readlines()
    f.close()
    log_dict = dict()
    text = [i.strip() for i in text if i.strip()]
    for line in text[1:]:
        (job_no, jobname_log, status_log,
         started, ended,submitted,command,user) = line.split('\t')
        log_dict[job_no] = [jobname_log,
                            status_log,started, ended,submitted,command,user]
    update_dict = check_status()
    if job_number and status and str(job_number) in log_dict:
        log_dict[str(job_number)] = [log_dict[str(job_number)][0],
                                     status]+log_dict[str(job_number)][2:]
        
    for key in log_dict:
        if key in update_dict:
            log_dict[key] = [log_dict[key][0],update_dict[key]]+log_dict[key][2:]
        elif not log_dict[key][1] == 'killed':
            log_dict[key] = [log_dict[key][0],'completed']+log_dict[key][2:]
    f = file(Log_file,'w')
    f.write('\t'.join(['Jobnumber','Jobname',
                       'Status','Started','Ended','Submitted','Command','User'])+'\n')
    #update the start time and end time, write to log.txt
    keys = log_dict.keys()
    keys.sort()
    check_time = datetime.datetime.now()
    check_time = check_time.strftime('%m/%d%l:%M%p')
    for key in keys:
        job_status = log_dict[key][1]
        job_start = log_dict[key][2]
        job_end = log_dict[key][3]
        if job_status == 'running' and job_start == ' ':
            log_dict[key][2]= check_time
        if job_status in ['completed','killed'] and job_end == ' ':
            log_dict[key][3] = check_time
        f.write('\t'.join([key]+log_dict[key])+'\n')
    f.close()


def queue_status():
    """read the log.txt and show the job status"""
    update_log()
    f=file(Log_file,'r')
    text=f.readlines()
    f.close()
    for line in text:
        print line.strip()


def job_name2job_number(job_name_or_job_number):
    f = file(Log_file,'r')
    text = f.readlines()
    f.close()
    log_dict = dict()
    text = [i.strip() for i in text if i.strip()]
    job_numbers = []
    for line in text[1:]:
        (job_no, jobname_log, status_log,
         started,ended,submitted, command) = line.split('\t')
        log_dict[jobname_log] = job_no
        job_numbers.append(job_no)
    #if input is a job number
    if job_name_or_job_number in job_numbers:
        return job_name_or_job_number
    #if input is a job name
    job_number = log_dict[job_name_or_job_number]
    return job_number

