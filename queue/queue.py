#!/usr/bin/env python
#queue.py

import queue_engine 
import argparse
import string
import random
import getpass

def main():
    parser = argparse.ArgumentParser(description='run the queue system')
    parser.add_argument('job_command', type=str, nargs='*',
                        help='submit the job command to the queue',
                        default=None)
    parser.add_argument('-j',dest='jobname', type=str, nargs='+',
                        help='submit the jobname and command to the queue')
    parser.add_argument('--user',dest='user', type=str,
                        default=getpass.getuser(),
                        help='name of the user who submit the job')
    parser.add_argument('--list', '-l',dest='status',
                        action='store_true', default=False,
                        help='show the job status')
    parser.add_argument('--kill', dest='kill_jobname', default=None,
                        type=str,help='kill the job')
    parser.add_argument('--killall', dest='kill_all', default=False,
                        action='store_true',help='kill the job')
    args = parser.parse_args()
    if args.job_command:
        jobname = ''.join(random.choice(string.ascii_uppercase+string.digits)
                          for x in range(6))
        command_line = ' '.join(args.job_command)
        job_number = queue_engine.job_submit(jobname,command_line,args.user)
    if args.jobname:
        assert len(args.jobname)>1,'jobname and command are required'
        job_name = args.jobname[0]
        command_line = ' '.join(args.jobname[1:])
        job_number = queue_engine.job_submit(job_name,command_line,args.user)
    if args.status:
        queue_engine.queue_status()
    if args.kill_jobname:
        job_number = queue_engine.job_name2job_number(args.kill_jobname)
        queue_engine.job_kill(job_number)
    if args.kill_all:
        queue_engine.job_kill(None,True)
        
if __name__=='__main__':
    main()

