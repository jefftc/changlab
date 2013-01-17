#!/usr/bin/env python

import os
import argparse
import sys
from genomicode import jmath,config
import zipfile


def main():

    parser = argparse.ArgumentParser(description = 'run the gene pattern module')
    
    parser.add_argument('--parameters',dest = 'parameters',
                        action = 'append',default = None,
                        type = str,help='key:value')
    parser.add_argument('-o',dest = 'outpath',type = str,
                        help = 'specify the outpath',default='.')
    parser.add_argument(dest = 'module_name',type = str,
                        help = 'specify the module_name',default = None)
    args = parser.parse_args()
    parameters = dict()
    if args.parameters:
        for i in args.parameters:
            assert ':' in i, 'parameters should be in key:value format'
            key,value = i.split(':')
            assert ':' not in value,'parameters should be in key:value format'
            parameters[key] = value
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)
    """given the module_name and the module parameters
       in dict, call module in Genepatttern"""
    R = jmath.start_R()
    username = config.gp_user
    password = config.gp_passwd 
    servername = config.gp_server
    jmath.R_equals(password, 'password')
    jmath.R_equals(servername, 'servername')
    jmath.R_equals(username, 'username')
    command = "\'" + args.module_name + "\'"
    for key in parameters.keys():
        command = command + ',' + key + '=' + '\"' + parameters[key] + '\"'
    cwd = os.getcwd()
    try:
        os.chdir(args.outpath)
        fullcommand = 'result<-run.analysis(gp.client,' + command + ')'
        R('library(GenePattern)')
        R('gp.client <- gp.login(servername, username, password)')
        R(fullcommand)
        jmath.R_equals(args.outpath,'outpath')
        R('job.result.download.files(result, outpath)')
    finally:
        os.chdir(cwd)

if __name__=='__main__':
    main()
