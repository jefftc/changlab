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
    os.chdir(args.outpath)
    try:
        fullcommand = 'result<-run.analysis(gp.client,' + command + ')'
        R('library(GenePattern)')
        R('gp.client <- gp.login(servername, username, password)')
        R(fullcommand)
        R('download.directory <- job.result.get.job.number(result)')
        R('download.directory <- as.character(download.directory)')
        R('job.result.download.files(result, download.directory)')
        download_directory = os.path.realpath(R('download.directory')[0])
        print download_directory
    finally:
        os.chdir(cwd)

if __name__=='__main__':
    main()
