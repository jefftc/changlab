#!/usr/bin/env python

import argparse
from genomicode import jmath, config
import os

def main():
    parser = argparse.ArgumentParser(description='run the gene pattern module')
    
    parser.add_argument(
        '--parameters', default=[], action='append', help='key:value')
    parser.add_argument(
        '-o', dest='outpath', default=".",
        help='Directory to the save the results.')
    parser.add_argument("module_name", nargs=1)
    parser.add_argument(
        "--id_and_version", 
        help='specify the lsid and version in id:verison format')
    args = parser.parse_args()
    module_name = args.module_name[0]
    
    parameters = dict()
    for i in args.parameters:
        assert ':' in i, 'parameters should be in key:value format'
        key, value = i.split(':', 1)
        assert ':' not in value, 'parameters should be in key:value format'
        parameters[key] = value

    # given the module_name and the module parameters in dict, call
    # module in Genepattern
    R = jmath.start_R()
    jmath.R_equals(config.gp_user, 'username')
    jmath.R_equals(config.gp_passwd, 'password')
    jmath.R_equals(config.gp_server, 'servername')
    R('library(GenePattern)')
    R('gp.client <- gp.login(servername, username, password)')
    
    params = []
    params.append("gp.client")
    if args.id_and_version:
        params.append("'urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:%s'" % args.id_and_version)
    else:
        params.append("'%s'" % module_name)
    for (key, value) in parameters.iteritems():
        params.append("%s='%s'" % (key, value))
    params_str = ", ".join(params)
    x = "result <- run.analysis(%s)" % params_str=
    R(x)

    # Download the files to outpath.
    jmath.R_equals(args.outpath, 'outpath')
    R('job.result.download.files(result, outpath)')
    assert os.path.exists(args.outpath), \
           "Missing output directory for: %s" % module_name

    # Look for "stderr.txt".
    result_files = os.listdir(args.outpath)
    assert 'stderr.txt' not in result_files, (
        "Run failed.  GenePattern generated an error:\n%s" %
        file(os.path.join(args.outpath,'stderr.txt')).read()
        )

if __name__=='__main__':
    main()
