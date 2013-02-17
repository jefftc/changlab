#!/usr/bin/env python

import argparse
from genomicode import jmath, config


def main():
    parser = argparse.ArgumentParser(description='run the gene pattern module')
    
    parser.add_argument('--parameters', dest='parameters',
                        action='append', default=None,
                        type=str, help='key:value')
    parser.add_argument('-o', dest='outpath', type=str,
                        help='specify the outpath', default='.')
    parser.add_argument("module_name", nargs=1)
    args = parser.parse_args()
    module_name = args.module_name[0]
    
    parameters = dict()
    if args.parameters:
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
    params.append("'%s'" % module_name)
    for (key, value) in parameters.iteritems():
        params.append("%s='%s'" % (key, value))
    params_str = ", ".join(params)
    #print params_str
    R("result <- run.analysis(%s)" % params_str)
    jmath.R_equals(args.outpath, 'outpath')
    R('job.result.download.files(result, outpath)')


if __name__=='__main__':
    main()
