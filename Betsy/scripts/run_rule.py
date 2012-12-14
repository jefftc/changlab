#! /usr/bin/env python
#run_rule.py
import os
import sys
import argparse
from Betsy import module_utils,config
from Betsy import rule_engine

def main():
    parser = argparse.ArgumentParser(description='Run the engine')
    parser.add_argument('--plin', dest='plin', default=None, action='append',
                        type=str, help='prolog command for pipeline start')
    parser.add_argument('--plout', dest='plout', default=None, type=str,
                        help='prolog command for querying')
    parser.add_argument('--id', dest='id', default=None, action='append',
                        type=str, help='sepecify the file or database id')
    parser.add_argument('--dry_run', dest='dry_run', action='store_const',
                        const=True, default=False,
                        help='show the test case procedure')
    parser.add_argument('--network', dest='network', type=str,
                        default=None,
                        help='generate the new network file')
    args = parser.parse_args()
    assert len(args.id) == len(args.plin), (
        'the length of id and plin should be the equal')
    pl_inputs = args.plin
    pl_output = args.plout
    identifiers = args.id
    for i in range(len(identifiers)):
        if os.path.exists(identifiers[i]):
            identifiers[i] = os.path.realpath(identifiers[i])
    objects = rule_engine.plstring2dataobject(pl_inputs, identifiers)
    print 'Start generate pipelines'
    pipelines = rule_engine.make_pipelines(pl_output, pl_inputs)
    assert pipelines, 'no pipelines can be generated'
    print '%d pipelines has been generated' %len(pipelines)
    if args.dry_run:
        print len(pipelines)
        for pipeline in pipelines:
            for analysis in pipeline:
                print 'module', analysis.name, '\r'
                print 'parameters', analysis.parameters
            print '------------------------'
        print len(pipelines)
    else:
        print 'Start running pipelines'
        k = 1
        for pipeline in pipelines:
            print  'pipeline' + str(k) + ':', '\r'
            rule_engine.run_pipeline(pipeline, objects)
            print '\r'
            k = k + 1
        print 'All pipelines has completed successfully'
    if args.network:
        network_file = os.path.join(os.getcwd(), args.network)
        original_network = config.NETWORKFILE
        for pipeline in pipelines:
            analysis_list = [analysis.name for analysis in pipeline]
            module_utils.high_light_path(
                original_network, analysis_list, network_file)
            original_network = network_file
        print ('new network file has been generated in %s'
               % os.path.realpath(network_file))


if __name__ == '__main__':
    main()
