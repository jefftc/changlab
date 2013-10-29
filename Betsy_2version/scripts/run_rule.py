#! /usr/bin/env python
#run_rule.py
import os
import sys
import argparse
from Betsy import rule_engine_bie
from Betsy import bie
from Betsy import rulebase


def main():
    parser = argparse.ArgumentParser(description='Run the engine')
    parser.add_argument('--in_data', dest='in_data', default=None, action='append',
                        type=str, help='input data object for network')
    parser.add_argument('--goal_datatype', dest='goal_datatype', default=None, type=str,
                        help='goal_datatype')
    parser.add_argument('--goal_attributes', dest='goal_attr', default=None, type=str,
                        help='goal_datatype')
    parser.add_argument('--dry_run', dest='dry_run', action='store_const',
                        const=True, default=False,
                        help='show the test case procedure')
    parser.add_argument('--network', dest='network',type=str, default=None,
                        help='generate the new network file')
    args = parser.parse_args()
    in_data = args.in_data
    goal_datatype = args.goal_datatype
    goal_attributes = args.goal_attr
    in_data = [eval(i) for i in in_data]
    goal_datatype = eval(goal_datatype)
    goal_attributes = eval(goal_attributes)
    print 'Generating network...'
    network = bie.backchain(rulebase.all_modules, goal_datatype, goal_attributes)
    network = bie.optimize_network(network)
    network = bie.prune_network_by_start(network, in_data)
    assert network, ('No pipeline has been generated, '
                       'please check your command.')
    if args.network:
        bie._plot_network_gv(args.network, network)
    if args.dry_run:
        bie._print_network(network)
    else:
        rule_engine_bie.run_pipeline(network,in_data)
        print 'All pipelines have completed successfully.'
     


if __name__ == '__main__':
    main()

