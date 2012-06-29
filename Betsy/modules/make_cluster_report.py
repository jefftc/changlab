#make_cluster_report.py
import protocol_utils

def run(outfiles,parameters,pipeline):
    protocol_utils.get_result_folder('cluster_genes',outfiles,
                                     parameters,pipeline,'cluster_report')
