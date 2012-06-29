#make_diffgenes_report.py
import protocol_utils

def run(outfiles,parameters,pipeline):
    protocol_utils.get_result_folder('differential_expressed_gene_analysis',outfiles,
                                     parameters,pipeline,'diffgenes_report')
