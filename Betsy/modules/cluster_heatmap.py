#cluster_heatmap.py
import protocol_utils

def run(outfiles,parameters,pipeline):
    protocol_utils.get_result_folder('make_heatmap',outfiles,
                                     parameters,pipeline,'heatmap_report')
