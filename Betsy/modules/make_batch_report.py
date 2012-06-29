#make_batch_report.py
import protocol_utils

def run(outfiles,parameters,pipeline):
    protocol_utils.get_result_folder('batch_effect_remove',outfiles,
                                     parameters,pipeline,'batch_effect_remove_report')
