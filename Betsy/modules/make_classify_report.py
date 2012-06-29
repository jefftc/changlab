#make_classify_report.py
import protocol_utils

def run(outfiles,parameters,pipeline):
    protocol_utils.get_result_folder('classification',outfiles,
                                     parameters,pipeline,'classify_report')
