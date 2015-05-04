#convert_signal_to_pcl.py
import os
from Betsy import module_utils
import shutil
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """convert signal file to pcl format"""
    import arrayio
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    f = file(outfile, 'w')
    M = arrayio.choose_format(in_data.identifier)
    if M.__name__[8:-7] == 'pcl':
        shutil.copyfile(in_data.identifier, outfile)
        f.close()
    else:
        M = arrayio.read(in_data.identifier)
        M_c = arrayio.convert(M, to_format=arrayio.pcl_format)
        arrayio.pcl_format.write(M_c, f)
        f.close()
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for convert_signal_to_pcl fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Normalize, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_' + original_file + '.pcl'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)

    return data_node
