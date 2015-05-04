#normalize_samples_with_bfrm.py
import os
import subprocess
import arrayio
from genomicode import config
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    bfrm_path = config.bfrmnorm
    bfrm_BIN = module_utils.which(bfrm_path)
    assert bfrm_BIN, 'cannot find the %s' % bfrm_path
    num_factor = 1
    #num_factor = 10
    if 'num_factors' in user_options.keys():
        num_factor = int(user_options['num_factors'])
        assert num_factor >= 1, 'the num_factor should be >=1'
        # What is single_object?
        #M = arrayio.read(single_object.identifier)
        M = arrayio.read(in_data.identifier)
        col_num = M.ncol()
        assert num_factor <= col_num, (
            'the num_factor should be less than %d' % col_num
        )
    
    tmp = 'tmp_dir'
    command = ['python', bfrm_BIN, in_data.identifier, '-f', str(num_factor),
               '-o', tmp]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    
    assert module_utils.exists_nz(tmp), (
        'the output dir %s for bfrm_normalize fails' % tmp
    )
    assert module_utils.exists_nz(os.path.join(tmp, 'normalized.gct')), (
        'the output gct file for bfrm_normalize fails'
    )
    out = os.path.join(tmp, 'normalized.gct')
    M = arrayio.read(out)
    M_new = arrayio.convert(M, to_format=arrayio.pcl_format)
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M_new, f)
    f.close()
    out_node = bie3.Data(rulebase._SignalFile_Merge, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_bfrm_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def get_out_attributes(antecedents, out_attributes):
    new_parameters = antecedents.data.attributes.copy()
    new_parameters['bfrm_norm'] = 'yes'
    return new_parameters


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
