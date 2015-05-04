#add_crossplatform_probeid.py

import os
import shutil
import subprocess
import arrayio
from Betsy import module_utils
from genomicode import jmath, Matrix, arrayplatformlib, config
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    DATA = arrayio.read(in_data.identifier)
    chipname = arrayplatformlib.identify_platform_of_matrix(DATA)
    platform = user_options['platform_name']
    assert arrayplatformlib.get_bm_attribute(platform), (
        'the desire platform %s is not recognized by Betsy' % platform
    )
    if chipname == platform:
        shutil.copyfile(in_data.identifier, outfile)
    else:
        Annot_path = config.annotate_matrix
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN, 'cannot find the %s' % Annot_path
        command = ['python', Annot_BIN, in_data.identifier, "--platform",
                   platform]
        f = file(outfile, 'w')
        try:
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=f,
                                       stderr=subprocess.PIPE)
        finally:
            f.close()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
    assert module_utils.exists_nz(outfile), (
        'the output file %s for add_crossplatform_probeid fails' % outfile
    )
    out_node = bie3.Data(rulebase._SignalFile_Annotate, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.get_identifier(network, module_id, pool,
                                            user_attributes)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'signal_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
