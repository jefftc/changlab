#annotate_genes_with_gather.py
import os
from Betsy import module_utils
import urllib
import urllib2
from Betsy import bie3, rulebase


def run(data_node, parameters, user_input, network, num_cores):
    """run GATHER"""
    outfile = name_outfile(data_node, user_input)
    kwargs = {
        'cmd': "report",
        'tax_id': '9606',
        'network': 0,
        'homologs': 0,
        'annot_type': 'gene_ontology'
    }
    code_dict = {'yes': 1, 'no': 0}
    kwargs['network'] = code_dict[parameters['network']]
    kwargs['homologs'] = code_dict[parameters['homologs']]
    kwargs['annot_type'] = parameters['annot_type']
    f = file(data_node.identifier, 'r')
    text = f.read()
    f.close()
    if ',' in text:
        raise ValueError('the gene list file cannot contain ","')
    words = text.split()
    kwargs['gene_box'] = ','.join(words)
    url = 'http://gather.genome.duke.edu/'
    data = urllib.urlencode(kwargs)
    req = urllib2.Request(url, data)
    result_file = urllib2.urlopen(req)
    result_text = result_file.read()
    fout = file(outfile, 'w')
    fout.write(result_text)
    fout.close()
    assert module_utils.exists_nz(outfile), (
        'the outfile for run_gather %s does not exist' % outfile
    )
    out_node = bie3.Data(rulebase.GatherFile, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes)
    return data_node


def name_outfile(data_node, user_input):
    original_file = module_utils.get_inputid(data_node.identifier)
    filename = 'gather_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(parameters, data_nodes):
    return parameters


def make_unique_hash(data_node, pipeline, parameters, user_input):
    identifier = data_node.identifier
    return module_utils.make_unique_hash(identifier, pipeline, parameters,
                                         user_input)
