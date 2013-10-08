#annotate_genes_with_gather.py
import os
#from Betsy
import module_utils
import urllib
import urllib2
from time import strftime,localtime
import bie
import rulebase

def run(data_node,parameters):
    """run GATHER"""
    outfile = name_outfile(data_node)
    kwargs = {'cmd': "report",
              'tax_id': '9606',
              'network': 0,
              'homologs': 0,
              'annot_type': 'gene_ontology'}
    code_dict = {'yes': 1, 'no': 0}
    kwargs['network'] = code_dict[parameters['network']]
    kwargs['homologs'] = code_dict[parameters['homologs']]
    kwargs['annot_type'] = parameters['annot_type']
    f = file(data_node.attributes['filename'], 'r')
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
        'the outfile for run_gather %s does not exist' % outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.GatherFile,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                    data_nodes)
    return data_node

def name_outfile(data_node):
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'gather_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,x):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):

    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)



