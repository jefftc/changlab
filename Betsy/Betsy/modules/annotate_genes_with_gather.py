#annotate_genes_with_gather.py
import os
from Betsy import module_utils
import urllib
import urllib2
from Betsy import bie3, rulebase


def run(network, antecedents, out_attributes, user_options, num_cores):
    """run GATHER"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    kwargs = {
        'cmd': "report",
        'tax_id': '9606',
        'network': 0,
        'homologs': 0,
        'annot_type': 'gene_ontology'
    }
    code_dict = {'yes': 1, 'no': 0}
    kwargs['network'] = code_dict[out_attributes['network']]
    kwargs['homologs'] = code_dict[out_attributes['homologs']]
    kwargs['annot_type'] = out_attributes['annot_type']
    f = file(in_data.identifier, 'r')
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
    out_node = bie3.Data(rulebase.GatherFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(antecedents.identifier)
    filename = 'gather_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)
