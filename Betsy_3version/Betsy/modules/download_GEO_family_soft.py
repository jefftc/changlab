#download_GEO_family_soft.py
import os
from Betsy import module_utils,bie3,rulebase
import sys

def limit_query_rate(DELAY):
    global GEO_TIMER
    import time
    time.sleep(DELAY)

def download_series_family(**params):
    outhandle = sys.stdout
    for name in params:
        exec "%s = params['%s']" % (name, name)
    import urllib
    limit_query_rate(DELAY)
    url = ("ftp://ftp.ncbi.nih.gov/pub/geo/" +
           "DATA/SOFT/by_series/%s/" % GSEID +
           "%s_family.soft.gz" % GSEID)
    handle = urllib.urlopen(url)
    for line in handle:
        print >>outhandle, line,
        
def run(data_node, parameters, user_input, network):
    """given a GEOID  get the family soft file"""
    outfile = name_outfile(data_node,user_input)
    GSEID = user_input['GSEID']
    
    assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
    download_series_family(
    GSEID=GSEID, DELAY=300, outhandle=open(outfile, 'w'))
    assert module_utils.exists_nz(outfile), (
        'the output file %s for download_GEO_family_soft fails' % outfile)
    out_node = bie3.Data(rulebase.GEOfamily, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(data_node,user_input):
    original_file = module_utils.get_inputid(user_input['GSEID'])
    filename = original_file+'_family.soft.gz'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(data_node,pipeline,parameters,user_input):
    identifier = user_input['GSEID']
    return identifier

def get_out_attributes(parameters,data_node):
    return parameters

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes)
    return data_node



