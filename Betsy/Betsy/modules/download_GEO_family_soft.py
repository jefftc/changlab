#download_GEO_family_soft.py
import os
from Betsy import module_utils, bie3, rulebase


def limit_query_rate(DELAY):
    global GEO_TIMER
    import time
    time.sleep(DELAY)


def download_series_family(GSEID, DELAY, outhandle):
    import urllib
    
    limit_query_rate(DELAY)
    url = ("ftp://ftp.ncbi.nih.gov/pub/geo/" +
           "DATA/SOFT/by_series/%s/" % GSEID + "%s_family.soft.gz" % GSEID)
    handle = urllib.urlopen(url)
    for line in handle:
        print >> outhandle, line,


def run(network, antecedents, out_attributes, user_options, num_cores):
    """given a GEOID  get the family soft file"""
    in_data = antecedents
    outfile = name_outfile(in_data, user_options)
    GSEID = user_options['GSEID']

    assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
    download_series_family(GSEID, 300, open(outfile, 'w'))
    assert module_utils.exists_nz(outfile), (
        'the output file %s for download_GEO_family_soft fails' % outfile
    )
    out_node = bie3.Data(rulebase.GEOFamilySoftFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    original_file = module_utils.get_inputid(user_options['GSEID'])
    filename = original_file + '_family.soft.gz'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    identifier = user_options['GSEID']
    return identifier


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node = module_utils.find_antecedents(network, module_id, user_attributes,
                                            pool)
    return data_node
