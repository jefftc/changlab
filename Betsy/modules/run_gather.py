#run_gather.py
import os
import module_utils
import urllib
import urllib2

def run(parameters,objects,pipeline):
    """run GATHER"""
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    kwargs={'cmd':"report",
        'tax_id':'9606',
        'network':0,
        'homologs':1,
        'annot_type':'gene_ontology'}
    network = {'yes_network':1,'no_network':0}
    homologs = {'yes_homologs':1,'no_homologs':0}
    if 'network' in parameters.keys():
        kwargs['network'] = network[parameters['network']]
    if 'homologs' in parameters.keys():
        kwargs['homologs'] = homologs[parameters['homologs']]
    if  'annot_type' in parameters.keys():
        kwargs['annot_type'] = parameters['annot_type']
    f = file(single_object.identifier,'r')
    text = f.read()
    f.close()
    if ',' in text:
        raise ValueError('the gene list file cannot contain ","')
    words = text.split()
    kwargs['gene_box'] = ','.join(words)
    url = 'http://gather.genome.duke.edu/'
    data = urllib.urlencode(kwargs)
    req = urllib2.Request(url,data)
    result_file = urllib2.urlopen(req)
    result_text=result_file.read()
    fout = file(outfile,'w')
    fout.write(result_text)
    fout.close()
    assert module_utils.exists_nz(outfile),(
        'the outfile for run_gather %s does not exist'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
        
def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'gather_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'gene_list_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for run_gather does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'gather_file',parameters,objects,single_object)
    return new_objects
