from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import urllib
        import urllib2
        from genomicode import filelib

        in_data = antecedents
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
        assert filelib.exists_nz(outfile), (
            'the outfile for run_gather %s does not exist' % outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        return 'gather_' + original_file + '.txt'
