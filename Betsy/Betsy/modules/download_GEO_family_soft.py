from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """given a GEOID  get the family soft file"""
        from Betsy import module_utils
        #in_data = antecedents
        GSEID = user_options['GSEID']

        assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
        download_series_family(GSEID, 300, open(outfile, 'w'))
        assert module_utils.exists_nz(outfile), (
            'the output file %s for download_GEO_family_soft fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(user_options['GSEID'])
        filename = original_file + '_family.soft.gz'
        return filename

    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        return user_options['GSEID']


def download_series_family(GSEID, DELAY, outhandle):
    import urllib
    from genomicode import timer

    timer.wait(DELAY)
    url = ("ftp://ftp.ncbi.nih.gov/pub/geo/" +
           "DATA/SOFT/by_series/%s/" % GSEID + "%s_family.soft.gz" % GSEID)
    handle = urllib.urlopen(url)
    for line in handle:
        print >>outhandle, line,
