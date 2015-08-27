from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        # Given a GEOID and GPLID, get the series matrix file.
        from genomicode import geolib
    
        GSEID = user_options['GSEID']
        GPLID = user_options.get("GPLID")
        assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
        assert not GPLID or GPLID.startswith('GPL'), \
               'GPLID %s is not correct' % GPLID
    
        outhandle = open(outfile, 'w')
        geolib.download_seriesmatrix_file(outhandle, GSEID, GPLID)
        outhandle.close()
        #if not os.path.exists(outfile):
        #    os.mkdir(outfile)
        #matrix_files = get_seriesmatrix_file(GSEID, GPLID)
        #for matrix_file in matrix_files:
        #    newmatrix_filename = os.path.split(matrix_file)[-1]
        #    shutil.copyfile(matrix_file, os.path.join(outfile, newmatrix_filename))
        #assert module_utils.exists_nz(outfile), (
        #    'the output file %s for download_geo_dseriesmatrix fails' % outfile
        #)
    


        # Not sure exactly what this does.  Looks like it makes a unique name
        # for the output of this module.
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils

        # BUG: what if there are multiple GPLIDs?
        original_file = module_utils.get_inputid(user_options['GSEID'])
        filename = "download_geo_seriesmatrix_%s" % original_file
        return filename

    
