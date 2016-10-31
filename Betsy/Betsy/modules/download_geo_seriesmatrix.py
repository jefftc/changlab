from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        # Given a GEOID and GPLID, get the series matrix file.
        from genomicode import geolib

        metadata = {}
        
        GSEID = user_options['GSEID']
        GPLID = user_options.get("GPLID")
        assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
        assert not GPLID or GPLID.startswith('GPL'), \
               'GPLID %s is not correct' % GPLID
        # Don't need to save user_options.
        #metadata["GSEID"] = GSEID
        #if GPLID:
        #    metadata["GPLID"] = GPLID
    
        outhandle = open(outfile, 'w')
        geolib.download_seriesmatrix_file(outhandle, GSEID, GPLID)
        outhandle.close()
        filelib.assert_exists_nz(outfile)
        #metadata["filesize"] = filelib.filesize(outfile)
        #if not os.path.exists(outfile):
        #    os.mkdir(outfile)
        #matrix_files = get_seriesmatrix_file(GSEID, GPLID)
        #for matrix_file in matrix_files:
        #    newmatrix_filename = os.path.split(matrix_file)[-1]
        #    shutil.copyfile(matrix_file, os.path.join(outfile, newmatrix_filename))
        #assert filelib.exists_nz(outfile), (
        #    'the output file %s for download_geo_dseriesmatrix fails' % outfile
        #)
        return metadata

    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        ## BUG: what if there are multiple GPLIDs?
        #original_file = module_utils.get_inputid(user_options['GSEID'])
        #filename = "download_geo_seriesmatrix_%s" % original_file
        #return filename
        return "seriesmatrix.txt"

    
