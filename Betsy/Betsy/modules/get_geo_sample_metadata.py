from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import geolib
        from genomicode import jmath


        convert_NA = user_options.get("set_NA_to")

        # Input should be a GEOSeriesMatrixFile.
        idata_in = antecedents
        filename = idata_in.identifier

        # Get the sample data and write it out.
        matrix = geolib._extract_sm_sample_meta(filename)
        matrix = geolib._clean_sm_sample_meta(matrix)
        matrix = geolib._prettify_sm_sample_meta(matrix)
        matrix = jmath.transpose(matrix)
           # each column is an annotation
        if convert_NA != "NA":
            for i in range(1, len(matrix)):
                for j in range(len(matrix[i])):
                    # Do case sensitive?
                    if matrix[i][j] == "NA":
                        matrix[i][j] = convert_NA
        
        outhandle = open(outfile, 'w')
        for x in matrix:
            print >>outhandle, "\t".join(map(str, x))
        
        outhandle.close()



    # Not sure exactly what this does.  Looks like it makes a unique name
    # for the output of this module.
    def name_outfile(self, antecedents, user_options):
        # BUG: what if there are multiple GPLIDs?
        #original_file = module_utils.get_inputid(user_options['GSEID'])
        #filename = "download_geo_seriesmatrix_%s" % original_file
        filename = "sample_metadata.tsv"
        return filename

