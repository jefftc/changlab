from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import geolib
        from genomicode import jmath

        # Input should be a GEOSeriesMatrixFile.
        filename = in_data.identifier

        convert_NA = user_options.get("set_NA_to")

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


    def name_outfile(self, antecedents, user_options):
        return "sample_metadata.tsv"

